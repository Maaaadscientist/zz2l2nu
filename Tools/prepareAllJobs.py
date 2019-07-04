#!/usr/bin/env python

from __future__ import division, print_function
import argparse
from collections import defaultdict
import copy
from glob import glob
import math
import os
import re
import shutil
import subprocess
import sys

import yaml

from dataset import Dataset


def parse_command_line():
    """Parse the command line parser"""
    parser = argparse.ArgumentParser(description='Launch baobab nutple production.')

    parser.add_argument('--listDataset', action='store', required=True,
                        help='Specifies the file containing the list of task to submit. The file must contain one task specification per line. The specification consist of three space-separated value: catalog of input files, pruner selection, pruner subselection. The catalog can be an eos path (/store/...).')
    parser.add_argument('--suffix', action='store', default=None,
                        help='suffix that will be added to the output directory')
    parser.add_argument('--harvest', action='store_true', default=None,
                        help='harvest the root files from the last submission')
    parser.add_argument('--config', default='2016.yaml',
                        help='Master configuration for the analysis.')
    parser.add_argument('--isPhotonDatadriven', action='store_true', default=None,
                        help='Launch HZZ with Instr. MET DY estimated from photon')
    parser.add_argument('--doInstrMETAnalysis', action='store_true', default=None,
                        help='Launch InstrMETAnalysis')
    parser.add_argument('--doNRBAnalysis', action='store_true', default=None,
                        help='Launch NRB Analysis')
    parser.add_argument('--express', action='store_true', default=None,
        help='Launch on the express queue (better if: very fast jobs "<10min" in a small amount "<16"). This queue has a wallTime of 32min and can only take 8 jobs per user')
    parser.add_argument('--localCopy', action='store_true', default=None,
        help='Copy the ROOT files locally on the node to improve speed.')
    parser.add_argument('--syst', action='store', default=None,
        help='Specify the systematic on which you need to run. If you dont specify _up or _down, both will be run at the same time. Use "all" to run on all the systematics defined in the systList.txt file.')

    return parser.parse_args()


def hadd(sources, output_path, overwrite=False):
    """Merge ROOT files with hadd.

    Arguments:
        sources:      Sequence of paths of source files to be merged.
            Glob-like masks are supported.
        output_path:  Name for output file with results of the merge.
        overwrite:    Specifies whether the output file should be
            overwritten if it already exists.

    Return value:
        None.

    If there is a single source file, copy it as output instead of
    calling hadd.
    """

    expanded_sources = []

    for mask in sources:
        expanded_sources += glob(mask)

    if len(expanded_sources) == 0:
        raise RuntimeError('Empty list of source files.')

    if os.path.exists(output_path):
        os.remove(output_path)

    if len(expanded_sources) == 1:
        # No need to call hadd, just copy the file
        shutil.copyfile(expanded_sources[0], output_path)
    else:
        subprocess.check_call(['hadd', output_path] + expanded_sources)


def parse_datasets_file(path):
    """Parse file with a list of dataset definition files.

    Paths to dataset definition files (DDF) are given in a text file,
    one per line.  A common directory with respect to which these paths
    are resolved can be specified optionally.  Empty lines and lines
    that only contain comments are skipped.

    Arguments:
        path:  Path to a file with a list of dataset definition files.

    Return value:
        List of constructed datasets.
    """

    # Extract paths to dataset definition files
    ddfs = []

    blank_regex = re.compile(r'^\s*$')
    comment_regex = re.compile(r'^\s*#')
    directory_regex = re.compile(r'^\s*catalogPath\s*=\s*(.+)\s*$')

    datasets_file = open(path, 'r')
    directory = ''

    for line in datasets_file:
        if blank_regex.match(line) or comment_regex.search(line):
            continue

        match = directory_regex.match(line)

        if match:
            directory = match.group(1)
        else:
            ddfs.append(line.strip())

    datasets_file.close()
    ddfs = [os.path.join(directory, ddf) for ddf in ddfs]


    # Read stem dataset definitions if available
    config_dir = os.path.join(os.environ['HZZ2L2NU_BASE'], 'config/')

    with open(config_dir + args.config) as f:
        config = yaml.safe_load(f)

    if 'dataset_stems' in config:
        with open(config_dir + config['dataset_stems']) as f:
            stem_list = yaml.safe_load(f)
            stems = {stem['name']: stem for stem in stem_list}
    else:
        stems = {}


    datasets = [Dataset(ddf, stems) for ddf in ddfs]
    return datasets


def parse_syst_file():
    global base_path
    try:
      systFile = open(base_path+"/systList.txt")
    except KeyError:
      sys.stderr.write("cannot open syst file")
    systLines = systFile.readlines()
    theSystDict={}
    for aLine in systLines:
      if aLine.startswith("//"): continue
      if aLine.startswith("#"): continue
      systKey=(aLine.split(" ")[0])
      thisLineAsList = aLine.split()
      theSystDict[systKey] = []
      for i in thisLineAsList:
        if i==thisLineAsList[0]: continue
        if i.startswith("//"): break
        if i.startswith("#"): break
        theSystDict[systKey].append(i)
    return theSystDict

def find_syst_in_file(currentSyst):
    global base_path
    try:
      systFile = open(base_path+"/systList.txt")
    except KeyError:
      sys.stderr.write("cannot open syst file")
    systLines = systFile.readlines()
    thisSystList=[]
    for aLine in systLines:
      if aLine.startswith("//"): continue
      if aLine.startswith("#"): continue
      if not currentSyst in aLine: continue
      thisLineAsList = aLine.split()
      for i in thisLineAsList:
        if i==thisLineAsList[0]: continue
        if i.startswith("//"): break
        if i.startswith("#"): break
        thisSystList.append(i)
    return thisSystList
    

def prepare_local_copy(dataset, skip_files, max_files, ddf_save_path):
    """Prepare local copying of input ROOT files.

    Construct shell commands to copy selected part of the given dataset
    to the working directory of the runnig job.  Write a new dataset
    definition file that contains only the selected input files.

    Arguments:
        dataset:  Datasets from which to copy files.
        skip_files:  Number of input files to skip.
        max_files:   Maximal number of files to process.  A value of -1
            means all remaining files.
        ddf_save_path:  Where to save the dataset definition file with
            local copies of the selected input files.

    Return value:
        List of strings with shell commands to copy input files.
    """

    if max_files < 0:
        max_files = len(dataset.files)

    script_commands = []
    local_files = []

    for path in dataset.files[skip_files:skip_files + max_files]:
        script_commands.append('dccp {} .'.format(path))
        local_files.append(os.path.basename(path))

    dataset_clone = copy.copy(dataset)
    dataset_clone.files = local_files
    dataset_clone.save(ddf_save_path)

    return script_commands


def extract_list_of_systs(syst):
    dictOfSysts = {}
    if not syst or syst=="no":
      dictOfSysts[None] = [""]
    elif ("_up" in syst) or ("_down" in syst):
      dictOfSysts[syst] = find_syst_in_file(syst)
    elif syst=="all":
      dictOfSysts=parse_syst_file()
      dictOfSysts[None] = [""]
    else:
      dictOfSysts[syst+"_up"] = find_syst_in_file(syst+"_up")
      dictOfSysts[syst+"_down"] = find_syst_in_file(syst+"_down")
      dictOfSysts[None] = [""] # Run also on the nominal shape
    return dictOfSysts


def prepare_job_script(dataset, syst, job_id=0, skip_files=0, max_files=-1):
    """Create a script to be run on a batch system

    The script performs set-up and executes runHZZanalysis.  It is saved
    in the jobs directory.

    Arguments:
        dataset:  Dataset to be processed.
        syst:     Label of requested systematic variation.
        job_id:   Number for the current job among all jobs in the given
            dataset.
        skip_files:  Number of input files from the dataset to skip.
        max_files:   Maximal number of input files to process in this
            job.  A value of -1 means all remaining files.

    Return value:
        None.
    """

    global base_path
    global thisSubmissionDirectory
    global outputDirectory
    global jobsDirectory

    if syst:
        job_name = '{}_{}_{}'.format(dataset.name, syst, job_id)
    else:
        job_name = '{}_{}'.format(dataset.name, job_id)

    script_commands = [
      'export INITDIR={}'.format(base_path),
      'cd $INITDIR',
      '. ./env.sh',
      'cd -',
      'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi',
      'hostname',
      'date'
    ]

    if args.localCopy:
        ddf_path = '{}/scripts/ddf_{}{}.yaml'.format(
            jobsDirectory, outputPrefixName, job_name
        )
        script_commands += prepare_local_copy(
            dataset, skip_files, max_files, ddf_path
        )
    else:
        ddf_path = dataset.path

    # Construct options for runHZZanalysis program
    options = [
        '--config={}'.format(args.config),
        '--catalog={}'.format(ddf_path),
        '--output={}{}.root'.format(outputPrefixName, job_name),
        '--skip-files={}'.format(
            0 if args.localCopy else skip_files
        ),
        '--max-files={}'.format(max_files), '--max-events=-1'
    ]

    if doInstrMETAnalysis:
        options.append('--analysis=InstrMET')
    elif doNRBAnalysis:
        options.append('--analysis=NRB')
    else:
        options.append('--analysis=Main')

    if isPhotonDatadriven:
        options.append('--dd-photon')

    if syst:
        options.append('--syst={}'.format(syst))

    if args.syst != 'all':
        options.append('--all-control-plots')

    # Use debug-level verbosity when running on the batch system since
    # one is expected to consult logs only in case of a problem
    options.extend(['-v', '2'])

    run_application_command = ' '.join(['runHZZanalysis'] + options)
    script_commands.append('echo ' + run_application_command)
    script_commands.append(run_application_command)

    script_commands.append(
        'cp {}{}.root {}'.format(outputPrefixName, job_name, outputDirectory)
    )


    script_path = '{}/scripts/runOnBatch_{}{}.sh'.format(
        jobsDirectory, outputPrefixName, job_name
    )

    with open(script_path, 'w') as f:
        for command in script_commands:
            f.write(command)
            f.write('\n')

    with open(
        '{}/sendJobs_{}.cmd'.format(thisSubmissionDirectory, args.suffix), 'a'
    ) as jobs_file:
        jobs_file.write(
            'qsub {} -j oe -o {}/logs/ {}\n'.format(
                doExpress, jobsDirectory, script_path
            )
        )


def prepare_jobs(dataset, syst):
    """Construct jobs for given dataset.

    Decide on the job splitting and delegate the construction of
    individual jobs to the dedicated function.

    Arguments:
        dataset:  Dataset for which to construct jobs.
        syst:     Label of requested systematic variation.

    Return value:
        None.
    """

    print(
        'Preparing scripts for dataset '
        '\033[1;33m{}\033[0;m'.format(dataset.name)
    )

    # For some systematic variations all files within a dataset have to
    # be processed within a single job
    if syst and ('pdf' in syst or 'QCDscale' in syst):
        prepare_job_script(dataset, syst)
    else:
        # Choose the number of files per job based on the average number
        # of events per file, if relevant information is available
        num_events = dataset.parameters.get('num_selected_events', None)

        if num_events is not None:
            events_per_job = 1000000
            job_splitting = max(int(round(num_events / events_per_job)), 1)
        else:
            job_splitting = 25

        num_jobs = int(math.ceil(len(dataset.files) / job_splitting))

        for job_id in range(num_jobs):
            prepare_job_script(
                dataset, syst, job_id,
                skip_files=job_id * job_splitting,
                max_files=job_splitting
            )


def runHarvesting():
    global thisSubmissionDirectory
    global outputDirectory

    datasets = parse_datasets_file(args.listDataset)
    merge_dir = os.path.join(thisSubmissionDirectory, 'MERGED')

    if not os.path.isdir(merge_dir):
        print('\033[1;34m Will create directory "{}"\033[0;m'.format(
            merge_dir
        ))
        os.mkdir(merge_dir)

    listForFinalPlots = defaultdict(list)
    listForFinalPlots_data = []
    dictOfSysts = extract_list_of_systs(args.syst)

    for currentSyst in dictOfSysts:
        dataSamplesList = []
        dataForThisSyst = None

        if not currentSyst:
            systString = ''
        else:
            systString = '_' + currentSyst

        for dataset in datasets:
            ddf_filename = os.path.basename(dataset.path)
            harvestForThisSyst = None

            for key in dictOfSysts[currentSyst]:
                if key in ddf_filename:
                    harvestForThisSyst = True

            if not harvestForThisSyst:
                continue

            if not 'Bonzais' in ddf_filename:
                continue

            print('\033[1;32m merging {}{}\033[0;m'.format(
                dataset.name, systString
            ))
            full_name = outputPrefixName + dataset.name + systString
            hadd(
                ['{}/{}_[0-9]*.root'.format(outputDirectory, full_name)],
                '{}/{}.root'.format(merge_dir, full_name), overwrite=True
            )

            if not dataset.is_sim:
                dataForThisSyst = True
                dataSamplesList.append('{}/{}.root'.format(
                    merge_dir, full_name
                ))
            else:
                listForFinalPlots[dataset.name].append('{}/{}.root'.format(
                    merge_dir, full_name
                ))

        if dataForThisSyst:
            listForFinalPlots_data.append('{}/{}Data{}.root'.format(
                merge_dir, outputPrefixName, systString
            ))

        if currentSyst:
            print(
                '\033[1;32m merging all Data (Single* and Double*) together '
                'for {}\033[0;m'.format(currentSyst)
            )
        else:
            print(
                '\033[1;32m merging all Data (Single* and Double*) together '
                'for nominal shapes\033[0;m'
            )

        if dataForThisSyst:
            hadd(
                dataSamplesList,
                '{}/{}Data{}.root'.format(
                    merge_dir, outputPrefixName, systString
                ),
                overwrite=True
            )

    if args.syst and not args.syst == 'no':
        print('\033[1;32m producing final ROOT files with all shapes\033[0;m')
        for key in listForFinalPlots:
            hadd(
                listForFinalPlots[key],
                '{}/{}{}_final.root'.format(merge_dir, outputPrefixName, key),
                overwrite=True
            )
        hadd(
            listForFinalPlots_data,
            '{}/{}Data_final.root'.format(merge_dir, outputPrefixName),
            overwrite=True
        )


def main():
    global args
    global base_path
    global thisSubmissionDirectory
    global outputDirectory
    global jobsDirectory
    global doInstrMETAnalysis
    global isPhotonDatadriven
    global outputPrefixName
    global doNRBAnalysis
    global doExpress
    #create the directories if needed
    base_path=os.path.expandvars('$HZZ2L2NU_BASE')
    if not os.path.isdir(base_path+"/OUTPUTS"):
        print("\033[1;31m OUTPUTS directory does not exist: will create it \033[0;m")
        os.mkdir(base_path+"/OUTPUTS")

    args = parse_command_line()
    
    if type(args.suffix) != type("txt"):
        thisSubmissionDirectory=base_path+"/OUTPUTS/Test"
    else:
        thisSubmissionDirectory=base_path+"/OUTPUTS/"+args.suffix
    
    outputDirectory=thisSubmissionDirectory+"/OUTPUTS"
    jobsDirectory=thisSubmissionDirectory+"/JOBS"
    
    if not os.path.isdir(thisSubmissionDirectory):
        print("\033[1;31m will create the directory "+thisSubmissionDirectory+"\033[0;m")
        os.mkdir(thisSubmissionDirectory)
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)
    if not os.path.isdir(jobsDirectory):
        os.mkdir(jobsDirectory)
        os.mkdir(jobsDirectory+"/scripts")
        os.mkdir(jobsDirectory+"/logs")

    #options
    outputPrefixName="outputHZZ_"
    if args.doInstrMETAnalysis:
        print("Preparing InstrMET analysis...\n")
        doInstrMETAnalysis = 1
        outputPrefixName="outputInstrMET_"
    else:
        doInstrMETAnalysis = 0

    if args.isPhotonDatadriven:
        print("Datadriven estimation of the Instr. MET option found...\n")
        isPhotonDatadriven = 1
        outputPrefixName="outputPhotonDatadriven_"
    else:
        isPhotonDatadriven = 0


    if args.doNRBAnalysis:
        print("Praparing Non-resonant Bkg. Analysis...\n")
        doNRBAnalysis = 1
        outputPrefixName="outputNRB_"
    else:
        doNRBAnalysis = 0
    if args.harvest:
        print("will harvest")
        runHarvesting()
        return

    if args.express:
        print("Will be launched on the express queue (NB: only do this for small and fast jobs)\n")
        doExpress = " -q express -l walltime=00:30:00 "
    else:
        print("WallTime is set to 20h. If you need more, please update the script. If you need to send only a small number of very short jobs, please consider using the express queue (--express)\n")
        doExpress = " -l walltime=20:00:00 "


    dictOfSysts = extract_list_of_systs(args.syst)

    for thisSyst in dictOfSysts:
      for dataset in parse_datasets_file(args.listDataset):
        runOnThisCatalog = None
        for key in dictOfSysts[thisSyst]:
          if key in os.path.basename(dataset.path): runOnThisCatalog = True
        if runOnThisCatalog: prepare_jobs(dataset, thisSyst)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt, e:
        print("\nBye!")
        pass
