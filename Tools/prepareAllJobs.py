#!/usr/bin/env python

from __future__ import division, print_function
import argparse
from collections import defaultdict
import copy
from glob import glob
import math
import itertools
import os
import re
import shutil
import subprocess
import sys

import yaml

from dataset import Dataset
from util import SystDatasetSelector


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

    args = parser.parse_args()

    if not args.syst or args.syst == 'no':
        args.syst = ''

    return args


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
        expanded_sources_cur_mask = glob(mask)

        if not expanded_sources_cur_mask:
            raise RuntimeError(
                'No files found matching mask "{}".'.format(mask)
            )

        expanded_sources += expanded_sources_cur_mask

    if os.path.exists(output_path):
        os.remove(output_path)

    if len(expanded_sources) == 1:
        # No need to call hadd, just copy the file
        shutil.copyfile(expanded_sources[0], output_path)
    else:
        subprocess.check_output(['hadd', output_path] + expanded_sources)


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

    merge_dir = os.path.join(thisSubmissionDirectory, 'MERGED')

    if not os.path.isdir(merge_dir):
        print('\033[1;34m Will create directory "{}"\033[0;m'.format(
            merge_dir
        ))
        os.mkdir(merge_dir)

    datasets = parse_datasets_file(args.listDataset)


    # Merge all data files
    data_masks = []

    for dataset in datasets:
        if dataset.is_sim:
            continue

        data_masks.append('{}/{}{}_[0-9]*.root'.format(
            outputDirectory, outputPrefixName, dataset.name
        ))

    print('\033[1;32m Merging all data files...\033[0;m')
    hadd(
        data_masks,
        '{}/{}Data{}.root'.format(
            merge_dir, outputPrefixName,
            '_final' if args.syst else ''
        ),
        overwrite=True
    )


    # Merge simulation, including different systematic variations
    merge_paths = defaultdict(list)  # All variations for each dataset
    dataset_selector = SystDatasetSelector(
        os.path.join(base_path, 'config/syst.yaml')
    )

    for variation, dataset in itertools.chain(
        [('', d) for d in datasets if d.is_sim],  # Nominal variation
        dataset_selector(datasets, args.syst)
    ):
        print('\033[1;32m Merging dataset "{}" variation "{}"...'
              '\033[0;m'.format(dataset.name, variation))
        syst_postfix = '_' + variation if variation else ''
        full_name = outputPrefixName + dataset.name + syst_postfix
        source_mask = '{}/{}{}{}_[0-9]*.root'.format(
            outputDirectory, outputPrefixName, dataset.name, syst_postfix
        )
        merge_path = '{}/{}{}{}.root'.format(
            merge_dir, outputPrefixName, dataset.name, syst_postfix
        )
        hadd([source_mask], merge_path, overwrite=True)
        merge_paths[dataset.name].append(merge_path)

    if args.syst:
        print(
            '\033[1;32m Merging all systematic variations...\033[0;m')

        for dataset_name, sources in merge_paths.items():
            hadd(
                sources, '{}/{}{}_final.root'.format(
                    merge_dir, outputPrefixName, dataset_name
                ),
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


    datasets = parse_datasets_file(args.listDataset)

    if not args.syst or args.syst == 'all':
        # Nominal configuration for systematic variations
        for dataset in datasets:
            prepare_jobs(dataset, '')

    if args.syst:
        # Systematic variations
        dataset_selector = SystDatasetSelector(
            os.path.join(base_path, 'config/syst.yaml')
        )

        for variation, dataset in dataset_selector(datasets, args.syst):
            prepare_jobs(dataset, variation)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt, e:
        print("\nBye!")
        pass
