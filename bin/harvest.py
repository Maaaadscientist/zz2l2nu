#!/usr/bin/env python

from __future__ import division, print_function
import argparse
from collections import defaultdict
from glob import glob
import itertools
import os
import shutil
import subprocess

import yaml

from hzz import Dataset, SystDatasetSelector, parse_datasets_file


def parse_command_line():
    """Parse the command line parser"""
    parser = argparse.ArgumentParser(description='Launch baobab nutple production.')

    parser.add_argument('--listDataset', action='store', required=True,
                        help='Specifies the file containing the list of task to submit. The file must contain one task specification per line. The specification consist of three space-separated value: catalog of input files, pruner selection, pruner subselection. The catalog can be an eos path (/store/...).')
    parser.add_argument('--suffix', action='store', default=None,
                        help='suffix that will be added to the output directory')
    parser.add_argument('--config', default='2016.yaml',
                        help='Master configuration for the analysis.')
    parser.add_argument('--isPhotonDatadriven', action='store_true', default=None,
                        help='Launch HZZ with Instr. MET DY estimated from photon')
    parser.add_argument('--doInstrMETAnalysis', action='store_true', default=None,
                        help='Launch InstrMETAnalysis')
    parser.add_argument('--doNRBAnalysis', action='store_true', default=None,
                        help='Launch NRB Analysis')
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


def runHarvesting():
    global thisSubmissionDirectory
    global outputDirectory

    merge_dir = os.path.join(thisSubmissionDirectory, 'MERGED')

    if not os.path.isdir(merge_dir):
        print('\033[1;34m Will create directory "{}"\033[0;m'.format(
            merge_dir
        ))
        os.mkdir(merge_dir)

    datasets = parse_datasets_file(args.listDataset, args.config)


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
    global outputPrefixName
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
    

    outputPrefixName="outputHZZ_"

    if args.doInstrMETAnalysis:
        outputPrefixName="outputInstrMET_"

    if args.isPhotonDatadriven:
        outputPrefixName="outputPhotonDatadriven_"

    if args.doNRBAnalysis:
        outputPrefixName="outputNRB_"

    runHarvesting()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt, e:
        print("\nBye!")
        pass

