#!/usr/bin/env python

"""Merges results of individual jobs."""

from __future__ import division, print_function
import argparse
from collections import defaultdict
from glob import glob
import itertools
import os
import shutil
import subprocess

from hzz import Dataset, SystDatasetSelector, parse_datasets_file


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


def harvest(datasets, source_dir, merge_dir, prefix='', syst=''):
    """Merge outputs produced for given datasets.

    Arguments:
        datasets:    Datasets to be processed.
        source_dir:  Directory with results of individual jobs, which
            are to be merged.
        merge_dir:   Directory in which merged files will be placed.
        syst:  Group of systematic variations to process.

    If syst == 'all', files with all systematic vairations for the same
    dataset will be merged together.
    """

    # Merge all data files
    data_masks = []

    for dataset in datasets:
        if dataset.is_sim:
            continue

        data_masks.append('{}/{}{}_[0-9]*.root'.format(
            source_dir, prefix, dataset.name
        ))

    print('\033[1;32m Merging all data files...\033[0;m')
    hadd(
        data_masks,
        '{}/{}Data{}.root'.format(
            merge_dir, prefix, '_final' if syst else ''
        ),
        overwrite=True
    )


    # Merge simulation, including different systematic variations
    merge_paths = defaultdict(list)  # All variations for each dataset
    dataset_selector = SystDatasetSelector(
        os.path.join(os.environ['HZZ2L2NU_BASE'], 'config/syst.yaml')
    )

    for variation, dataset in dataset_selector(datasets, syst):
        print('\033[1;32m Merging dataset "{}" variation "{}"...'
              '\033[0;m'.format(dataset.name, variation))
        syst_postfix = '_' + variation if variation else ''
        full_name = prefix + dataset.name + syst_postfix
        source_mask = '{}/{}{}{}_[0-9]*.root'.format(
            source_dir, prefix, dataset.name, syst_postfix
        )
        merge_path = '{}/{}{}{}.root'.format(
            merge_dir, prefix, dataset.name, syst_postfix
        )
        hadd([source_mask], merge_path, overwrite=True)
        merge_paths[dataset.name].append(merge_path)

    if syst:
        print(
            '\033[1;32m Merging all systematic variations...\033[0;m')

        for dataset_name, sources in merge_paths.items():
            hadd(
                sources, '{}/{}{}_final.root'.format(
                    merge_dir, prefix, dataset_name
                ),
                overwrite=True
            )


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument(
        'datasets', help='File with a list of datasets to process.'
    )
    arg_parser.add_argument(
        '-d', '--task-dir', default='task',
        help='Directory for scripts and results of this task.'
    )
    arg_parser.add_argument(
        '--config', default='2016.yaml',
        help='Master configuration for the analysis.'
    )
    arg_parser.add_argument(
        '-a', '--analysis', default='Main',
        help='Analysis that was run as given to runHZZanalysis.'
    )
    arg_parser.add_argument(
        '--dd-photon', action='store_true',
        help='Whether data-driven estimation for photon+jets as used.'
    )
    arg_parser.add_argument(
        '--syst', default='',
        help='Requested systematic variation or a group of them.'
    )
    args = arg_parser.parse_args()

    if args.syst == 'no':
        args.syst = ''
    

    source_dir = os.path.join(args.task_dir, 'output')
    merge_dir = os.path.join(args.task_dir, 'merged')

    if not os.path.exists(merge_dir):
        print('\033[1;34m Will create directory "{}"\033[0;m'.format(
            merge_dir
        ))
        os.mkdir(merge_dir)


    if args.analysis == 'Main':
        prefix = 'outputHZZ_'
    elif args.analysis == 'InstrMET':
        prefix = 'outputInstrMET_'
    elif args.analysis == 'NRB':
        prefix = 'outputNRB_'
    else:
        raise RuntimeError('Unrecognized analysis "{}".'.format(args.analysis))

    if args.dd_photon:
        prefix = 'outputPhotonDatadriven_'


    datasets = parse_datasets_file(args.datasets, args.config)
    harvest(datasets, source_dir, merge_dir, prefix=prefix, syst=args.syst)

