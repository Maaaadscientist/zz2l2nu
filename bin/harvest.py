#!/usr/bin/env python

"""Merges results of individual jobs."""


import argparse
from collections import defaultdict
from glob import glob
import itertools
import os
import shutil
import subprocess

import yaml
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


class Grouper:
    """Implements grouping of datasets."""

    def __init__(self):
        config_path = os.path.join(
            os.environ['HZZ2L2NU_BASE'], 'config/dataset_groups.yaml')
        with open(config_path) as f:
            group_defs = yaml.safe_load(f)

        self.group_map = {}
        for group, dataset_names in group_defs.items():
            for name in dataset_names:
                if name in self.group_map:
                    raise RuntimeError(
                        f'Dataset "{name}" is included in multiple groups.')
                self.group_map[name] = group

    def __call__(self, dataset_name):
        """Return name of group to which given dataset belogs.

        If the dataset is not included in any group, return the
        dataset's name.
        """
        return self.group_map.get(dataset_name, dataset_name)


def harvest(datasets, source_dir, merge_dir, prefix='', syst='',
            tree_analysis=False):
    """Merge outputs produced for given datasets.

    Arguments:
        datasets:    Datasets to be processed.
        source_dir:  Directory with results of individual jobs, which
            are to be merged.
        merge_dir:   Directory in which merged files will be placed.
        prefix:      Prefix used in the names of source and merged ROOT
            files.
        syst:  Group of systematic variations to process.
        tree_analysis:  Indicates whether this is a tree-based or
            histogram-based analysis.

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

    if data_masks:
        print('\033[1;32mMerging all data files...\033[0;m')
        hadd(
            data_masks,
            '{}/{}Data{}.root'.format(
                merge_dir, prefix,
                '_final' if syst and not tree_analysis else ''
            ),
            overwrite=True
        )


    # Determine how files should be merged in simulation, while
    # respecting grouping of datasets.  Describe merging rules with a
    # mapping from group and variation names to lists of included
    # datasets and masks with input files.
    merge_rules = {}
    dataset_selector = SystDatasetSelector(
        os.path.join(os.environ['HZZ2L2NU_BASE'], 'config/syst.yaml')
    )
    grouper = Grouper()

    for variation, dataset in dataset_selector(
        datasets, syst, combine_weights=tree_analysis
    ):
        group = grouper(dataset.name)
        syst_postfix = '_' + variation if variation else ''
        source_mask = '{}/{}{}{}_[0-9]*.root'.format(
            source_dir, prefix, dataset.name, syst_postfix
        )
        if (group, variation) not in merge_rules:
            merge_rules[group, variation] = ([dataset.name], [source_mask])
        else:
            v = merge_rules[group, variation]
            v[0].append(dataset.name)
            v[1].append(source_mask)


    # Now the files according to the constructed rules
    for (group, variation), (dataset_names, sources) in merge_rules.items():
        print('\033[1;32mMerging datasets {} variation "{}"...'
              '\033[0;m'.format(
            ', '.join([f'"{name}"' for name in dataset_names]),
            variation
        ))
        merge_path = '{}/{}{}{}.root'.format(
            merge_dir, prefix, group, '_' + variation if variation else ''
        )
        hadd(sources, merge_path, overwrite=True)


    # Merge all variations within every group if this is not a
    # tree-based analysis
    if syst and not tree_analysis:
        merge_paths = defaultdict(list)
        for (group, _), (_, sources) in merge_rules:
            merge_paths[group].append(sources)

        print(
            '\033[1;32mMerging all systematic variations...\033[0;m')

        for group, sources in merge_paths.items():
            hadd(
                sources, '{}/{}{}_final.root'.format(
                    merge_dir, prefix, group
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
        '--syst', default='',
        help='Requested systematic variation or a group of them.'
    )
    arg_parser.add_argument(
        '--prefix', default='',
        help='Prefix for names of output files.'
    )
    arg_parser.add_argument(
        '--hist-analysis', action='store_true',
        help='Indicates that results of a histogram-based analysis are going '
        'to be harvested.'
    )
    args = arg_parser.parse_args()

    if args.syst == 'no':
        args.syst = ''


    source_dir = os.path.join(args.task_dir, 'output')
    merge_dir = os.path.join(args.task_dir, 'merged')

    if not os.path.exists(merge_dir):
        print('\033[1;34mWill create directory "{}"\033[0;m'.format(
            merge_dir
        ))
        os.mkdir(merge_dir)


    datasets = parse_datasets_file(args.datasets, args.config)
    harvest(datasets, source_dir, merge_dir, prefix=args.prefix,
            syst=args.syst, tree_analysis=not args.hist_analysis)

