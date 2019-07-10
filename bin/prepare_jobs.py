#!/usr/bin/env python

"""Prepares PBS jobs for the analysis."""

from __future__ import division, print_function
import argparse
import copy
import math
import os

from hzz import Dataset, SystDatasetSelector, parse_datasets_file


class JobBuilder:
    def __init__(
        self, task_dir, config_path, analysis_options='', output_prefix='',
        local_copy=False
    ):
        """Initialize the builder.
        
        Arguments:
            task_dir:  Directory for the task.  Created if needed.
            config_path:  Path to master configuration for the analysis.
                Forwarded to runHZZanalysis.
            analysis_options:  List with additional options to be
                forwared to runHZZanalysis.
            output_prefix:  Prefix to be added to names of output files.
            local_copy:  Requests that input ROOT files are copied to
                the node on which the job is running.
        """

        # Make sure task_dir is an absolute path because it is used to
        # define the output path for jobs
        self.task_dir = os.path.abspath(task_dir)

        self.config_path = config_path
        self.analysis_options = analysis_options
        self.output_prefix = output_prefix
        self.local_copy = local_copy
        self.install_path = os.environ['HZZ2L2NU_BASE']

        self._create_directories()
        self.submit_commands = []


    def prepare_jobs(self, dataset, syst):
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
            self._prepare_job_script(dataset, syst)
        else:
            # Choose the number of files per job based on the average number
            # of events per file, if relevant information is available
            num_events = dataset.parameters.get('num_selected_events', None)

            if num_events is not None:
                events_per_job = 1000000

                if num_events > events_per_job:
                    job_splitting = int(math.ceil(
                        len(dataset.files) / (num_events / events_per_job)
                    ))
                else:
                    job_splitting = len(dataset.files)
            else:
                job_splitting = 25

            num_jobs = int(math.ceil(len(dataset.files) / job_splitting))

            for job_id in range(num_jobs):
                self._prepare_job_script(
                    dataset, syst, job_id,
                    skip_files=job_id * job_splitting,
                    max_files=job_splitting
                )


    def write_submit_script(self, name='send_jobs.sh'):
        """Write script to submit jobs.
        
        All jobs created so far by method prepare_jobs are included.
        The script is created in the task directory.
        """

        with open(os.path.join(self.task_dir, name), 'w') as f:
            for command in self.submit_commands:
                f.write(command)
                f.write('\n')


    def _create_directories(self):
        """Create directory structure for the task if needed."""

        if not os.path.exists(self.task_dir):
            print('\033[1;31m Will create task directory '
                  '"{}"\033[0;m'.format(self.task_dir))
            os.makedirs(self.task_dir)

        sub_dirs = ['output', 'jobs', 'jobs/scripts', 'jobs/logs']

        if self.local_copy:
            sub_dirs.append('jobs/ddfs')

        for sub_dir in sub_dirs:
            try:
                os.makedirs(os.path.join(self.task_dir, sub_dir))
            except OSError:
                pass


    def _prepare_local_copy(
        self, dataset, skip_files, max_files, ddf_save_path
    ):
        """Prepare local copying of input ROOT files.

        Construct shell commands to copy selected part of the given
        dataset to the working directory of the runnig job.  Write a new
        dataset definition file that contains only the selected input
        files.

        Arguments:
            dataset:  Datasets from which to copy files.
            skip_files:  Number of input files to skip.
            max_files:   Maximal number of files to process.  A value of
                -1 means all remaining files.
            ddf_save_path:  Where to save the dataset definition file
                with local copies of the selected input files.

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


    def _prepare_job_script(
        self, dataset, syst, job_id=0, skip_files=0, max_files=-1
    ):
        """Create script defining the PBS job.

        The script performs set-up and executes runHZZanalysis.  It is
        saved in subdirectory jobs/scripts.

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

        if syst:
            job_name = '{}_{}_{}'.format(dataset.name, syst, job_id)
        else:
            job_name = '{}_{}'.format(dataset.name, job_id)

        script_commands = [
          'export INITDIR={}'.format(self.install_path),
          'cd $INITDIR',
          '. ./env.sh',
          'cd -',
          'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi',
          'hostname',
          'date'
        ]

        if self.local_copy:
            ddf_path = '{}/jobs/ddfs/{}{}.yaml'.format(
                self.task_dir, self.output_prefix, job_name
            )
            script_commands += self._prepare_local_copy(
                dataset, skip_files, max_files, ddf_path
            )
        else:
            ddf_path = dataset.path

        # Construct options for runHZZanalysis program
        options = [
            '--config={}'.format(self.config_path),
            '--catalog={}'.format(ddf_path),
            '--output={}{}.root'.format(self.output_prefix, job_name),
            '--skip-files={}'.format(
                0 if self.local_copy else skip_files
            ),
            '--max-files={}'.format(max_files), '--max-events=-1'
        ]

        options += self.analysis_options

        if syst:
            options.append('--syst={}'.format(syst))

        # Use debug-level verbosity when running on the batch system
        # since one is expected to consult logs only in case of a
        # problem
        options.extend(['-v', '2'])

        run_application_command = ' '.join(['runHZZanalysis'] + options)
        script_commands.append('echo ' + run_application_command)
        script_commands.append(run_application_command)

        script_commands.append('cp {}{}.root {}/output'.format(
            self.output_prefix, job_name, self.task_dir
        ))


        script_path = '{}/jobs/scripts/runOnBatch_{}{}.sh'.format(
            self.task_dir, self.output_prefix, job_name
        )

        with open(script_path, 'w') as f:
            for command in script_commands:
                f.write(command)
                f.write('\n')

        self.submit_commands.append(
            'qsub -l walltime=20:00:00 -j oe -o {}/jobs/logs/ {}\n'.format(
                self.task_dir, script_path
            )
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
        help='Analysis to be run. Forwarded to runHZZanalysis.'
    )
    arg_parser.add_argument(
        '--dd-photon', action='store_true',
        help='Use data-driven photon+jets. Forwarded to runHZZanalysis.'
    )
    arg_parser.add_argument(
        '--local-copy', action='store_true',
        help='Copy input files to node on which the job is running.'
    )
    arg_parser.add_argument(
        '--syst', default='',
        help='Requested systematic variation or a group of them.'
    )
    args = arg_parser.parse_args()

    if args.syst == 'no':
        args.syst = ''


    analysis_options = ['--analysis=' + args.analysis]

    if args.analysis == 'Main':
        output_prefix = 'outputHZZ_'
    elif args.analysis == 'InstrMET':
        output_prefix = 'outputInstrMET_'
    elif args.analysis == 'NRB':
        output_prefix = 'outputNRB_'
    else:
        raise RuntimeError('Unrecognized analysis "{}".'.format(args.analysis))

    if args.dd_photon:
        output_prefix = 'outputPhotonDatadriven_'
        analysis_options.append('--dd-photon')

    if args.syst != 'all':
        analysis_options.append('--all-control-plots')


    job_builder = JobBuilder(
        args.task_dir, args.config, analysis_options=analysis_options,
        output_prefix=output_prefix, local_copy=args.local_copy
    )
    datasets = parse_datasets_file(args.datasets, args.config)

    if not args.syst or args.syst == 'all':
        # Nominal configuration for systematic variations
        for dataset in datasets:
            job_builder.prepare_jobs(dataset, '')

    if args.syst:
        # Systematic variations
        dataset_selector = SystDatasetSelector(
            os.path.join(job_builder.install_path, 'config/syst.yaml')
        )

        for variation, dataset in dataset_selector(datasets, args.syst):
            job_builder.prepare_jobs(dataset, variation)

    job_builder.write_submit_script()
