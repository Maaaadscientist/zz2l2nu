#!/usr/bin/env python

import argparse
import math
import os


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument(
        'path', help='Path to the merged trees, from HZZ2L2NU_BASE'
    )
    arg_parser.add_argument(
        '-o', '--output', help='name of the output template (without extension)'
    )
    arg_parser.add_argument(
        '-s', '--selection', help='"dilepton" or "photon"'
    )
    arg_parser.add_argument(
        '-y', '--year', help='year'
    )
    args = arg_parser.parse_args()

    channels = ["eq0jets", "eq1jets", "geq2jets_discrbin1",
                "geq2jets_discrbin2", "geq2jets_discrbin3",
                "geq2jets_discrbin4", "geq2jets_discrbin5",
                "geq2jets_discrbin6", "geq2jets_discrbin7"]

    nrb_option = ''
    if args.selection == 'dilepton':
        nrb_option = '--nrb'

    path_to_scripts = os.environ['HZZ2L2NU_BASE'] + '/templates/scripts/' + args.output + '/'
    if os.path.exists(path_to_scripts):
        raise RuntimeError(f'Path already exists.')
    else:
        os.makedirs(path_to_scripts)
        os.makedirs(path_to_scripts + 'logs/')

    for channel in channels:
        script_name = path_to_scripts + 'run_' + channel + '.sh'
        script_commands = [
          'export INITDIR={}'.format(os.environ['HZZ2L2NU_BASE']),
          'cd $INITDIR',
          '. ./env.sh',
          'cd -',
          'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi',
          'hostname',
          'date'
        ]
        script_commands.append('build_templates.py -o {}.root -a {} -c {} -y {} {} {}'.format(os.environ['HZZ2L2NU_BASE'] + '/templates/' + args.output + "_" + channel, args.selection, channel, args.year, nrb_option, os.environ['HZZ2L2NU_BASE'] + '/' + args.path))

        with open(script_name, 'w') as f:
            for command in script_commands:
                f.write(command)
                f.write('\n')
    submit_script = path_to_scripts + 'submit.sh'
    submission_lines = []
    for channel in channels:
        submission_lines.append('qsub -l walltime=01:00:00 -j oe -o {}logs/ {}'.format(path_to_scripts, path_to_scripts + 'run_' + channel + '.sh'))

    with open(submit_script, 'w') as f:
        for line in submission_lines:
            f.write(line)
            f.write('\n')
