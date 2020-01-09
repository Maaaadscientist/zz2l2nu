#!/usr/bin/env python

"""Creates analysis histograms from trees."""

import argparse
from array import array
import itertools
import os
import re

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


CHANNELS = {'eq0jets': 0, 'geq1jets': 1, 'vbf': 2}
MT_BINNING = {
    'geq1jets': array('d', [150, 225, 300, 375, 450, 525, 600, 725, 850, 975,
                            1100, 1350, 1600, 2100, 3000]),
    'vbf': array('d', [150, 225, 300, 375, 450, 600, 750, 1100, 3000])
}
MT_BINNING['eq0jets'] = MT_BINNING['geq1jets'][1:]


def fill_hists(path):
    """Construct templates from a ROOT file.
    
    Arguments:
        path:  Path to a ROOT file produced by DileptonTrees analysis.

    Return value:
        ROOT histograms of the analysis observable for the nominal case
        and systematic variations represented by alternative weights.
        The histograms are organized into a mapping from pairs of labels
        (channel, syst).  The label for the central variation is ''.
    """

    input_file = ROOT.TFile(path)
    tree = input_file.Get('Vars')

    syst_branches = []
    for branch in tree.GetListOfBranches():
        branch_name = branch.GetName()
        if branch_name == 'weight':
            syst_branches.append(('', branch_name))
        elif branch_name.startswith('weight_'):
            syst = branch_name[len('weight_'):]
            syst_branches.append((syst, branch_name))
    if not syst_branches:
        # There are no weight branches.  This must be real data.
        syst_branches.append(('', None))

    data_frame = ROOT.RDataFrame(tree).Define('ptmiss', 'p4Miss.Pt()')\
        .Filter('leptonCat != 2 && ptmiss > 125.')
    proxies = {}
    for channel, channel_index in CHANNELS.items():
        df_channel = data_frame.Filter(f'jetCat == {channel_index}')
        hist_model = ROOT.RDF.TH1DModel(
            '', '', len(MT_BINNING[channel]) - 1, MT_BINNING[channel])
        for syst, weight_branch in syst_branches:
            if weight_branch:
                proxy = df_channel.Histo1D(hist_model, 'mT', weight_branch)
            else:
                proxy = df_channel.Histo1D(hist_model, 'mT')
            proxies[channel, syst] = proxy

    hists = {}
    for (channel, syst), proxy in proxies.items():
        # Clone the histogram because the one given by the proxy is
        # owned by the RDataFrame and will be deleted
        hist = proxy.GetValue().Clone()
        hist.SetDirectory(None)
        hists[channel, syst] = hist

    input_file.Close()
    return hists


def collect_hists(directory, processes):
    """Combine templates for all processes in given group.

    For each systematic variation, add all processes together.  If a
    variation is missing for a process, use the nominal template for it
    instead.

    Arguments:
        directory:   Path to directory containing ROOT files with trees.
        processes:   Names of processes included in this group.

    Return value:
        Mapping (channel, syst) -> template.
    """

    # Process each file independently since different processes are not
    # guaranteed to have the same sets of systematic variations.  Fill
    # mapping (process, channel, syst) -> histogram.
    process_hists = {}
    filenames = os.listdir(directory)
    for process in processes:
        # Construct nominal histogram and weight-only variations
        if process + '_weights.root' in filenames:
            filename = process + '_weights.root'
        elif process + '.root' in filenames:
            filename = process + '.root'
        else:
            raise RuntimeError(
                f'Nominal file not found for process {process}.')
        hists = fill_hists(os.path.join(directory, filename))
        for (channel, syst), hist in hists.items():
            process_hists[process, channel, syst] = hist

        # Construct histograms for variations saved in separate files
        filename_regex = re.compile(f'^{process}_(.+_(up|down))\.root$')
        for filename in filenames:
            match = filename_regex.match(filename)
            if not match:
                continue
            syst = match.group(1)
            hists = fill_hists(os.path.join(directory, filename))
            if len(hists) != len(CHANNELS):
                raise RuntimeError(
                    'More than one systematic variation found in '
                    'file "{}".'.format(os.path.join(directory, filename)))
            for (channel, _), hist in hists.items():
                process_hists[process, channel, syst] = hist

    # For each (channel, syst), combine all processes allowing for
    # missing systematic variations in some of them
    combined_hists = {}
    all_systs = {key[2] for key in process_hists.keys()}
    for channel, syst in itertools.product(CHANNELS, all_systs):
        combined_hist = None
        for process in processes:
            hist = process_hists.get((process, channel, syst), None)
            if not hist:
                # The current systematic variation is not available for
                # this process.  Take the nominal template for it
                # instead.
                hist = process_hists[process, channel, '']
            if not combined_hist:
                combined_hist = hist.Clone()
                combined_hist.SetName('')
            else:
                combined_hist.Add(hist)
        combined_hists[channel, syst] = combined_hist

    return combined_hists


class SystRename:
    """Class implementing renaming rules for systematic variations."""

    def __init__(self):
        # Map from (template_name, syst_base) to new base labels for
        # systematic variations.  Base labels are defined by removing
        # the "_up"/"_down" postfix.
        self.rename_map = {}

        # Define correlations for renormalization and factorization
        # scales
        for correlation_group, template_names in [
            ('V', ['DYJets', 'WJets']),
            ('VV', ['WW', 'WZ', 'ZZ']),
            ('TT', ['TT', 'TTV']),
            ('ST', ['ST']),
            ('VVV', ['VVV'])
        ]:
            for template_name, syst in itertools.product(
                template_names, ['me_renorm', 'factor']
            ):
                self.rename_map[template_name, syst] = \
                    f'{syst}_{correlation_group}'

        self._syst_regex = re.compile('^(.+)_(up|down)$')


    def __call__(self, template_name, syst):
        """Find new name for given variation.

        Arguments:
            template_name:  Name of process group, i.e. the template.
            syst:           Original label of systematic variation.

        Return value:
            New label for the given systematic variation.  The name is
            returned unchanged if there is no renaming rule for it.
        """

        match = self._syst_regex.match(syst)
        if not match:
            return syst

        syst_base = self.rename_map.get(
            (template_name, match.group(1)), match.group(1))
        return f'{syst_base}_{match.group(2)}'



if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('directory', help='Directory with trees.')
    arg_parser.add_argument('-o', '--output', default='templates.root',
                            help='Name for output file with histograms.')
    args = arg_parser.parse_args()

    output_file = ROOT.TFile(args.output, 'recreate')
    channel_dirs = {}
    for channel in CHANNELS:
        channel_dirs[channel] = output_file.mkdir(channel)
    output_hists = []

    ROOT.ROOT.EnableImplicitMT()
    syst_rename = SystRename()
    for group_name, processes in [
        ('data_obs', ['Data']),
        ('DYJets', ['DYJetsToLL_M-50']),
        ('GGToZZ_S', ['GGToHToZZTo2E2Nu', 'GGToHToZZTo2Mu2Nu']),
        ('GGToZZ_BSI', ['GGToZZTo2E2Nu_BSI', 'GGToZZTo2Mu2Nu_BSI']),
        ('GGToZZ_B', ['GluGluToContinToZZTo2e2nu',
                      'GluGluToContinToZZTo2mu2nu']),
        ('ST', ['ST_s-channel', 'ST_t-channel_antitop', 'ST_t-channel_top',
                'ST_tW_antitop', 'ST_tW_top']),
        ('TT', ['TT']),
        ('TTV', ['TTWJetsToLNu', 'TTZToLLNuNu_M-10']),
        ('WJets', [
            'WJetsToLNu_LO', 'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-200To400', 'WJetsToLNu_HT-400To600',
            'WJetsToLNu_HT-600To800', 'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-1200To2500', 'WJetsToLNu_HT-2500ToInf'
        ]),
        ('WW', ['WWTo2L2Nu']),
        ('WZ', ['WZTo2L2Q', 'WZTo3LNu']),
        ('ZZ', ['ZZTo2L2Nu', 'ZZTo2L2Q', 'ZZTo4L']),
        ('VVV', ['WWW', 'WWZ', 'WZZ', 'ZZZ'])
    ]:
        hists = collect_hists(args.directory, processes)
        for channel, syst in sorted(hists.keys()):
            hist = hists[channel, syst]
            if not syst:
                hist.SetName(group_name)
            else:
                hist.SetName('{}_{}'.format(
                    group_name, syst_rename(group_name, syst)))
            hist.SetDirectory(channel_dirs[channel])
            output_hists.append(hist)

    output_file.Write()
    output_file.Close()

