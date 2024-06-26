#!/usr/bin/env python

"""Creates analysis histograms from trees."""

import argparse
from array import array
import itertools
import os
import re

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class Channel:
    """A single channel for statistical analysis.

    Attributes:
        name:             String representing name of this channel.
        selection:        String with event selection.
        mt_binning:       Binning in mt represented with an array.array.
        reweight_formula: Formula for additional reweighting to apply to both
                          data and MC in this channel, in terms of tree
                          variables (e.g., reweighting of the photon CR).
                          Default is "1".
    """

    def __init__(self, name, selection, mt_binning, reweight_formula='1'):
        """Initialize from full specification.

        The binning can be given in the form of an array.array or a
        generic iterable.  Other arguments are added directly as the
        attributes.
        """

        self.name = name
        self.selection = selection
        if isinstance(mt_binning, array):
            self.mt_binning = mt_binning
        else:
            self.mt_binning = array('d', mt_binning)
        self.reweight_formula = reweight_formula


def fill_hists(path, channels):
    """Construct templates from a ROOT file.

    Arguments:
        path:           Path to a ROOT file produced by DileptonTrees
                        analysis.
        channels:       Channels to include.

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

    data_frame = ROOT.RDataFrame(tree)
    proxies = {}
    for channel in channels:
        df_channel = data_frame.Filter(channel.selection)
        hist_model = ROOT.RDF.TH1DModel(
            '', '', len(channel.mt_binning) - 1, channel.mt_binning)
        for syst, weight_branch in syst_branches:
            if weight_branch:
                df_channel_sim = df_channel.Define(
                    'weight_sim',
                    weight_branch + '*' + channel.reweight_formula)
                proxy = df_channel_sim.Histo1D(hist_model, 'mT', 'weight_sim')
            else:
                df_channel_data = df_channel.Define(
                    'weight_data', channel.reweight_formula)
                proxy = df_channel_data.Histo1D(
                    hist_model, 'mT', 'weight_data')
            proxies[channel.name, syst] = proxy

    hists = {}
    for (channel_name, syst), proxy in proxies.items():
        # Clone the histogram because the one given by the proxy is
        # owned by the RDataFrame and will be deleted
        hist = proxy.GetValue().Clone()
        # Clip values to a minimum of 1e-6 to avoid bugs with Combine
        for i in range(1, hist.GetNbinsX()+1):
            if hist.GetBinContent(i) <= 0:
                hist.SetBinContent(i, 1e-6)
        hist.SetDirectory(ROOT.nullptr)
        hists[channel_name, syst] = hist

    input_file.Close()
    return hists


def collect_hists(directory, processes, channels):
    """Combine templates for all processes in given group.

    For each systematic variation, add all processes together.  If a
    variation is missing for a process, use the nominal template for it
    instead.

    Arguments:
        directory:   Path to directory containing ROOT files with trees.
        processes:   Names of processes included in this group.
        channels:    Channels to include.

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
            continue
        hists = fill_hists(os.path.join(directory, filename), channels)
        for (channel_name, syst), hist in hists.items():
            process_hists[process, channel_name, syst] = hist

        # Construct histograms for variations saved in separate files
        filename_regex = re.compile(f'^{process}_(.+_(up|down))\\.root$')
        filename_regex = re.compile(
            f'^{process}_(?:(?!PtG\\-130_))(?:(?!BSI_))(.+_(up|down))\\.root$')
        for filename in filenames:
            match = filename_regex.match(filename)
            if not match:
                continue
            syst = match.group(1)
            hists = fill_hists(os.path.join(directory, filename),
                               channels)
            if len(hists) != len(channels):
                raise RuntimeError(
                    'More than one systematic variation found in '
                    'file "{}".'.format(os.path.join(directory, filename)))
            for (channel_name, _), hist in hists.items():
                process_hists[process, channel_name, syst] = hist

    # For each (channel, syst), combine all processes allowing for
    # missing systematic variations in some of them
    combined_hists = {}
    all_systs = {key[2] for key in process_hists.keys()}
    for channel, syst in itertools.product(channels, all_systs):
        combined_hist = None
        for process in processes:
            hist = process_hists.get((process, channel.name, syst), None)
            if not hist:
                # The current systematic variation is not available for
                # this process.  Take the nominal template for it
                # instead.
                hist = process_hists[process, channel.name, '']
            if not combined_hist:
                combined_hist = hist.Clone()
                combined_hist.SetName('')
            else:
                combined_hist.Add(hist)
        combined_hists[channel.name, syst] = combined_hist

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
            ('V', ['DYJets', 'WJets', 'ZJets']),
            ('VV', ['WW', 'WZ', 'ZZ']),
            ('TT', ['TT', 'TTV']),
            ('ST', ['ST']),
            ('VVV', ['VVV']),
            ('WG', ['WG']),
            ('ZG', ['ZG']),
            ('GGToZZ', ['GGToZZ_S', 'GGToZZ_B', 'GGToZZ_BSI'])
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
            New label for the given systematic variation.  If there is
            no renaming rule for the given systematic variation, its
            base label is left unchanged.  The "_up"/"_down" postfix is
            always replaced with "Up"/"Down" respectively.
        """

        match = self._syst_regex.match(syst)
        if not match:
            raise RuntimeError(
                f'Failed to parse systematic variation "{syst}".')

        syst_base = self.rename_map.get(
            (template_name, match.group(1)), match.group(1))
        if syst_base in ['jec', 'jer', 'metuncl']:
            syst_base = syst_base + '_' + args.year
        return f'{syst_base}{match.group(2).capitalize()}'


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('directory', help='Directory with trees.')
    arg_parser.add_argument('-o', '--output', default='templates.root',
                            help='Name for output file with histograms.')
    arg_parser.add_argument('-a', '--analysis', choices={'dilepton', 'photon'},
                            default='dilepton',
                            help='Analysis name. Default is "dilepton".')
    arg_parser.add_argument('-c', '--channel', default='all',
                            help='channel on which to run. Put "all"'
                            '(default) if you want to run on all channels.')
    arg_parser.add_argument('-y', '--year', default='2016',
                            help='Year. Used to decorrelate systematics.')
    arg_parser.add_argument('--nrb', default=False, action='store_true',
                            help='results will be NRB data-driven.')
    args = arg_parser.parse_args()

    geq1jets_binning = [
        100, 200, 300, 350, 400, 450, 500, 550, 600, 700, 850, 1000,
        1250, 1500, 3000]
    vbf_binning = [150, 225, 300, 375, 450, 600, 750, 1100, 3000]
    if args.analysis == 'dilepton':
        channels = [
            Channel(
                'eq0jets', 'lepton_cat != 2 && jet_cat == 0 && ptmiss > 125.',
                geq1jets_binning),
            Channel(
                'eq1jets', 'lepton_cat != 2 && jet_cat == 1 && ptmiss > 125.',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin1', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0. && sm_DjjVBF < 0.05',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin2', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0.05 && sm_DjjVBF < 0.1',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin3', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0.1 && sm_DjjVBF < 0.2',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin4', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0.2 && sm_DjjVBF < 0.8',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin5', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0.8 && sm_DjjVBF < 0.9',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin6', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0.9 && sm_DjjVBF < 0.95',
                geq1jets_binning),
            Channel(
                'geq2jets_discrbin7', 'lepton_cat != 2 && jet_cat == 2 && '
                'ptmiss > 125. && sm_DjjVBF >= 0.95', geq1jets_binning),
            # Event counting in the emu control region.  Use a finite range
            # instead of (-inf, inf) to allow inspection in TBrowser.
            Channel('emu', 'lepton_cat == 2 && ptmiss > 80.', [0., 1e4])
        ]
        if args.nrb:
            channels_nrb = [
                Channel(
                    'eq0jets',
                    'lepton_cat == 2 && jet_cat == 0 && ptmiss > 125.',
                    geq1jets_binning),
                Channel(
                    'eq1jets',
                    'lepton_cat == 2 && jet_cat == 1 && ptmiss > 125.',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin1',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0. && sm_DjjVBF < 0.05',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin2',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0.05 && sm_DjjVBF < 0.1',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin3',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0.1 && sm_DjjVBF < 0.2',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin4',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0.2 && sm_DjjVBF < 0.8',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin5',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0.8 && sm_DjjVBF < 0.9',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin6',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0.9 && sm_DjjVBF < 0.95',
                    geq1jets_binning),
                Channel(
                    'geq2jets_discrbin7',
                    'lepton_cat == 2 && jet_cat == 2 && ptmiss > 125. && '
                    'sm_DjjVBF >= 0.95',
                    geq1jets_binning),
            ]
    elif args.analysis == 'photon':
        weights_photon = 'photon_reweighting * trigger_weight / mean_weight'
        channels = [
            Channel(
                'eq0jets', 'jet_cat == 0 && ptmiss > 125.',
                geq1jets_binning, weights_photon),
            Channel(
                'eq1jets', 'jet_cat == 1 && ptmiss > 125.',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin1', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0. && sm_DjjVBF < 0.05',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin2', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0.05 && sm_DjjVBF < 0.1',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin3', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0.1 && sm_DjjVBF < 0.2',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin4', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0.2 && sm_DjjVBF < 0.8',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin5', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0.8 && sm_DjjVBF < 0.9',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin6', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0.9 && sm_DjjVBF < 0.95',
                geq1jets_binning, weights_photon),
            Channel(
                'geq2jets_discrbin7', 'jet_cat == 2 && ptmiss > 125. && '
                'sm_DjjVBF >= 0.95', geq1jets_binning, weights_photon),
        ]
    else:
        raise RuntimeError(
            f'Unaccepted "analysis" option.')

    if args.channel != 'all':
        for channel in channels:
            if channel.name == args.channel:
                good_channel = channel
                channels = [good_channel]
        if args.nrb:
            for channel in channels_nrb:
                if channel.name == args.channel:
                    good_channel = channel
                    channels_nrb = [good_channel]

    output_file = ROOT.TFile(args.output, 'recreate')
    for channel in channels:
        output_file.mkdir(channel.name)
    output_hists = []

    ROOT.ROOT.EnableImplicitMT()
    syst_rename = SystRename()
    all_processes = [
        ('data_obs', ['Data']),
        ('DYJets', ['DYJetsToLL_M-50']),
        ('GGToZZ_S', ['GGToHToZZ']),
        ('GGToZZ_B', ['GGToZZ']),
        ('GGToZZ_BSI', ['GGToZZ_BSI']),
        ('ST', ['ST_s-channel', 'ST_t-channel', 'ST_tW']),
        ('TT', ['TT']),
        ('TTV', ['TTWJetsToLNu', 'TTZToLLNuNu_M-10']),
        ('WJets', ['WJets']),
        ('WW', ['WWTo2L2Nu']),
        ('WZ', ['WZTo2L2Q', 'WZTo3LNu']),
        ('ZZ', ['ZZTo2L2Nu', 'ZZTo2L2Q', 'ZZTo4L']),
        ('VVV', ['WWW', 'WWZ', 'WZZ', 'ZZZ']),
        ('WG', ['WGToLNuG']),
        ('ZG', ['ZGTo2NuG', 'ZGTo2NuG_PtG-130']),
        ('ZJets', ['ZJetsToNuNu'])
    ] if not args.nrb else[
        ('data_obs', ['Data']),
        ('DYJets', ['DYJetsToLL_M-50']),
        ('GGToZZ_S', ['GGToHToZZ']),
        ('GGToZZ_B', ['GGToZZ']),
        ('GGToZZ_BSI', ['GGToZZ_BSI']),
        ('TTV', ['TTZToLLNuNu_M-10']),
        ('WZ', ['WZTo2L2Q', 'WZTo3LNu']),
        ('ZZ', ['ZZTo2L2Nu', 'ZZTo2L2Q', 'ZZTo4L']),
        ('VVV', ['WWZ', 'WZZ', 'ZZZ']),
        ('WG', ['WGToLNuG']),
        ('ZG', ['ZGTo2NuG', 'ZGTo2NuG_PtG-130']),
        ('ZJets', ['ZJetsToNuNu']),
    ]
    timer = ROOT.TStopwatch()
    timer.Start()
    for group_name, processes in all_processes:
        totaltime = timer.RealTime()/60.
        timer.Continue()
        fraction = (all_processes.index((group_name,
                    processes))+1) / len(all_processes)
        remaining = totaltime*(1-fraction)/fraction
        print("Done {}/{} Processing {}  "
              "(elapsed {:.1f} min, remaining {:.1f} min)".format(
                int(all_processes.index((group_name, processes))),
                len(all_processes), group_name, totaltime, remaining))
        dirs = {}
        hists = collect_hists(args.directory, processes, channels)
        for channel, syst in sorted(hists.keys()):
            path = '/'.join([channel, group_name])
            if not output_file.GetDirectory(path):
                output_file.mkdir(path)
            dirs[channel] = output_file.Get(path)
            hist = hists[channel, syst]
            if not syst:
                hist.SetName('nominal')
            else:
                hist.SetName(syst_rename(group_name, syst))
            hist.SetDirectory(dirs[channel])
            output_hists.append(hist)

    if args.nrb:
        hists_nrb = collect_hists(args.directory, ['NRB'], channels_nrb)
        for channel, syst in sorted(hists_nrb.keys()):
            path = '/'.join([channel, "NRB"])
            if not output_file.GetDirectory(path):
                output_file.mkdir(path)
            dir_nrb = output_file.Get(path)
            hist = hists_nrb[channel, syst]
            if not syst:
                hist.SetName('nominal')
            else:
                hist.SetName(syst_rename(group_name, syst))
            hist.SetDirectory(dir_nrb)
            output_hists.append(hist)

    output_file.Write()
    output_file.Close()
