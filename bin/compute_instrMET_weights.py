#!/usr/bin/env python

"""Computes weights in nvtx and pT from analysis trees"""

import argparse
from array import array

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class Channel:
    """A single channel to compute the weights.

    Attributes:
        name:               String representing name of this channel.
        selection_dilepton: String with event selection for dilepton data.
        selection_photon:   String with event selection for photon data.
    """

    def __init__(self, name, selection_dilepton, selection_photon):
        """Initialize from full specification.
        """

        self.name = name
        self.selection_dilepton = selection_dilepton
        self.selection_photon = selection_photon


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('dir_dilepton',
                            help='Directory with dilepton trees.')
    arg_parser.add_argument('dir_photon', help='Directory with photon trees.')
    arg_parser.add_argument('-o', '--output', default='weights.root',
                            help='Name for output file.')
    arg_parser.add_argument('-s', '--step', choices={'nvtx', 'pt'},
                            default='nvtx',
                            help='Reweighting step to use. "nvtx" must be run'
                            'before "pt" (it is necessary to re-run the photon'
                            'trees with nvtx weights already applied).')
    args = arg_parser.parse_args()

    channels = [
        Channel(
            'eq0jets_ll', 'lepton_cat != 2 && jet_cat == 0 && ptmiss <=125',
            'jet_cat == 0 && ptmiss <=125'),
        Channel(
            'eq1jets_ll', 'lepton_cat != 2 && jet_cat == 1 && ptmiss <=125',
            'jet_cat == 1 && ptmiss <=125'),
        Channel(
            'geq2jets_ll', 'lepton_cat != 2 && jet_cat == 2 && ptmiss <=125',
            'jet_cat == 2 && ptmiss <=125'),
        Channel(
            'll', 'lepton_cat != 2 && ptmiss <=125', 'ptmiss <=125')
    ]

    output_file = ROOT.TFile(args.output, 'recreate')
    dilepton_file = ROOT.TFile(args.dir_dilepton + '/Data.root')
    photon_file = ROOT.TFile(args.dir_photon + '/Data.root')

    dilepton_tree = dilepton_file.Get('Vars')
    photon_tree = photon_file.Get('Vars')
    dilepton_data_frame = ROOT.RDataFrame(dilepton_tree)
    photon_data_frame = ROOT.RDataFrame(photon_tree)

    proxies = {}

    for channel in channels:
        dilepton_df_channel = dilepton_data_frame.Filter(
            channel.selection_dilepton)
        photon_df_channel = photon_data_frame.Filter(channel.selection_photon)

        if args.step == 'nvtx':
            hist_model = ROOT.RDF.TH1DModel(
                '', ';number of good vertices;weight', 100, 0, 100)
            proxy_dilepton = dilepton_df_channel.Histo1D(
                hist_model, 'num_pv_good')
            proxy_photon = photon_df_channel.Histo1D(
                hist_model, 'num_pv_good', 'trigger_weight')
            proxies[channel.name, 'dilepton'] = proxy_dilepton
            proxies[channel.name, 'photon'] = proxy_photon

        if args.step == 'pt':
            pt_binning = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165,
                          180, 195, 210, 225, 240, 255, 270, 285, 300, 315,
                          330, 345, 360, 375, 390, 405, 435, 465, 495, 525,
                          555, 585, 615, 675, 735, 795, 855, 975, 1500]
            pt_binning = array('d', pt_binning)
            hist_model = ROOT.RDF.TH1DModel(
                '', ';boson p_{T} (GeV);weight', len(pt_binning) - 1,
                pt_binning)
            proxy_dilepton = dilepton_df_channel.Histo1D(hist_model, 'll_pt')
            photon_df_channel = photon_df_channel.Define(
                "photon_weight", 'photon_nvtx_reweighting * trigger_weight')
            proxy_photon = photon_df_channel.Histo1D(
                hist_model, 'photon_pt', 'photon_weight')
            proxies[channel.name, 'dilepton'] = proxy_dilepton
            proxies[channel.name, 'photon'] = proxy_photon

    hists = {}
    for (channel_name, selection), proxy in proxies.items():
        hist = proxy.GetValue().Clone()
        hists[channel_name, selection] = hist

    output_file.cd()

    for channel in channels:
        if args.step == 'nvtx':
            hists[channel.name, 'dilepton'].Scale(
                1./hists[channel.name, 'dilepton'].Integral())
            hists[channel.name, 'photon'].Scale(
                1./hists[channel.name, 'photon'].Integral())
        h_reweighting = hists[channel.name, 'dilepton'].Clone()
        h_reweighting.Divide(hists[channel.name, 'photon'])
        h_reweighting.SetName("WeightHisto_" + channel.name)

        h_reweighting.Write()
        h_reweighting.Delete()

    output_file.Close()
