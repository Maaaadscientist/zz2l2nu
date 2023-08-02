#!/usr/bin/env python

"""Computes weights in nvtx, eta and pT from analysis trees"""

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


def smoothen_histo(histo):
    """Smoothen a histogram by making a fit.

    A simple linear function is taken for nvtx, and a 3rd order polynomial is
    used for pT, above 180 GeV only (the rest is kept at its original value).

    Arguments:
        histo: histogram with the weights
    """
    func_nvtx = ROOT.TF1("func_nvtx", "[0]+[1]*x", 0, 100.)
    histo.Fit("func_nvtx")
    for i in range(0, histo.GetNbinsX() + 1):
        histo.SetBinContent(i, func_nvtx.Eval(histo.GetBinCenter(i)))
    return histo


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('dir_dilepton',
                            help='Directory with dilepton trees.')
    arg_parser.add_argument('dir_photon', help='Directory with photon trees.')
    arg_parser.add_argument('-o', '--output', default='weights.root',
                            help='Name for output file.')
    arg_parser.add_argument('-s', '--step', choices={'nvtx', 'eta', 'pt'},
                            default='nvtx',
                            help='Reweighting step to use. "nvtx" must be run'
                            'before "eta" and "pt" (it is necessary to re-run'
                            'the photon trees with other weights already'
                            'applied).')
    args = arg_parser.parse_args()

    channels = [
        Channel(
            'eq0jets_ee', 'lepton_cat == 0 && jet_cat == 0 && ptmiss <=120',
            'jet_cat == 0 && ptmiss <=120'),
        Channel(
            'eq1jets_ee', 'lepton_cat == 0 && jet_cat == 1 && ptmiss <=120',
            'jet_cat == 1 && ptmiss <=120'),
        Channel(
            'geq2jets_ee', 'lepton_cat == 0 && jet_cat == 2 && ptmiss <=120',
            'jet_cat == 2 && ptmiss <=120'),
        Channel(
            'ee', 'lepton_cat == 0 && ptmiss <=120', 'ptmiss <=120'),
        Channel(
            'eq0jets_mumu', 'lepton_cat == 1 && jet_cat == 0 && ptmiss <=120',
            'jet_cat == 0 && ptmiss <=120'),
        Channel(
            'eq1jets_mumu', 'lepton_cat == 1 && jet_cat == 1 && ptmiss <=120',
            'jet_cat == 1 && ptmiss <=120'),
        Channel(
            'geq2jets_mumu', 'lepton_cat == 1 && jet_cat == 2 && ptmiss <=120',
            'jet_cat == 2 && ptmiss <=120'),
        Channel(
            'mumu', 'lepton_cat == 1 && ptmiss <=120', 'ptmiss <=120'),
        Channel(
            'eq0jets_ll', 'lepton_cat != 2 && jet_cat == 0 && ptmiss <=120',
            'jet_cat == 0 && ptmiss <=120'),
        Channel(
            'eq1jets_ll', 'lepton_cat != 2 && jet_cat == 1 && ptmiss <=120',
            'jet_cat == 1 && ptmiss <=120'),
        Channel(
            'geq2jets_ll', 'lepton_cat != 2 && jet_cat == 2 && ptmiss <=120',
            'jet_cat == 2 && ptmiss <=120'),
        Channel(
            'll', 'lepton_cat != 2 && ptmiss <=120', 'ptmiss <=120')
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
            photon_df_channel = photon_df_channel.Define(
                "photon_weight", 'trigger_weight * beam_halo_weight')
            proxy_photon = photon_df_channel.Histo1D(
                hist_model, 'num_pv_good', 'photon_weight')
            proxies[channel.name, 'dilepton'] = proxy_dilepton
            proxies[channel.name, 'photon'] = proxy_photon

        if args.step == 'eta':
            hist_model = ROOT.RDF.TH1DModel(
                '', ';boson #eta;weight', 24, 0, 2.4)
            dilepton_df_channel = dilepton_df_channel.Define(
                "ll_abs_eta", 'fabs(ll_eta)')
            proxy_dilepton = dilepton_df_channel.Histo1D(
                hist_model, 'll_abs_eta')
            photon_df_channel = photon_df_channel.Define(
                "photon_weight", 'photon_nvtx_reweighting * trigger_weight'
                ' * beam_halo_weight')
            photon_df_channel = photon_df_channel.Define(
                "photon_abs_eta", 'fabs(photon_eta)')
            proxy_photon = photon_df_channel.Histo1D(
                hist_model, 'photon_abs_eta', 'photon_weight')
            proxies[channel.name, 'dilepton'] = proxy_dilepton
            proxies[channel.name, 'photon'] = proxy_photon

        if args.step == 'pt':
            pt_binning = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165,
                          180, 250, 1500]
            pt_binning = array('d', pt_binning)
            hist_model = ROOT.RDF.TH1DModel(
                '', ';boson p_{T} (GeV);weight', len(pt_binning) - 1,
                pt_binning)
            proxy_dilepton = dilepton_df_channel.Histo1D(hist_model, 'll_pt')
            photon_df_channel = photon_df_channel.Define(
                "photon_weight", 'photon_nvtx_reweighting'
                ' * photon_eta_reweighting * trigger_weight'
                ' * beam_halo_weight')
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
        if args.step == 'nvtx' or args.step == 'eta':
            hists[channel.name, 'dilepton'].Scale(
                1./hists[channel.name, 'dilepton'].Integral())
            hists[channel.name, 'photon'].Scale(
                1./hists[channel.name, 'photon'].Integral())
        h_reweighting = hists[channel.name, 'dilepton'].Clone()
        h_reweighting.Divide(hists[channel.name, 'photon'])
        h_reweighting.SetName("WeightHisto_" + channel.name)

        if args.step == 'nvtx':
            h_reweighting = smoothen_histo(h_reweighting)

        h_reweighting.Write()
        h_reweighting.Delete()

    output_file.Close()
