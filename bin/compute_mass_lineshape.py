#!/usr/bin/env python

"""Computes mass lineshape from analysis trees, to give a mass to photons"""

import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class Channel:
    """A single channel to compute the weights.

    Attributes:
        name:      String representing name of this channel.
        selection: String with event selection.
    """

    def __init__(self, name, selection):
        """Initialize from full specification.
        """

        self.name = name
        self.selection = selection


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('dir_dilepton',
                            help='Directory with dilepton trees.')
    arg_parser.add_argument('-o', '--output', default='lineshape_mass.root',
                            help='Name for output file.')
    args = arg_parser.parse_args()

    channels = [
        Channel(
            'eq0jets_ll', 'lepton_cat != 2 && jet_cat == 0 && ptmiss <=125'),
        Channel(
            'eq1jets_ll', 'lepton_cat != 2 && jet_cat == 1 && ptmiss <=125'),
        Channel(
            'geq2jets_ll', 'lepton_cat != 2 && jet_cat == 2 && ptmiss <=125'),
        Channel(
            'll', 'lepton_cat != 2 && ptmiss <=125')
    ]

    output_file = ROOT.TFile(args.output, 'recreate')
    dilepton_file = ROOT.TFile(args.dir_dilepton + '/Data.root')

    dilepton_tree = dilepton_file.Get('Vars')
    dilepton_data_frame = ROOT.RDataFrame(dilepton_tree)

    proxies = {}

    for channel in channels:
        dilepton_df_channel = dilepton_data_frame.Filter(channel.selection)
        hist_model = ROOT.RDF.TH1DModel('', '', 100, 76, 106)
        proxy_dilepton = dilepton_df_channel.Histo1D(hist_model, 'll_mass')
        proxies[channel.name] = proxy_dilepton

    output_file.cd()
    for (channel_name), proxy in proxies.items():
        hist = proxy.GetValue().Clone()
        hist.SetName("WeightHisto_" + channel_name)
        hist.Write()
        hist.Delete()

    output_file.Close()
