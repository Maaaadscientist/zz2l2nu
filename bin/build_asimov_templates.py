#!/usr/bin/env python

"""Makes Asimov templates from templates in the CR and SR."""

import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

CHANNELS = ["eq0jets", "geq1jets", "vbf"]  # Drop emu for now
PROCESSES_TO_EXCLUDE = ["data_obs", "GGToZZ_S", "GGToZZ_B"]


def copy_dir(source):
    """Copy (recursively) the contents of a given source (sub-directory)"""
    adir = ROOT.gDirectory.mkdir(source.GetName())
    adir.cd()

    for key in source.GetListOfKeys():
        classname = key.GetClassName()
        cl = ROOT.gROOT.GetClass(classname)
        if (not cl):
            continue
        if (cl.InheritsFrom(ROOT.TDirectory.Class())):
            subdir = source.Get(key.GetName())
            adir.cd()
            copy_dir(subdir)
            adir.cd()
        else:
            source.cd()
            obj = key.ReadObj()
            adir.cd()
            obj.Write()
            obj.Delete()

    adir.SaveSelf(ROOT.kTRUE)
    adir.GetMotherDir().cd()


def construct_asimov(channel, template):
    """Construct a histogram with Asimov data in a given channel.

    These Asimov data do not include the "InstrMET" for the CR template.

    Arguments:
        channel:  Current channel
        template: Template (for SR or CR)
    """
    template.cd(channel)
    h_asimov = ROOT.gDirectory.Get("data_obs/nominal")
    h_asimov.Reset()
    for key in ROOT.gDirectory.GetListOfKeys():
        subdir = key.GetName()
        if subdir in PROCESSES_TO_EXCLUDE:
            print("Not taken: " + subdir)
            continue
        else:
            h_nominal = ROOT.gDirectory.Get(subdir+"/nominal")
            h_asimov += h_nominal
            h_nominal.Delete()
    return h_asimov


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('input_SR', help='Input template for SR.')
    arg_parser.add_argument('input_CR', help='Input template for CR.')
    arg_parser.add_argument('-w', '--input_mean_weights', default='\
                            data/InstrMetReweighting/meanWeights_2017.root',
                            help='Input file for mean weights.')
    arg_parser.add_argument('-s', '--output_SR',
                            default='templates_SR_Asimov.root',
                            help='Output Asimov template for SR.')
    arg_parser.add_argument('-c', '--output_CR',
                            default='templates_CR_Asimov.root',
                            help='Output Asimov template for CR.')
    args = arg_parser.parse_args()

    template_input_sr = ROOT.TFile(args.input_SR)
    template_input_cr = ROOT.TFile(args.input_CR)
    template_output_sr = ROOT.TFile(args.output_SR, 'recreate')
    template_output_cr = ROOT.TFile(args.output_CR, 'recreate')
    mean_weights_file = ROOT.TFile(args.input_mean_weights)

    template_output_sr.cd()
    for key in template_input_sr.GetListOfKeys():
        subdir = template_input_sr.Get(key.GetName())
        copy_dir(subdir)
    template_output_cr.cd()
    for key in template_input_cr.GetListOfKeys():
        subdir = template_input_cr.Get(key.GetName())
        copy_dir(subdir)

    for channel in CHANNELS:
        h_asimov_sr = construct_asimov(channel, template_output_sr)
        h_asimov_sr.SetName("nominal")
        template_output_sr.cd(channel+"/data_obs")
        ROOT.gDirectory.Delete("nominal;1")
        h_asimov_sr.Write()
        h_asimov_sr.Delete()

        h_asimov_cr = construct_asimov(channel, template_output_cr)

        # Inject instrumental pT-miss in the CR template
        template_output_sr.cd(channel)
        h_instrmet = ROOT.gDirectory.Get("DYJets/nominal")
        h_mean_weights = mean_weights_file.Get("mean_weights_tot_"+channel)
        print("channel " + channel + ": h_instrmet has "
              + str(h_instrmet.GetNbinsX()) + " bins and h_mean_weights has "
              + str(h_mean_weights.GetNbinsX()) + " bins.")
        h_instrmet.Divide(h_mean_weights)
        template_output_cr.cd(channel)
        h_asimov_cr += h_instrmet
        h_instrmet.Delete()
        h_mean_weights.Delete()

        h_asimov_cr.SetName("nominal")
        template_output_cr.cd(channel+"/data_obs")
        ROOT.gDirectory.Delete("nominal;1")
        h_asimov_cr.Write()
        h_asimov_cr.Delete()

    template_output_sr.Close()
    template_output_cr.Close()
