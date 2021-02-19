#!/usr/bin/env python

"""Makes Asimov templates from templates in the CR and SR."""

import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

CHANNELS = ["eq0jets", "eq1jets", "geq2jets_discrbin1", "geq2jets_discrbin2",
            "geq2jets_discrbin3", "geq2jets_discrbin4", "geq2jets_discrbin5",
            "geq2jets_discrbin6", "geq2jets_discrbin7"]
# Drop emu for now
PROCESSES_TO_EXCLUDE = ["data_obs", "GGToZZ_S", "GGToZZ_B", "DYJets"]


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
    arg_parser.add_argument('-o', '--output_SR',
                            default='templates_SR_Asimov.root',
                            help='Output Asimov template for SR.')
    args = arg_parser.parse_args()

    template_input_sr = ROOT.TFile(args.input_SR)
    template_input_cr = ROOT.TFile(args.input_CR)
    template_output_sr = ROOT.TFile(args.output_SR, 'recreate')
    mean_weights_file = ROOT.TFile(args.input_mean_weights)

    template_cr_substraction = ROOT.TFile("templates_CR_temp.root", 'recreate')

    template_output_sr.cd()
    for key in template_input_sr.GetListOfKeys():
        subdir = template_input_sr.Get(key.GetName())
        copy_dir(subdir)
    template_cr_substraction.cd()
    for key in template_input_cr.GetListOfKeys():
        subdir = template_input_cr.Get(key.GetName())
        copy_dir(subdir)

    for channel in CHANNELS:
        template_cr_substraction.cd(channel)
        h_instrmet = ROOT.gDirectory.Get("data_obs/nominal").Clone()
        h_genuine_met = construct_asimov(
            channel, template_cr_substraction)

        h_asimov_sr = construct_asimov(channel, template_output_sr)
        h_asimov_sr.SetName("nominal")

        # Inject instrumental pT-miss in the SR template
        h_mean_weights = mean_weights_file.Get("mean_weights_tot_"+channel)
        print("channel " + channel + ": h_instrmet has "
              + str(h_instrmet.GetNbinsX()) + " bins and h_mean_weights has "
              + str(h_mean_weights.GetNbinsX()) + " bins.")
        h_instrmet.Add(h_genuine_met, -1.)
        h_instrmet.Multiply(h_mean_weights)
        # Clip values to a minimum of 0 as a process needs to be positive.
        for i in range(1, h_instrmet.GetNbinsX()+1):
            if h_instrmet.GetBinContent(i) <= 0:
                h_instrmet.SetBinContent(i, 0)
        h_asimov_sr.Add(h_instrmet, 1.)
        h_instrmet.Delete()
        h_genuine_met.Delete()
        h_mean_weights.Delete()

        h_asimov_sr.SetName("nominal")
        template_output_sr.cd(channel+"/data_obs")
        ROOT.gDirectory.Delete("nominal;1")
        h_asimov_sr.Write()
        h_asimov_sr.Delete()

    template_output_sr.Close()
    template_cr_substraction.Close()
