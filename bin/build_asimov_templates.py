#!/usr/bin/env python

"""Makes Asimov templates from templates in the CR and SR."""

import argparse
from array import array

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

CHANNELS = ["eq0jets", "geq1jets", "vbf"]  # Drop emu for now
PROCESSES_TO_EXCLUDE = ["data_obs", "GGToZZ_S", "GGToZZ_B"]

def copy_dir(source):
    """Copy (recursively) the contents of a given source (sub-directory)
    """
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

    Aruments:
        channel:  Current channel
        template: Template (for SR or CR)
    """
    template.cd(channel)
    h_asimov = ROOT.gDirectory.Get("data_obs/nominal")
    h_asimov.Reset()
    for key in ROOT.gDirectory.GetListOfKeys():
        subdir = key.GetName()
        if subdir in PROCESSES_TO_EXCLUDE:
            print ("Not taken: " + subdir)
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
    arg_parser.add_argument('-w', '--input_mean_weights', default=
                            'data/InstrMetReweighting/meanWeights_2017.root',
                            help='Input file for mean weights.')
    arg_parser.add_argument('-s', '--output_SR', 
                            default='templates_SR_Asimov.root',
                            help='Output Asimov template for SR.')
    arg_parser.add_argument('-c', '--output_CR', 
                            default='templates_CR_Asimov.root',
                            help='Output Asimov template for CR.')
    args = arg_parser.parse_args()

    template_input_SR = ROOT.TFile(args.input_SR)
    template_input_CR = ROOT.TFile(args.input_CR)
    template_output_SR = ROOT.TFile(args.output_SR, 'recreate')
    template_output_CR = ROOT.TFile(args.output_CR, 'recreate')
    mean_weights_file = ROOT.TFile(args.input_mean_weights)

    template_output_SR.cd()
    for key in template_input_SR.GetListOfKeys():
        subdir = template_input_SR.Get(key.GetName())
        copy_dir(subdir)
    template_output_CR.cd()
    for key in template_input_CR.GetListOfKeys():
        subdir = template_input_CR.Get(key.GetName())
        copy_dir(subdir)

    for channel in CHANNELS:
        h_asimov_SR = construct_asimov(channel, template_output_SR)
        h_asimov_SR.SetName("nominal")
        template_output_SR.cd(channel+"/data_obs")
        ROOT.gDirectory.Delete("nominal;1")
        h_asimov_SR.Write()
        h_asimov_SR.Delete()

        h_asimov_CR = construct_asimov(channel, template_output_CR)

        # Inject instrumental pT-miss in the CR template
        template_output_SR.cd(channel)
        h_instrMET = ROOT.gDirectory.Get("DYJets/nominal")
        h_mean_weights = mean_weights_file.Get("mean_weights_tot_"+channel)
        print("channel " + channel + ": h_instrMET has " + str(h_instrMET.GetNbinsX()) + " bins and h_mean_weights has " + str(h_mean_weights.GetNbinsX()) + " bins.")
        h_instrMET.Divide(h_mean_weights)
        template_output_CR.cd(channel)
        h_asimov_CR += h_instrMET
        h_instrMET.Delete()
        h_mean_weights.Delete()

        h_asimov_CR.SetName("nominal")
        template_output_CR.cd(channel+"/data_obs")
        ROOT.gDirectory.Delete("nominal;1")
        h_asimov_CR.Write()
        h_asimov_CR.Delete()

    template_output_SR.Close()
    template_output_CR.Close()
