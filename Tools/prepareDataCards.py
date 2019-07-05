#!/usr/bin/env python
from __future__ import print_function



import sys
import re
import os
import argparse
import shutil
import itertools
import ROOT
import numpy as np
import errno

import yaml

processesDictionnary={
#structure =
# name of the process : [list of the samples contributing to the given process]
"Data": {"samples":["Data"],},
#"ggH":  ["ggHname"],
#"qqH":  ["qqHname"],
"InstrMET":    {"samples":["DYJetsToLL_M-50"],"processID":4},
"ZZ":   {"samples":["ZZTo4L","ZZTo2L2Nu","ZZTo2L2Q"],"processID":1},
"WZ":   {"samples":["WZTo2L2Q","WZTo3LNu"],"processID":2},
"TopWW":   {"samples":["TTJets_DiLept","WWTo2L2Nu"],"processID":3},
#"TopWW":   {"samples":["TT_TuneCUETP8M2T4","WWTo2L2Nu","ST_s-channel_4f_leptonDecays","ST_t-channel_top_4f_inclusiveDecays","ST_t-channel_antitop_4f_inclusiveDecays","ST_tW_top_5f_inclusiveDecays","ST_tW_antitop_5f_inclusiveDecays","TTWJetsToLNu","TTZToLLNuNu_M-10","WJetsToLNu"],"processID":3},
"ggHsigOnly": {"samples":["GluGluHToZZTo2L2Nu_M800"],"processID":0}#signal should have a process ID =0 or negative (convention in combine)
}

sampleInfosDictionary={}
systInfoDictionary={}
dataDrivenMode=False
integratedLumi = 0
pathToHistos=""
outputPath=""
contentDictionnary={}
#observable/variable used to extract the limit
observable="mT_final"

jetCategories=["eq0jets","geq1jets","vbf"]
leptonsCategories=["mumu","ee"]
#file containing the samples infos
sampleFile="samples.h"
#file containing the syst description
systFile="config/syst.yaml"
listRankedSample=["ZZ","WZ","TopWW","InstrMET","Total Bkgd.","Data"]

#The leading genuine MET samples used in the Instr.MET building and considered for systematics. Names follow conventions from the Tools/harvestInstrMET.C script
genuineMETsamplesInInstrMET=["WJets","WGamma","ZGamma"]

#Normalization uncertainties. FIXME: just done for yields, not for datacards for now.
normalization_uncertainty = {}
normalization_uncertainty["ZZ"] = 0.00
normalization_uncertainty["WZ"] = 0.00
normalization_uncertainty["TopWW"] = 0.00
normalization_uncertainty["InstrMET"] = 0.10
unc_lumi = 0.026
unc_eff_e = 0.072124 #Comes from old fwk, not checked or recomputed.
unc_eff_mu = 0.061788 #Comes from old fwk, not checked or recomputed.

def parse_command_line():
    """Parse the command line parser"""
    parser = argparse.ArgumentParser(description='Launch baobab nutple production.')

    parser.add_argument('--suffix', action='store', default=None,
                        help='suffix that will be added to the output directory')
    parser.add_argument('--dataDriven', action='store_true', default=None,
                        help='will use the data driven background estimation when possible')

    return parser.parse_args()

def load_the_samples_dict(base_path):
    global sampleInfosDictionary
    global integratedLumi
    try:
      sampleInfoFile = open(base_path+"/"+sampleFile,'r')
    except IOError:
        raise NameError("\033[1;31m "+base_path+"/"+sampleFile+" not found !\n\033[0;m")
    samplelines = sampleInfoFile.readlines()
    for sample in samplelines:
        if "instLumi" in sample:
            integratedLumi=float(sample.split("=")[1][:-2])
        if not ("allMCsamples.push_back" in sample or "GluGlu" in sample or "VBF" in sample):
            continue
        crossSectionString=sample.split(",")[-3]
        crossSection = 0
        if not "*" in crossSectionString:
            crossSection = float(crossSectionString)
        else:
            crossSection=float(crossSectionString.split("*")[0])*float(crossSectionString.split("*")[1])
        sampleInfosDictionary[sample.split(",")[-4].replace(" ","")[1:-1]] = crossSection

def load_systs_list(base_path):
    global systInfoDictionary

    with open(os.path.join(base_path, systFile)) as f:
        syst_mapping = yaml.safe_load(f)

    for syst_stem, dataset_masks in syst_mapping.items():
        for direction in ['up', 'down']:
            key = '{}_{}'.format(syst_stem, direction)
            systInfoDictionary[key] = dataset_masks

        if 'MC' in dataset_masks:
            for sample, direction in itertools.product(
                genuineMETsamplesInInstrMET, ['up', 'down']
            ):
                key = '{}_{}_{}'.format(sample, syst_stem, direction)
                systInfoDictionary[key] = ['InstrMET']


def load_a_sample(categories,sample):
    subsamples = processesDictionnary[sample]["samples"]
    histosData={}
    histoCentral={}
    for aSubSample in subsamples:
        fdata = ROOT.TFile(pathToHistos+"outputHZZ_"+aSubSample+".root")
        ROOT.gROOT.cd()
        nbEventsInBaobabs = fdata.Get("totEventInBaobab_tot").Integral()
        sampleWeight = 1
        #if not aSubSample=="Data" and not aSubSample=="InstrMET" and not "GluGluH" in aSubSample: #Why not GluGluH???
        if not aSubSample=="Data" and not aSubSample=="InstrMET":
            sampleWeight = integratedLumi*sampleInfosDictionary[aSubSample]/nbEventsInBaobabs
        for aCategory in categories:
            histoName=observable+"_"+aCategory[1]+"_"+aCategory[0]
            histo = fdata.Get(histoName)
            if not histo:
                raise NameError("\033[1;31m "+histoName+ " not found in "+pathToHistos+"outputHZZ_"+aSubSample+".root\033[0;m")
            localHisto = histo.Clone(histoName+"_"+aSubSample)#need clone the histo otherwise it will be overriden
            localHisto.Scale(sampleWeight)
            if not aCategory in histoCentral.keys():
                histoCentral[aCategory]= localHisto
            else:
                histoCentral[aCategory].Add(localHisto)
            #print(aSubSample+" "+aCategory[1]+aCategory[0]+" "+str(histoCentral[aCategory].Integral()))

    histosData["nominal"] = histoCentral
    for aSyst in systInfoDictionary:
        #print("will try the syst "+aSyst)
        histoSyst={}
        for aSubSample in subsamples:
            doSyst=False
            if aSubSample in systInfoDictionary[aSyst] or ("MC" in systInfoDictionary[aSyst] and not sample=="Data" and not (sample=="InstrMET" and dataDrivenMode)):
                #print("the syst "+aSyst+" will be considered for the sample "+aSubSample)
                doSyst=True
            fdata = ROOT.TFile(pathToHistos+"outputHZZ_"+aSubSample+".root")
            ROOT.gROOT.cd()
            nbEventsInBaobabs = fdata.Get("totEventInBaobab_tot").Integral()
            sampleWeight = 1
            #if not aSubSample=="Data" and not aSubSample=="InstrMET" and not "GluGluH" in aSubSample: #Why not GluGluH???
            if not aSubSample=="Data" and not aSubSample=="InstrMET":
                sampleWeight = integratedLumi*sampleInfosDictionary[aSubSample]/nbEventsInBaobabs
            for aCategory in categories:
                if doSyst:
                    histoName=observable+"_"+aCategory[1]+"_"+aCategory[0]+"_"+aSyst
                else:
                    histoName=observable+"_"+aCategory[1]+"_"+aCategory[0]
                histo = fdata.Get(histoName)
                if not histo:
                    raise NameError("\033[1;31m "+histoName+ " not found in "+pathToHistos+"outputHZZ_"+aSubSample+".root\033[0;m")
                localHisto = histo.Clone(histoName+"_"+aSubSample)#need clone the histo otherwise it will be overriden
                localHisto.Scale(sampleWeight)
                if not aCategory in histoSyst.keys():
                    histoSyst[aCategory]= localHisto
                else:
                    histoSyst[aCategory].Add(localHisto)

        histosData[aSyst] = histoSyst
    contentDictionnary[sample]=histosData

def load_all_samples(categories):
    for theSample in processesDictionnary:
        print("loading "+theSample)
        load_a_sample(categories,theSample)

def sum_event_allCategories(categories, histosSample):
    totalSum = 0
    for aCat in categories:
        totalSum += histosSample[aCat].Integral()
    return totalSum

def sum_event_muons(categories, histosSample):
    totalSum = 0
    for aCat in categories:
        if (aCat[0] == "mumu"): totalSum += histosSample[aCat].Integral()
    return totalSum

def sum_event_electrons(categories, histosSample):
    totalSum = 0
    for aCat in categories:
        if (aCat[0] == "ee"): totalSum += histosSample[aCat].Integral()
    return totalSum

def square_inclusive_syst(categories, h_nominal, h_variation):
    square_inclusive_syst = 0.
    for aCat in categories: 
        variation = h_nominal[aCat].Integral()-h_variation[aCat].Integral()
        square_inclusive_syst += variation*variation
    return square_inclusive_syst

def stat_error_from_TH1(theInputHisto):
    nbBins = theInputHisto.GetNbinsX()
    sumSqrError=0
    for iBin in range(0,nbBins):
        theError = theInputHisto.GetBinError(iBin+1)
        sumSqrError += theError*theError
    return np.sqrt(sumSqrError)

def sum_stat_error_allCategories(categories, histosSample):
    totalSum = 0
    for aCat in categories:
        totalSum += np.power(stat_error_from_TH1(histosSample[aCat]), 2)
    return np.sqrt(totalSum)

def add_global_normalization_uncertainties_on_the_inclusive_bin(categories, squaredErrorOnSum, aSample, totalSumForThisSample):
    if dataDrivenMode and (aSample=="InstrMET" or aSample=="Total Bkgd."):
        sample="InstrMET"
        totalSumForThisSample = sum_event_allCategories(categories, contentDictionnary[sample]["nominal"])
        squaredErrorOnSum += np.power(totalSumForThisSample*normalization_uncertainty[sample], 2)
    if not (dataDrivenMode and aSample=="InstrMET"): 
        squaredErrorOnSum += np.power(totalSumForThisSample*unc_lumi, 2)
        squaredErrorOnSum += np.power(sum_event_electrons(categories, contentDictionnary[aSample]["nominal"])*unc_eff_e, 2)
        squaredErrorOnSum += np.power(sum_event_muons(categories, contentDictionnary[aSample]["nominal"])*unc_eff_mu, 2)
    
    return squaredErrorOnSum

def add_global_normalization_uncertainties_on_the_sample(aSample, squaredError, aCategory):
    if aSample=="Total Bkgd.":
        for sample in listRankedSample:
          if (sample =="Total Bkgd.") or (sample=="Data"): continue
          squaredError += np.power(contentDictionnary[sample]["nominal"][aCategory].Integral()*normalization_uncertainty[sample], 2)
          if not (dataDrivenMode and aSample=="InstrMET"): squaredError += np.power(contentDictionnary[sample]["nominal"][aCategory].Integral()*unc_lumi, 2)
          if (not (dataDrivenMode and aSample=="InstrMET") and aCategory[0] == "ee"): squaredError += np.power(contentDictionnary[sample]["nominal"][aCategory].Integral()*unc_eff_e, 2)
          if (not (dataDrivenMode and aSample=="InstrMET") and aCategory[0] == "mumu"): squaredError += np.power(contentDictionnary[sample]["nominal"][aCategory].Integral()*unc_eff_mu, 2)
    else:
        squaredError += np.power(contentDictionnary[aSample]["nominal"][aCategory].Integral()*normalization_uncertainty[aSample], 2)
        if not (dataDrivenMode and aSample=="InstrMET"): squaredError += np.power(contentDictionnary[aSample]["nominal"][aCategory].Integral()*unc_lumi, 2)
        if (not (dataDrivenMode and aSample=="InstrMET") and aCategory[0] == "ee"): squaredError += np.power(contentDictionnary[aSample]["nominal"][aCategory].Integral()*unc_eff_e, 2)
        if (not (dataDrivenMode and aSample=="InstrMET") and aCategory[0] == "mumu"): squaredError += np.power(contentDictionnary[aSample]["nominal"][aCategory].Integral()*unc_eff_mu, 2)
    return squaredError

def web_latex_header():
    header  = "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 3.2 Final//EN\">\n"
    header += "<html>\n"
    header += "<head>\n"
    header += "<title>Yields</title>\n"
    header += "<script type=\"text/x-mathjax-config\">\n"
    header += "  MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});\n"
    header += "</script>\n"
    header += "<script type=\"text/javascript\"\n"
    header += "  src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML\">\n"
    header += "</script>\n"
    header += "</head>\n"
    header += "\n"
    header += "<body>\n"
    header += "<h2>Yields</h2>\n"
    header += "\n"

    return header

def web_latex_footer():
    footer  ="\n"
    footer += "</body>\n"
    footer += "</html>"

def print_number_table(categories):
    tableFile = open(outputPath+"yields.tex","w")

    tableFile.write("\\begin{array}{|c|c|c|c|c|c|c|} \n")
    tableFile.write("\\hline \n")
    tableFile.write("\\text{Channel} & \\text{Inclusive} ")
    for aCategory in categories:
        if(aCategory[1] == "eq0jets"): jetCatText = " =\\text{0jets} "
        elif(aCategory[1] == "geq1jets"): jetCatText = " \\geq\\text{1jets} "
        elif(aCategory[1] == "vbf"): jetCatText = " \\text{VBF} "

        if(aCategory[0] == "ee"): lepCatText = " ee "
        elif(aCategory[0] == "mumu"): lepCatText = " \mu\mu "
        tableFile.write(" & "+lepCatText+"\\  "+jetCatText)
    tableFile.write(" \\\\ \n")
    tableFile.write("\\hline \n")
    for aSample in listRankedSample:
        listSystForThisSample=list(contentDictionnary[aSample].keys())
        listSystForThisSample.remove("nominal")
        listTypeOfSyst=[]
        for aSyst in listSystForThisSample:
            systShort=aSyst.rsplit("_", 1)[0]
            if not systShort in listTypeOfSyst:
                listTypeOfSyst.append(systShort)
        if aSample in "Total Bkgd.": tableFile.write("\\hline \n")
        tableFile.write(" \\text{"+ aSample +"} ")
        if not aSample in processesDictionnary.keys() and not "Total Bkgd.":
            raise NameError("\033[1;31m "+aSample+" not in the list of samples: check the listRankedSample list ! \033[0;m")
        histosSample=contentDictionnary[aSample]["nominal"]
        totalSumForThisSample = sum_event_allCategories(categories, histosSample)
        totalStatErrorForThisSample = sum_stat_error_allCategories(categories, histosSample)
        #error on the sum
        squaredErrorOnSumUp=0
        squaredErrorOnSumDown=0
        for aSyst in listTypeOfSyst:
            #UP variation
            histosSampleSystUp = contentDictionnary[aSample][aSyst+"_up"]
            histosSampleSystDown = contentDictionnary[aSample][aSyst+"_down"]
            totalSumForThisSampleAndSystUp = sum_event_allCategories(categories, histosSampleSystUp)
            upVariation = sum_event_allCategories(categories, histosSampleSystUp)-totalSumForThisSample
            downVariation = sum_event_allCategories(categories, histosSampleSystDown)-totalSumForThisSample
            if upVariation<0 and downVariation>=0:
                squaredErrorOnSumUp+= square_inclusive_syst(categories, histosSample, histosSampleSystDown)
                squaredErrorOnSumDown+= square_inclusive_syst(categories, histosSample, histosSampleSystUp)
                #print("Sample "+aSample+", syst "+aSyst+": nominal = "+str(totalSumForThisSample)+". Uncertainty up = "+str(np.sqrt(square_inclusive_syst(categories, histosSample, histosSampleSystDown)))+" and uncertainty down = "+str(np.sqrt(square_inclusive_syst(categories, histosSample, histosSampleSystUp)))+". Anticorrelated.")
            elif (upVariation<0 and downVariation<0) or (upVariation>0 and downVariation>0):
                raise NameError("\033[1;31m Both the up and down variations are in the same direction for sample "+aSample+" and syst "+aSyst+". Stopping here. \033[0;m")
            else:
                squaredErrorOnSumUp+= square_inclusive_syst(categories, histosSample, histosSampleSystUp)
                squaredErrorOnSumDown+= square_inclusive_syst(categories, histosSample, histosSampleSystDown)
                #print("Sample "+aSample+", syst "+aSyst+": nominal = "+str(totalSumForThisSample)+". Uncertainty up = "+str(np.sqrt(square_inclusive_syst(categories, histosSample, histosSampleSystUp)))+" and uncertainty down = "+str(np.sqrt(square_inclusive_syst(categories, histosSample, histosSampleSystDown)))+".")
        if aSample=="Data": #don't print syst for Data
            #tableFile.write(" & "+str(round(totalSumForThisSample,2))+" \pm "+str(round(totalStatErrorForThisSample,2))) # No sense to add an uncertainty to an observed number.
            tableFile.write(" & "+str(int(round(totalSumForThisSample,2))))
        else:
            #Add the uncertainty on the method for Instr.MET. FIXME: This is just done here and not (yet) propagated to datacards or anywhere else!
            squaredErrorOnSumUp = add_global_normalization_uncertainties_on_the_inclusive_bin(categories, squaredErrorOnSumUp, aSample, totalSumForThisSample)
            squaredErrorOnSumDown = add_global_normalization_uncertainties_on_the_inclusive_bin(categories, squaredErrorOnSumDown, aSample, totalSumForThisSample)
            tableFile.write(" & "+str(round(totalSumForThisSample,2))+" \pm "+str(round(totalStatErrorForThisSample,2))+"^{+"+str(round(np.sqrt(squaredErrorOnSumUp),2))+"}_{-"+str(round(np.sqrt(squaredErrorOnSumDown),2))+"}" )

        for aCategory in categories:
            squaredErrorUp=0
            squaredErrorDown=0
            for aSyst in listTypeOfSyst:
                #UP variation
                histosSampleSystUp = contentDictionnary[aSample][aSyst+"_up"][aCategory]
                histosSampleSystDown = contentDictionnary[aSample][aSyst+"_down"][aCategory]
                upVariation=histosSampleSystUp.Integral()-histosSample[aCategory].Integral()
                downVariation=histosSampleSystDown.Integral()-histosSample[aCategory].Integral()
                if upVariation<0 and downVariation>=0:
                    squaredErrorUp+=(downVariation*downVariation)
                    squaredErrorDown+=(upVariation*upVariation)
                    #print("Sample "+aSample+", category = "+str(aCategory)+", syst "+aSyst+": nominal = "+str(histosSample[aCategory].Integral())+". Uncertainty up = "+str(downVariation)+" and uncertainty down = "+str(upVariation)+". Anticorrelated.")
                #elif (upVariation<0 and downVariation<0) or (upVariation>0 and downVariation>0):
                #    raise NameError("\033[1;31m Both the up and down variations are in the same direction for sample "+aSample+" and syst "+aSyst+". Stopping here. \033[0;m")
                else:
                    squaredErrorUp+=(upVariation*upVariation)
                    squaredErrorDown+=(downVariation*downVariation)
                    #print("Sample "+aSample+", category = "+str(aCategory)+", syst "+aSyst+": nominal = "+str(histosSample[aCategory].Integral())+". Uncertainty up = "+str(upVariation)+" and uncertainty down = "+str(downVariation)+".")
            if aSample=="Data":
                #tableFile.write(" & "+str(round(histosSample[aCategory].Integral(),2))+" \pm "+str(round(stat_error_from_TH1(histosSample[aCategory]),2))) # No sense to add an uncertainty to an observed number.
                tableFile.write(" & "+str(int(round(histosSample[aCategory].Integral(),2))))
            else:
                #Add the uncertainty on the method for Instr.MET. FIXME: This is just done here and not (yet) propagated to datacards or anywhere else!
                squaredErrorUp = add_global_normalization_uncertainties_on_the_sample(aSample, squaredErrorUp, aCategory)
                squaredErrorDown = add_global_normalization_uncertainties_on_the_sample(aSample, squaredErrorDown, aCategory)
                tableFile.write(" & "+str(round(histosSample[aCategory].Integral(),2))+" \pm "+str(round(stat_error_from_TH1(histosSample[aCategory]),2))+"^{+"+str(round(np.sqrt(squaredErrorUp),2))+"}_{-"+str(round(np.sqrt(squaredErrorDown),2))+"}" )
        tableFile.write(" \\\\ \n")
    tableFile.write("\\hline \n")
    tableFile.write("\\end{array} \n")

    tableFile.close()
    #View table on the web
    with open(outputPath+"yields.tex") as f:
        with open(outputPath+"webYields.html", "w") as f1:
            f1.write(web_latex_header())
            for line in f:
                    f1.write(line)


def add_the_sum(categories):
    listSamplesName = list(contentDictionnary.keys())
    listSamplesName.remove('Data')
    if 'ggHsigOnly' in listSamplesName: listSamplesName.remove('ggHsigOnly') #need to be done better in the futur ;) 
    allSumData = {}
    sumEntryNominal={}
    for aCategory in categories:
        histoSum = contentDictionnary[listSamplesName[0]]["nominal"][aCategory].Clone("sum_"+aCategory[1]+"_"+aCategory[0])
        for aSample in listSamplesName[1:]:
            histoSum.Add(contentDictionnary[aSample]["nominal"][aCategory])
        sumEntryNominal[aCategory] = histoSum
    allSumData["nominal"] = sumEntryNominal
    for aSyst in systInfoDictionary:
        sumEntrySyst={}
        for aCategory in categories:
            histoSumSyst = allSumData["nominal"][aCategory].Clone("sum_"+aCategory[1]+"_"+aCategory[0]+"_"+aSyst)
            histoSumSyst.Reset()
            for aSample in listSamplesName:
                if aSyst in list(contentDictionnary[aSample].keys()):
                    histoSumSyst.Add(contentDictionnary[aSample][aSyst][aCategory].Clone("sum_"+aCategory[1]+"_"+aCategory[0]+"_"+aSyst))
            sumEntrySyst[aCategory] = histoSumSyst
        allSumData[aSyst] = sumEntrySyst
    contentDictionnary["Total Bkgd."]= allSumData

def has_an_intersection(list1, list2):
    hasInter=False
    for aItem in list1:
        if aItem in list2:
            return True
    return False

def give_normError(histoDict, category,theSyst):
    nominalValue=histoDict["nominal"][category].Integral()
    upValue=histoDict[theSyst+"_up"][category].Integral()
    downValue=histoDict[theSyst+"_down"][category].Integral()
    if nominalValue==0:
        return "-"
    else:
        upRelat = (upValue-nominalValue)/nominalValue
        downRelat = (nominalValue-downValue)/nominalValue
        return str(round(1-downRelat,4)) + "/" + str(round(1+upRelat,4)) #Convention is down/up apparently. Allows anti-correlated uncertainties.
    return ""

def addSpaces(theString):
    outString=""
    for aSpace in range(0, 25-len(theString)):
        outString+=" "
    return outString

def printWithFixedSpacing(theString):
    outString=theString+addSpaces(theString)
    return outString

def create_datacard(categories):
    fWrite = ROOT.TFile(outputPath+"hzz2l2v_X_2016.root","RECREATE")
    allProcButData=list(processesDictionnary.keys())
    allProcButData.remove("Data")
    nbOfProcess = len(allProcButData)

    #get the process IDs
    procIDs = []
    for aProc in allProcButData:
        procIDs.append(processesDictionnary[aProc]["processID"])
    minProc=min(procIDs)
    maxProc=max(procIDs)
    #get the type of syst
    listTypeOfSyst=[]
    for aSyst in systInfoDictionary:
        systShort=aSyst.rsplit("_",1)[0]
        if not systShort in listTypeOfSyst:
            listTypeOfSyst.append(systShort)

    for theCat in categories:
        catName=theCat[0]+theCat[1]
        fWrite.mkdir(catName)
        fWrite.cd(catName)
        #print("creating the datacard for category "+catName)
        cardLines="imax 1 number of bins\n"
        cardLines+="jmax "+str(nbOfProcess-1)+" number of processes minus 1\n"
        cardLines+="kmax * number of nuisance parameters\n"
        cardLines+="------------------------\n"
        cardLines+="shapes * * hzz2l2v_X_2016.root "+catName+"/$PROCESS "+catName+"/$PROCESS_$SYSTEMATIC\n"
        cardLines+="------------------------\n"
        cardLines+="bin bin1 \n"
        cardLines+="observation "+str(contentDictionnary["Data"]["nominal"][theCat].Integral())+"\n"
        cardLines+="------------------------\n"
        cardLines+=printWithFixedSpacing("bin")+printWithFixedSpacing("")
        contentDictionnary["Data"]["nominal"][theCat].Write("data_obs")
        for aProc in allProcButData:
            cardLines+=printWithFixedSpacing("bin1")
        cardLines+="\n"
        cardLines+=printWithFixedSpacing("process")+printWithFixedSpacing("")
        for aProcIte in range(minProc,maxProc+1):
            cardLines+=printWithFixedSpacing(allProcButData[procIDs.index(aProcIte)])
        cardLines+="\n"
        cardLines+=printWithFixedSpacing("process")+printWithFixedSpacing("")
        for aProcIte in range(minProc,maxProc+1):
            cardLines+=printWithFixedSpacing(str(aProcIte))
        cardLines+="\n"
        cardLines+=printWithFixedSpacing("rate")+printWithFixedSpacing("")
        for aProcIte in range(minProc,maxProc+1):
            cardLines+=printWithFixedSpacing(str(round(contentDictionnary[allProcButData[procIDs.index(aProcIte)]]["nominal"][theCat].Integral(),6)))
            contentDictionnary[allProcButData[procIDs.index(aProcIte)]]["nominal"][theCat].Write(allProcButData[procIDs.index(aProcIte)])
        cardLines+="\n"
        cardLines+="------------------------\n"
        for aSyst in listTypeOfSyst:
            cardLines+=printWithFixedSpacing(aSyst)+printWithFixedSpacing("shape") #Previously shapeN2, switched to "shape" because more "standard".
            for aProcIte in range(minProc,maxProc+1):
                if has_an_intersection(processesDictionnary[allProcButData[procIDs.index(aProcIte)]]["samples"],systInfoDictionary[aSyst+"_down"]) or 'MC' in systInfoDictionary[aSyst+"_down"]:
                    contentDictionnary[allProcButData[procIDs.index(aProcIte)]][aSyst+"_up"][theCat].Write(allProcButData[procIDs.index(aProcIte)]+"_"+aSyst+"_shapeUp")
                    contentDictionnary[allProcButData[procIDs.index(aProcIte)]][aSyst+"_down"][theCat].Write(allProcButData[procIDs.index(aProcIte)]+"_"+aSyst+"_shapeDown")
                    cardLines+=printWithFixedSpacing("1.0")
                else:
                    cardLines+=printWithFixedSpacing("-")
            # This is a relic from before. Actually, it was a mistake: these norm uncertainties are taken by default when we take the "shape" uncertainties.
            #cardLines+="\n"
            #cardLines+=printWithFixedSpacing(aSyst+"_norm")+printWithFixedSpacing("lnN")
            #for aProcIte in range(minProc,maxProc+1):
            #    if has_an_intersection(processesDictionnary[allProcButData[procIDs.index(aProcIte)]]["samples"],systInfoDictionary[aSyst+"_down"]) or 'MC' in systInfoDictionary[aSyst+"_down"]:
            #        cardLines+=printWithFixedSpacing(give_normError(contentDictionnary[allProcButData[procIDs.index(aProcIte)]],theCat,aSyst))
            #    else:
            #        cardLines+=printWithFixedSpacing("-")
            cardLines+="\n"
        #Lumi (by hand)
        cardLines+=printWithFixedSpacing("lumi_13TeV")+printWithFixedSpacing("lnN")
        for aProcIte in range(minProc,maxProc+1):
            if not ("InstrMET" in allProcButData[procIDs.index(aProcIte)]):
                cardLines+=printWithFixedSpacing(str(unc_lumi+1.))
            else:
                cardLines+=printWithFixedSpacing("-")
        cardLines+="\n"
        #Lepton efficiencies (by hand)
        cardLines+=printWithFixedSpacing("CMS_eff_e")+printWithFixedSpacing("lnN")
        for aProcIte in range(minProc,maxProc+1):
            if ('ee' in catName) and not ("InstrMET" in allProcButData[procIDs.index(aProcIte)]): #Data-driven backgrounds do not need to apply the incertainties on lepton id+trigger.
              cardLines+=printWithFixedSpacing(str(unc_eff_e+1.)) 
            else:
              cardLines+=printWithFixedSpacing("-")
        cardLines+="\n"
        cardLines+=printWithFixedSpacing("CMS_eff_mu")+printWithFixedSpacing("lnN")
        for aProcIte in range(minProc,maxProc+1):
            if ('mumu' in catName) and not ("InstrMET" in allProcButData[procIDs.index(aProcIte)]):
                cardLines+=printWithFixedSpacing(str(unc_eff_mu+1.))
            else:
                cardLines+=printWithFixedSpacing("-")
        cardLines+="\n"
        #Instr. MET: uncertainty from closure test
        cardLines+=printWithFixedSpacing("InstrMET_closure")+printWithFixedSpacing("lnN")
        for aProcIte in range(minProc,maxProc+1):
            if 'InstrMET' in allProcButData[procIDs.index(aProcIte)]:
                cardLines+=printWithFixedSpacing(str(1.+normalization_uncertainty["InstrMET"]))
            else:
                cardLines+=printWithFixedSpacing("-")
        cardLines+="\n"
        cardLines+="* autoMCStats 0" #Handles automatically the bin-by-bin stat uncertainties.
        theCard = open(outputPath+"hzz2l2v_X_2016_"+catName+".txt","w")
        theCard.write(cardLines)
        theCard.close()




def main():
    global outputPath
    global pathToHistos
    global dataDrivenMode
    base_path=os.path.expandvars('$HZZ2L2NU_BASE')

    args = parse_command_line()
    try:
        pathToHistos = base_path+"/OUTPUTS/"+args.suffix+"/MERGED/"
    except TypeError:
        print("\033[1;31m you should specify from which output you can create datacarts with the --suffix option\033[0;m")
        raise
    outputPath = base_path+"/OUTPUTS/"+args.suffix+"/PLOTS/YIELDS/"
    if not os.path.exists(os.path.dirname(outputPath)):
        try:
            os.makedirs(os.path.dirname(outputPath))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


    #first the load the sample info
    load_the_samples_dict(base_path)
    #and also the list of systematics
    load_systs_list(base_path)
    print(systInfoDictionary)

    #second built the list of the categories
    #categories = list(itertools.product(jetCategories,leptonsCategories))
    categories = list(itertools.product(leptonsCategories, jetCategories))

    if args.dataDriven:
        print("Instr.MET will be data driven \n")
        processesDictionnary["InstrMET"]["samples"]=["InstrMET"]
        dataDrivenMode=True
    else:
        print("Instr.MET estimation from DY Monte Carlo")


    #third load all the samples
    load_all_samples(categories)
    #print(contentDictionnary)
    #youhou now all the samples are loaded in contentDictionnary !
    #we can use them and for example print the table
    add_the_sum(categories)

    print_number_table(categories)
    #print(contentDictionnary)

    create_datacard(categories)

#    fWrite = ROOT.TFile("histo.root","RECREATE")
#    for aData in histosData:
#        print(aData)
#        histosData[aData].Write("mT_final_"+aData[0]+"_"+aData[1])



if __name__ == '__main__':
    #try:
    main()
    #except KeyboardInterrupt, e:
    #    print("\nBye!")
    #    pass

