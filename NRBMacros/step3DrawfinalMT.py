import ROOT as r
import os,shutil
import plottery as ply

#def gethisto(name):
#    file = r.TFile.Open('outputs/'+name)
#    hdata = file.h_mttbar
#    print hdata
#    return hdata
#filenames= {'plots_datamuon.root','plots_dataelectron.root','plots_wjets2500toInf.root','plots_ttjets.root','plots_wjets600to800.root','plots_wjets1200to2500.root','plots_wjets800to1200.root','plots_rsgluon2TeV.root','plots_rsgluon3TeV.root'}
files = [
'Data',
'InstrMET',
'NRB',
'GluGluHToZZTo2L2Nu_M800',
'TTJets_DiLept',
'WJetsToLNu_HT-100To200',
'WJetsToLNu_HT-200To400',
'WJetsToLNu_HT-400To600',
'WJetsToLNu_HT-600To800',
'WJetsToLNu_HT-800To1200',
'WJetsToLNu_HT-1200To2500',
'WJetsToLNu_HT-2500ToInf',
'WWTo2L2Nu',
'WZTo2L2Q',
'WZTo3LNu',
'ZZTo2L2Nu',
'ZZTo2L2Q',
'ZZTo4L'
] # input Root file names
allprocess={
'Data':'Data',
'InstrMET':'InstrMET',
'NRB':'TopWW',
'GluGluHToZZTo2L2Nu_M800':'Signal800',
'TTJets_DiLept':'Top',
'TTWJetsToLNu':'Top',
'TTZToLLNuNu_M-10':'Top',
'WWTo2L2Nu':'WW',
#'WJetsToLNu':'W#rightarrow l#nu',
'WJetsToLNu_HT-100To200':'W#rightarrow l#nu',
'WJetsToLNu_HT-200To400':'W#rightarrow l#nu',
'WJetsToLNu_HT-400To600':'W#rightarrow l#nu',
'WJetsToLNu_HT-600To800':'W#rightarrow l#nu',
'WJetsToLNu_HT-800To1200':'W#rightarrow l#nu',
'WJetsToLNu_HT-1200To2500':'W#rightarrow l#nu',
'WJetsToLNu_HT-2500ToInf':'W#rightarrow l#nu',
'DYJetsToLL_M-50':'Z#rightarrow ee/#mu#mu',
'DYJetsToLL_M-10to50':'Z#rightarrow e/#mu#mu',
'DYJetsToTauTau_M-50':'Z#rightarrow #tau#tau',
'DYJetsToTauTau_M-10to50':'Z#rightarrow #tau#tau',
'ST_s-channel_4f_leptonDecays':'Top',
'ST_t-channel_antitop_4f_inclusiveDecays':'Top',
'ST_t-channel_top_4f_inclusiveDecays':'Top',
'ST_tW_antitop_5f_inclusiveDecays':'Top',
'ST_tW_top_5f_inclusiveDecays':'Top',
'WZTo2L2Q':'WZ',
'WZTo3LNu':'WZ',
#'WWW_4F',
'WWZ':'ZVV',
'WZZ':'ZVV',
'ZZZ':'ZVV',
'ZZTo2L2Nu':'ZZ',
'ZZTo2L2Q':'ZZ',
'ZZTo4L':'ZZ',
'ZZToTauTau2Nu':'ZZ#rightarrow Z#tau#tau',
'ZZToTauTau2Q':'ZZ#rightarrow Z#tau#tau'
}

histos=[
'mT_final'
]
jet_cats = ['eq0jets','geq1jets','vbf']
channel = ['ee','mumu']


def writeHisto(filename,isMC):
    histo = 'mT_final'
    for ch in channel:
        for jet in jet_cats:
            if 'WJets' in filename:
                file = r.TFile.Open("../OUTPUTS/NRB/MERGED/outputNRB_"+filename+".root")
            else:
                file = r.TFile.Open("../OUTPUTS/HZZdatadriven/MERGED/outputHZZ_"+filename+".root")
            if not 'mt_shapes_NRBctrl' in histo:
                h_Data =r.TH1F()
            else:
                h_Data =r.TH2F()
            pointer = file.FindObjectAny(histo+'_'+jet+'_'+ch)
            if not pointer == None:
                file.GetObject(histo+'_'+jet+'_'+ch,h_Data)
                h_Data = h_Data.DrawCopy()
                h_Data.SetName(histo+'_'+jet+'_'+ch+'_'+filename)
                fff=r.TFile.Open("normalized_mT.root","update")
                h_Data.Write()
                del fff
            del h_Data



def addHisto(file,process,histoname,prehisto):
    pointer = file.FindObjectAny(histoname+'_'+process)
    if not pointer == None:
        h_tmp =r.TH1F()
        file.GetObject(histoname+'_'+process,h_tmp)
        prehisto.Add(h_tmp)
        del h_tmp
def getHisto(file,process,histoname):
    pointer = file.FindObjectAny(histoname+'_'+process)
    if not pointer == None:
        h_tmp =r.TH1F()
        file.GetObject(histoname+'_'+process,h_tmp)
        return h_tmp
        del h_tmp
    else:
        return None
def finalHisto(proc):
    global histoname,file
    FirstLoop = True
    h_tmp =r.TH1F()
    for procname in allprocess:
        if allprocess[procname] == proc:
            if FirstLoop:
                h_tmp = getHisto(file,procname,histoname)
                if h_tmp is not None:
                    FirstLoop = False
            else:
                addHisto(file,procname,histoname,h_tmp)
    return h_tmp
    del h_tmp


if not os.path.isfile('../OUTPUTS/HZZdatadriven/MERGED/outputHZZ_NRB.root'):
    shutil.copy('outputHZZ_NRB.root','../OUTPUTS/HZZdatadriven/MERGED/outputHZZ_NRB.root')

f_tmp = r.TFile("normalized_mT.root","RECREATE")
h111= r.TH1F()
f_tmp.Close()
for i in range(len(files)):
    if i < 3:
        writeHisto(files[i],False)
    else:
        writeHisto(files[i],True)

if not os.path.exists('NRB_PLOTS'):
    os.makedirs('NRB_PLOTS')
for ch in channel:
    for jet in jet_cats:
        histoname = 'mT_final_'+jet+'_'+ch
        file = r.TFile.Open("normalized_mT.root")
        #h_Data = getHisto(file,'Data',histoname)
        #h_Top  = getHisto(file,'TTJets_DiLept',histoname)
        #file.GetObject(histoname+'_TTJets_DiLept',h_Top)
        #for i in range(2,8):
        #    addHisto(file,files[i],histoname,h_Top)
        #h_DY = r.TH1F()
        #file.GetObject(histoname+'_DYJetsToLL_M-50',h_DY)
        #addHisto(file,files[10],histoname,h_DY)
        h_Data = finalHisto("Data")
        h_Instr= finalHisto("InstrMET")
        h_TopWW= finalHisto("TopWW")
        h_W    = finalHisto("W#rightarrow l#nu")
        h_ZZ   = finalHisto("ZZ")
        h_WZ   = finalHisto("WZ")     
        h_Top  = finalHisto("Top")
        h_WW   = finalHisto('WW')
        h_sig  = finalHisto('Signal800')
        print h_sig
        print h_Data
        allproc = ['W#rightarrow l#nu','WZ','ZZ','TopWW','InstrMET']
        allcolor = [809,594,592,8,831]
        mylist = [h_W,h_WZ,h_ZZ,h_TopWW,h_Instr]
        if None in mylist:
            if type(mylist.index(None)) is not (int):
                for index in mylist.index(None):
                    del allcolor[index]
                    del allproc[index]
            else:
                del allcolor[mylist.index(None)]
                del allproc[mylist.index(None)]
            while None in mylist:
                mylist.remove(None)
        #print allproc,allcolor,mylist
        ply.plot_hist(
            data=h_Data,
            bgs= mylist,
            sigs = [h_sig] ,
            sig_labels = ["signal800"],

            colors = allcolor,
            legend_labels = allproc,
            options = {
                #"canvas_height" :400,
                #"canvas_width" :400,
                "do_stack": True,
                "legend_scalex": 1.5,
                "legend_scaley": 1.8,
                "legend_ncolumns": 3,
                "legend_border": True,
                "legend_coordinates":[0.4, 0.73, 0.93, 0.87],
                "legend_column_separation":0.2,
               # "extra_text": ["#slash{E}_{T} > 50 GeV","N_{jets} #geq 2","H_{T} > 300 GeV"],
                "yaxis_log": True,
                "ratio_range":[0.0,2.0],
                # "ratio_pull": True,
                #"hist_disable_xerrors": True,
                #"ratio_chi2prob": True,
                "output_name": "NRB_PLOTS/"+histoname+"_final.pdf",
                #"legend_percentageinbox": True,
                "yaxis_moreloglabels":False,
                "yaxis_range": [0.01,100000000],
                "cms_label": "Preliminary",
                "lumi_value": "35.9",
                "legend_border": False,
                "legend_rounded":False,
                "legend_smart": False,
                "hist_line_black":True,
                #"bin_text_size": 1.2
                "bkg_sort_method":"unsorted",
                "legend_percentageinbox":False,
                #"output_ic": True,
                "xaxis_label": "Transverse Mass(GeV)",
                "yaxis_label": "Events",
                #"us_flag": True,
                # "output_jsroot": True,
                # "output_diff_previous": True,
                }
            )
    
