import ROOT as r
import os
import plottery as ply

#def gethisto(name):
#    file = r.TFile.Open('outputs/'+name)
#    hdata = file.h_mttbar
#    print hdata
#    return hdata
#filenames= {'plots_datamuon.root','plots_dataelectron.root','plots_wjets2500toInf.root','plots_ttjets.root','plots_wjets600to800.root','plots_wjets1200to2500.root','plots_wjets800to1200.root','plots_rsgluon2TeV.root','plots_rsgluon3TeV.root'}
files = [
'Data',
'TTJets_DiLept',
'TTWJetsToLNu',
'TTZToLLNuNu_M-10',
'ST_s-channel_4f_leptonDecays',
'ST_t-channel_antitop_4f_inclusiveDecays',
'ST_t-channel_top_4f_inclusiveDecays',
'ST_tW_antitop_5f_inclusiveDecays',
'ST_tW_top_5f_inclusiveDecays',
'DYJetsToLL_M-50',
'DYJetsToLL_M-10to50',
'DYJetsToTauTau_M-50',
'DYJetsToTauTau_M-10to50',
'WWTo2L2Nu',
'WZTo2L2Q',
'WZTo3LNu',
#'WWW_4F',
'WWZ',
'WZZ',
'ZZZ',
#'WJetsToLNu',
'WJetsToLNu_HT-100To200',
'WJetsToLNu_HT-200To400',
'WJetsToLNu_HT-400To600',
'WJetsToLNu_HT-600To800',
'WJetsToLNu_HT-800To1200',
'WJetsToLNu_HT-1200To2500',
'WJetsToLNu_HT-2500ToInf',
'ZZTo2L2Nu',
'ZZTo2L2Q',
'ZZToTauTau2Nu',
'ZZToTauTau2Q',
'ZZTo4L'
] # input Root file names
xsec = {
'TTJets_DiLept':87.31,
'TTWJetsToLNu':0.2043,
'TTZToLLNuNu_M-10':0.2529,
'WWTo2L2Nu':12.178,
'WWToLNuQQ':49.997,
'WJetsToLNu':61526.7,
'WJetsToLNu_HT-100To200':1345.,
'WJetsToLNu_HT-200To400':359.7,
'WJetsToLNu_HT-400To600':48.91,
'WJetsToLNu_HT-600To800':12.05,
'WJetsToLNu_HT-800To1200':5.501,
'WJetsToLNu_HT-1200To2500':1.329,
'WJetsToLNu_HT-2500ToInf':0.03216,
'DYJetsToLL_M-50':5765,
'DYJetsToLL_M-10to50':18610.,
'DYJetsToTauTau_M-50':5765,
'DYJetsToTauTau_M-10to50':18610.,
'ST_s-channel_4f_leptonDecays':3.362,
'ST_t-channel_antitop_4f_inclusiveDecays':70.69,
'ST_t-channel_top_4f_inclusiveDecays':70.69,
'ZZTo2L2Nu':0.564,
'ZZTo2L2Q':3.22,
'ZZToTauTau2Nu':0.564,
'ZZToTauTau2Q':3.22,
'ZZTo4L':1.256,
'ST_tW_antitop_5f_inclusiveDecays':35.6,
'ST_tW_top_5f_inclusiveDecays':35.6,
'WZTo2L2Q':5.595,
'WZTo3LNu':4.42965,
'ZZZ':0.01398 ,
'WZZ':0.05565 ,
'WWZ':0.16510 
#'WWW_4F':0.1651 
}
allprocess={
'Data':'Data',
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
'mt_Outbveto125_tot',
'mt_Outbtag125_tot',
'mt_Outbtag50_tot',
'mt_Inbveto125_tot',
'mt_Inbveto80_tot',
'mt_Inbveto50_tot',
'mt_Inbtag50_tot',
'met_Outbveto_tot',
'met_Outbtag_tot',
'mt_Outbtag80_tot',
'zmass_bveto125_tot',
'mt_Inbtag125_tot',
'eventflow_tot',
'zmass_bveto50_tot',
'zmass_btag80_tot',
'mt_Outbveto50_tot',
'zmass_btag50_tot',
'zmass_btag50_tot',
'zmass_bveto80_tot',
'mt_Inbtag80_tot',
'zmass_btag125_tot',
'mt_Outbveto80_tot',
'met_Inbtag_tot',
'met_Inbveto_tot'
]
channel = ['ee','mumu','emu','ll']
instLumi = 35920.0


def writeHisto(filename,isMC):
    global histos
    for ch in channel:
        for histo in histos:
            file = r.TFile.Open("../OUTPUTS/firstTest_NRB/MERGED/outputNRB_"+filename+".root")
            h_Data =r.TH1F()
            pointer = file.FindObjectAny(histo+'_'+ch)
            if not pointer == None:
                file.GetObject(histo+'_'+ch,h_Data)
                h_Data = h_Data.DrawCopy()
                if isMC:
                    Nevent = r.TH1F()
                    file.GetObject("totEventInBaobab_tot",Nevent)
                    norm = instLumi*xsec[filename]/Nevent.Integral();
                    print filename," norm is:",norm
                    h_Data.Scale(norm)
                    del Nevent
                h_Data.SetName(histo+'_'+ch+'_'+filename)
                fff=r.TFile.Open("normalized.root","update")
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
    global file,histoname,allprocess
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

if not os.path.isfile('normalized.root'):
    f_tmp = r.TFile("normalized.root","RECREATE")
    h111= r.TH1F()
    f_tmp.Close()
    for i in range(len(files)):
        if i ==0:
            writeHisto(files[i],False)
        else:
            writeHisto(files[i],True)
for ch in channel:
    for histo in histos:
        histoname = histo+'_'+ch
        file = r.TFile.Open("normalized.root")
        #h_Data = getHisto(file,'Data',histoname)
        #h_Top  = getHisto(file,'TTJets_DiLept',histoname)
        #file.GetObject(histoname+'_TTJets_DiLept',h_Top)
        #for i in range(2,8):
        #    addHisto(file,files[i],histoname,h_Top)
        #h_DY = r.TH1F()
        #file.GetObject(histoname+'_DYJetsToLL_M-50',h_DY)
        #addHisto(file,files[10],histoname,h_DY)
        h_Data = finalHisto("Data")
        h_Ztt   = finalHisto("Z#rightarrow #tau#tau")
        h_Zll   = finalHisto("Z#rightarrow ee/#mu#mu")
        h_W    = finalHisto("W#rightarrow l#nu")
        h_ZZ   = finalHisto("ZZ")
        h_ZZtt   = finalHisto("ZZ#rightarrow Z#tau#tau")
        h_WZ   = finalHisto("WZ")
        h_ZVV  = finalHisto("ZVV")      
        h_Top  = finalHisto("Top")
        h_WW   = finalHisto('WW')
        allproc = ['Z#rightarrow ee/#mu#mu','Z#rightarrow #tau#tau','W#rightarrow l#nu','ZZ','ZZ#rightarrow Z#tau#tau','WZ','ZVV','Top','WW']
        allcolor = [831,833,809,592,595,594,869,8,590]
        mylist = [h_Zll,h_Ztt,h_W,h_ZZ,h_ZZtt,h_WZ,h_ZVV,h_Top,h_WW]
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
            colors = allcolor,
            legend_labels = allproc,
            options = {
                "do_stack": True,
                "legend_scalex": 0.7,
                "legend_scaley": 1.8,
               # "extra_text": ["#slash{E}_{T} > 50 GeV","N_{jets} #geq 2","H_{T} > 300 GeV"],
                "yaxis_log": True,
                "ratio_range":[0.4,1.6],
                # "ratio_pull": True,
                #"hist_disable_xerrors": True,
                #"ratio_chi2prob": True,
                "output_name": histoname+".pdf",
                #"legend_percentageinbox": True,
                "yaxis_moreloglabels":False,
                "yaxis_range": [0.01,100000000],
                "cms_label": "Preliminary",
                "lumi_value": "35.9",
                "legend_border": True,
                "legend_smart": False,
                "hist_line_black":True,
                #"bin_text_size": 1.2
                #"bkg_sort_method":"unsorted",
                "legend_percentageinbox":False,
                #"output_ic": True,
                "xaxis_label": "(GeV)",
                #"us_flag": True,
                # "output_jsroot": True,
                # "output_diff_previous": True,
                }
            )
    