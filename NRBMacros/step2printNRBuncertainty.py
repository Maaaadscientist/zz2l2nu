import ROOT as r
import os
import plottery as ply

files = [
'Data',
'TTJets_DiLept',
'TTJets_SingleLeptFromTbar',
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
'WWToLNuQQ',
'WZTo2L2Q',
'WZTo3LNu',
'WWZ',
'WZZ',
'ZZZ',
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

processes={
'Data' :[
          'Data'
        ],
'Top'  :[
          'TTJets_DiLept',
          'TTJets_SingleLeptFromTbar',
          'TTWJetsToLNu',
          'TTZToLLNuNu_M-10',
          'ST_s-channel_4f_leptonDecays',
          'ST_t-channel_antitop_4f_inclusiveDecays',
          'ST_t-channel_top_4f_inclusiveDecays',
          'ST_tW_antitop_5f_inclusiveDecays',
          'ST_tW_top_5f_inclusiveDecays'
        ],
'WW'   :[
          'WWTo2L2Nu',
          'WWToLNuQQ'
        ],
'TopWW':[
          'TTJets_DiLept',
          'TTJets_SingleLeptFromTbar',
          'TTWJetsToLNu','TTZToLLNuNu_M-10',
          'ST_s-channel_4f_leptonDecays',
          'ST_t-channel_antitop_4f_inclusiveDecays',
          'ST_t-channel_top_4f_inclusiveDecays',
          'ST_tW_antitop_5f_inclusiveDecays',
          'ST_tW_top_5f_inclusiveDecays',
          'WWTo2L2Nu',
          'WWToLNuQQ'
        ],
'W'    :[
          'WJetsToLNu_HT-100To200',
          'WJetsToLNu_HT-200To400',
          'WJetsToLNu_HT-400To600',
          'WJetsToLNu_HT-600To800',
          'WJetsToLNu_HT-800To1200',
          'WJetsToLNu_HT-1200To2500',
          'WJetsToLNu_HT-2500ToInf',
        ],
'ttbar':[
          'TTJets_DiLept',
          'TTJets_SingleLeptFromTbar'
        ]
}
histos=[
'mT_final'
]
jet_cats = [
'eq0jets',
'geq1jets',
'vbf'
]
channel = ['ee','mumu','emu']

alphaValue={
'ee':0.361,
'mumu':0.677,
'ee_err':0.006,
'mumu_err':0.009
}
def writeHisto(filename,isMC,histos):
    if not os.path.isfile("../OUTPUTS/NRB/MERGED/outputNRB_Data.root"):
        print 'no such file: "../OUTPUTS/NRB/MERGED/outputNRB_Data.root"'
    for jet_c in jet_cats:
        for ch in channel:
            for histo in histos:
                file = r.TFile.Open("../OUTPUTS/NRB/MERGED/outputNRB_"+filename+".root")
                if not 'mt_shapes_NRBctrl' in histo:
                    h_Data =r.TH1F()
                else:
                    h_Data =r.TH2F()
                pointer = file.FindObjectAny(histo+'_'+jet_c+'_'+ch)
                if not pointer == None:
                    file.GetObject(histo+'_'+jet_c+'_'+ch,h_Data)
                    h_Data = h_Data.Clone()
                    h_Data.SetName(histo+'_'+jet_c+'_'+ch+'_'+filename)
                    fff=r.TFile.Open("forNRBunc.root","update")
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
def finalHisto(file,proc,histoname):
    FirstLoop = True
    h_tmp =r.TH1F()
    for procname in processes[proc]:
        if FirstLoop:
            h_tmp = getHisto(file,procname,histoname)
            if h_tmp is not None:
                FirstLoop = False
        else:
            addHisto(file,procname,histoname,h_tmp)
    return h_tmp
    del h_tmp
def replaceBinContent(histo,ch):
    h_re = histo.Clone()
    for bi in range(1,h_re.GetXaxis().GetNbins()+1):
        val = h_re.GetBinContent(bi)
        err = h_re.GetBinError(bi)
        newval = val * alphaValue[ch]
        newerr = ((err * alphaValue[ch])**2+(val*alphaValue[ch+'_err'])**2) ** 0.5
        h_re.SetBinContent(bi,newval)
        h_re.SetBinError(bi,newerr)
    return h_re
def rescaleMC(histo_Data,histo_MC):
    h_return = histo_MC.Clone()
    h_return.Scale(histo_Data.Integral(1,histo_Data.GetXaxis().GetNbins()+1)/histo_MC.Integral(1,histo_MC.GetXaxis().GetNbins()+1))
    return h_return
if not os.path.isfile('forNRBunc.root'):
    f_tmp = r.TFile("forNRBunc.root","RECREATE")
    h111= r.TH1F()
    f_tmp.Close()
    for i in range(len(files)):
        if i ==0:
            writeHisto(files[i],False,histos)
        else:
            writeHisto(files[i],True,histos)
NormalizedFile = r.TFile.Open("forNRBunc.root")
f_tmp2 = r.TFile("outputHZZ_NRB.root","RECREATE")
h222= r.TH1F()
f_tmp2.Close()

if not os.path.exists('NRB_PLOTS'):
    os.makedirs('NRB_PLOTS')
sum_val_DD = 0.0
sum_stat_errsq_DD = 0.0
sum_syst_errsq_DD =0.0
sum_val_MC = 0.0
sum_syst_errsq_MC =0.0
sum_stat_errsq_MC =0.0
table_DD = ''
table_MC = ''
for ch in ["mumu","ee"]:
    for jet in jet_cats:
        histoname = 'mT_final'
        h_emu = finalHisto(NormalizedFile,'Data',histoname+'_'+jet+'_emu')
        h_Data= replaceBinContent(finalHisto(NormalizedFile,'Data',histoname+'_'+jet+'_emu'),ch)
        h_Data.SetName('mT_final_'+jet+'_'+ch)
        h_Data.SetTitle('mT_final_'+jet+'_'+ch)
        h_Data =h_Data.Clone()
        h_TopWW= finalHisto(NormalizedFile,'TopWW',histoname+'_'+jet+'_'+ch)
        h_Top = finalHisto(NormalizedFile,'Top',histoname+'_'+jet+'_'+ch)
        h_WW = finalHisto(NormalizedFile,'WW',histoname+'_'+jet+'_'+ch)
        valDD_err = r.Double()
        valMC_err = r.Double()
        valDD= h_Data.IntegralAndError(1,h_Data.GetXaxis().GetNbins()+1,valDD_err)
        valMC= replaceBinContent(rescaleMC(h_emu,h_TopWW),ch).IntegralAndError(1,h_Data.GetXaxis().GetNbins()+1,valMC_err)
        print jet,ch,valDD,valDD_err,valDD*0.13
        sum_val_DD+= valDD
        sum_val_MC+= valMC
        sum_stat_errsq_DD += valDD_err**2
        sum_stat_errsq_MC += valMC_err**2
        sum_syst_errsq_DD += (valDD *0.13)**2
        sum_syst_errsq_MC += (valMC *0.13)**2
        table_DD += '& '+'%.2f' % valDD +' \\pm '+'%.2f' %valDD_err +' \\pm '+'%.2f' %(valDD*0.13) +' '
        table_MC += '& '+'%.2f' % valMC +' \\pm '+'%.2f' %valMC_err +' \\pm '+'%.2f' %(valMC*0.13) +' '
        file_tobesaved=r.TFile.Open("outputHZZ_NRB.root","update")
        h_Data.Write()
        ply.plot_hist(
            data=h_Data,
            bgs= [h_TopWW],
            colors = [8],
            legend_labels = ['TopWW'],
            options = {
                "do_stack": True,
                #"legend_scalex": 0.7,
                #"legend_scaley": 1.8,
               # "extra_text": ["#slash{E}_{T} > 50 GeV","N_{jets} #geq 2","H_{T} > 300 GeV"],
                "yaxis_log": True,
                "ratio_range":[0.4,1.6],
                "legend_scalex": 1.0,
                "legend_scaley": 1.8,
                "legend_ncolumns": 2,
                "legend_border": True,
                "legend_coordinates":[0.6, 0.77, 0.93, 0.87],
                "legend_column_separation":0.2,
                # "ratio_pull": True,
                #"hist_disable_xerrors": True,
                #"ratio_chi2prob": True,
                "output_name": "NRB_PLOTS/"+histoname+'_'+jet+'_'+ch+"_TopWW.pdf",
                #"legend_percentageinbox": True,
                "yaxis_moreloglabels":False,
                "yaxis_range": [0.01,10000],
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
        del file_tobesaved
        del h_Data

text_DD = ' \\text{TopWW_DD}  & '+'%.2f' %(sum_val_DD)+ ' \\pm '+'%.2f' %(sum_stat_errsq_DD**0.5) +' \\pm '+'%.2f' %(sum_syst_errsq_DD**0.5)+' '
text_MC = ' \\text{TopWW_MC}e\\mu  & '  +'%.2f' %(sum_val_MC)+ ' \\pm '+'%.2f' %(sum_stat_errsq_MC**0.5) +' \\pm '+'%.2f' %(sum_syst_errsq_MC**0.5)+' '
print text_DD+table_DD+'\\\\'+'\n'+text_MC+table_MC+'\\\\'+'\n'
file_alpha = r.TFile("alphaValue.root","RECREATE")
h_alpha = r.TH1F("alphaValue","alphaValue",2,0,2)
h_alpha.SetBinContent(1,alphaValue['ee'])
h_alpha.SetBinContent(2,alphaValue['mumu'])
h_alpha.SetBinError(1,alphaValue['ee_err'])
h_alpha.SetBinError(2,alphaValue['mumu_err'])
h_alpha.Write()
file_alpha.Close()
