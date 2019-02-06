import ROOT as rt
import math
import numpy
import os

def draw_TGraphAsym(graph1,graph2,title,legend):
  graph1.GetYaxis().SetRangeUser(0,3)
  graph1.GetXaxis().SetRangeUser(46,140)
  graph1.SetMarkerColor(9)
  graph1.SetMarkerSize(1)
  graph1.SetMarkerStyle(20)
  graph2.GetYaxis().SetRangeUser(0,3)
  graph2.GetXaxis().SetRangeUser(46,140)
  graph2.SetMarkerColor(46)
  graph2.SetMarkerSize(1)
  graph2.SetMarkerStyle(20)
  legend.AddEntry(graph1,"ee channel("+title+")","p")
  legend.AddEntry(graph2,"mumu channel("+title+")","p")
  graph1.Draw("AP")
  graph2.Draw("P")
  legend.Draw()   
def draw_TGraphAsym_closure(graph1,graph2,legend):
   
  graph1.GetYaxis().SetRangeUser(-0.5,0.5)
  graph1.GetXaxis().SetRangeUser(46,140)
  graph1.SetMarkerColor(9)
  graph1.SetMarkerSize(1)
  graph1.SetMarkerStyle(22)
  graph2.GetYaxis().SetRangeUser(-0.5,0.5)
  graph2.GetXaxis().SetRangeUser(46,140)
  graph2.SetMarkerColor(46)
  graph2.SetMarkerSize(1)
  graph2.SetMarkerStyle(22)
  legend.AddEntry(graph1,"ee channel","p")
  legend.AddEntry(graph2,"mumu channel","p")
  graph1.Draw("AP")
  graph2.Draw("P")
  legend.Draw() 

   
def fillYields(sidebandtype,Btagveto,results_2D = {}):
  for i in range(18):
    for channel in ['ee','mumu']:
      for label in closurelabels:
        Num_allSB_EE       = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)]
        if channel == 'ee' and label == 'WW' and Btagveto == 'bveto' : 
          print 'Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5), Num_allSB_EE  
        Err_allSB_up_EE    = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+label+'_'+channel+'_METcut'+str(50+i*5)]
        Err_allSB_low_EE   = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+label+'_'+channel+'_METcut'+str(50+i*5)]
        Num_allSB_EMu      = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)]
        Err_allSB_up_EMu   = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+label+'_emu_METcut'+str(50+i*5)]
        Err_allSB_low_EMu  = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+label+'_emu_METcut'+str(50+i*5)]
        results_2D['alpha_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)]
        results_2D['Erralpha_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =math.sqrt(Err_allSB_up_EE*numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_EE))/numpy.square(Num_allSB_EMu)
        results_2D['Erralpha_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)]  =math.sqrt(Err_allSB_low_EE*numpy.square(Num_allSB_EMu)+Err_allSB_low_EMu*numpy.square(Num_allSB_EE))/numpy.square(Num_allSB_EMu)
        Num_exp_EE       = results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125']
        Err_exp_up_EE    = math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_up_'+label+'_'+channel+'_METcut125'])
        Err_exp_low_EE   = math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_low_'+label+'_'+channel+'_METcut125'])
        # N pre (or N est) = alpha * N[emu in] = = N[ll out]/N[emu out] *N[emu in]
        Num_pre_EE       = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)] *results_2D['Num_in_bveto_'+label+'_emu_METcut125']
        # E pre  = sqrt (E[emu in]^2 * N[ll in]^2 / N[emu in]^2 + E[ll out]^2 * N[ll out]^2 / N[emu out]^2 + E[emu out]^2 * N[ll out]^2 * N[emu in]^2 / N[emu out]^4 )
        Err_pre_up_EE    = math.sqrt((Err_allSB_up_EE*numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_EE))/numpy.square(numpy.square(Num_allSB_EMu))*results_2D['Num_in_bveto_'+label+'_emu_METcut125']*results_2D['Num_in_bveto_'+label+'_emu_METcut125']+results_2D['Errsquare_in_bveto_up_'+label+'_emu_METcut125']*numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)]))
        #Err_pre_up_EE    = math.sqrt(results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+label+'_'+channel+'_METcut'+str(50+i*5)]  *numpy.square(results_2D['Num_in_bveto_'+label+'_emu_METcut125'] )     /numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)])+
        #                             results_2D['Errsquare_in'+'_'+Btagveto+'_low_'+label+'_emu_METcut125']       *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)])/numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)])+
        #                             results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+label+'_emu_METcut'+str(50+i*5)]*numpy.square(results_2D['Num_in_bveto_'+label+'_emu_METcut125'])     *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)])/numpy.square(numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)])))
        Err_pre_low_EE   = math.sqrt((Err_allSB_low_EE*numpy.square(Num_allSB_EMu)+Err_allSB_low_EMu*numpy.square(Num_allSB_EE))/numpy.square(numpy.square(Num_allSB_EMu))*results_2D['Num_in_bveto_'+label+'_emu_METcut125']*results_2D['Num_in_bveto_'+label+'_emu_METcut125']+results_2D['Errsquare_in_bveto_low_'+label+'_emu_METcut125']*numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)]))
        #Err_pre_low_EE    = math.sqrt(results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+label+'_'+channel+'_METcut'+str(50+i*5)]  *numpy.square(results_2D['Num_in_bveto_'+label+'_emu_METcut125'] )     /numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)])+
         #                            results_2D['Errsquare_in'+'_'+Btagveto+'_low_'+label+'_emu_METcut125']       *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)])/numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)])+
         #                            results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+label+'_emu_METcut'+str(50+i*5)]*numpy.square(results_2D['Num_in_bveto_'+label+'_emu_METcut125'])     *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)])/numpy.square(numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+label+'_emu_METcut'+str(50+i*5)])))
        if not 'Data' in label:
          # N[pre] (or N est) = alpha * N[emu in] = = N[ll out]/N[emu out] *N[emu in]
          results_2D['Num_exp_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125']
          results_2D['ErrNum_exp_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_up_'+label+'_'+channel+'_METcut125'])
          results_2D['ErrNum_exp_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_low_'+label+'_'+channel+'_METcut125'])
          results_2D['Num_pre_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = Num_pre_EE
          results_2D['ErrNum_pre_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = Err_pre_up_EE
          results_2D['ErrNum_pre_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = Err_pre_low_EE
          results_2D['ratio_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =Num_pre_EE/results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125']
          results_2D['ratio_err_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =math.sqrt(numpy.square(Err_pre_up_EE)/numpy.square(results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125'])+numpy.square(math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_up_'+label+'_'+channel+'_METcut125']))*numpy.square(Num_pre_EE)/numpy.square(numpy.square(results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125'])))
          results_2D['ratio_err_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =math.sqrt(numpy.square(Err_pre_low_EE)/numpy.square(results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125'])+numpy.square(math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_low_'+label+'_'+channel+'_METcut125']))*numpy.square(Num_pre_EE)/numpy.square(numpy.square(results_2D['Num_in_bveto_'+label+'_'+channel+'_METcut125'])))
        if 'Data' in label:
          # N[pre] (or N est) = alpha * N[emu in] = = N[ll out]/N[emu out] *N[emu in]
          results_2D['Num_exp_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = results_2D['Num_in_bveto_allMC_'+channel+'_METcut125']
          results_2D['ErrNum_exp_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_up_allMC_'+channel+'_METcut125'])
          results_2D['ErrNum_exp_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_low_allMC_'+channel+'_METcut125'])
          results_2D['Num_pre_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = Num_pre_EE
          results_2D['ErrNum_pre_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = Err_pre_up_EE
          results_2D['ErrNum_pre_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] = Err_pre_low_EE
          results_2D['ratio_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =Num_pre_EE/results_2D['Num_in_bveto_allMC_'+channel+'_METcut125']
          if 'Num_pre_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5) == 'Num_pre_allSB_btag_outputNRB_Data_ee_METcut70':
            print Num_pre_EE, results_2D['Num_in_bveto_allMC_'+channel+'_METcut125'], 'ratio_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)
          results_2D['ratio_err_high_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =math.sqrt(numpy.square(Err_pre_up_EE)/numpy.square(results_2D['Num_in_bveto_allMC_'+channel+'_METcut125'])+numpy.square(math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_up_allMC_'+channel+'_METcut125']))*numpy.square(Num_pre_EE)/numpy.square(numpy.square(results_2D['Num_in_bveto_allMC_'+channel+'_METcut125'])))
          results_2D['ratio_err_low_'+sidebandtype+'_'+Btagveto+'_'+label+'_'+channel+'_METcut'+str(50+i*5)] =math.sqrt(numpy.square(Err_pre_low_EE)/numpy.square(results_2D['Num_in_bveto_allMC_'+channel+'_METcut125'])+numpy.square(math.sqrt(results_2D['Errsquare_in'+'_bveto'+'_low_allMC_'+channel+'_METcut125']))*numpy.square(Num_pre_EE)/numpy.square(numpy.square(results_2D['Num_in_bveto_allMC_'+channel+'_METcut125'])))
def fill_TGraphAsym(graph1,graph2,sidebandtype,Btagveto,DataorMC,results_2D = {}):
  for i in range(18):
  #print i, 50+i*5, results_2D['Num_allSB_EE_ee_METcut'+str(50+i*5)]/results_2D['Num_allSB_EMu_emu_METcut'+str(50+i*5)]
    if 'data' in DataorMC:
      Num_allSB_EE       = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_ee_METcut'+str(50+i*5)]
      Err_allSB_up_EE    = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_outputNRB_Data_ee_METcut'+str(50+i*5)]
      Err_allSB_low_EE   = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_outputNRB_Data_ee_METcut'+str(50+i*5)]
      Num_allSB_MuMu     = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_mumu_METcut'+str(50+i*5)]
      Err_allSB_up_MuMu  = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_outputNRB_Data_mumu_METcut'+str(50+i*5)]
      Err_allSB_low_MuMu = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_outputNRB_Data_mumu_METcut'+str(50+i*5)]
      Num_allSB_EMu      = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_emu_METcut'+str(50+i*5)]
      Err_allSB_up_EMu   = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_outputNRB_Data_emu_METcut'+str(50+i*5)]
      Err_allSB_low_EMu  = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_outputNRB_Data_emu_METcut'+str(50+i*5)]
      graph1.SetPoint(i,50+i*5,results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_ee_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_emu_METcut'+str(50+i*5)]) # alpha = N(ll out)/N(emu out)
      graph1.SetPointEYhigh(i,math.sqrt(Err_allSB_up_EE/numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_EE)/numpy.square(numpy.square(Num_allSB_EMu)))) # err = sqrt (Ell out^2 / Nemu out^2 + Eemu out^2 * Nll out^2/N emu out^4 )
      graph1.SetPointEYlow(i,math.sqrt(Err_allSB_low_EE/numpy.square(Num_allSB_EMu)+Err_allSB_low_EMu*numpy.square(Num_allSB_EE)/numpy.square(numpy.square(Num_allSB_EMu))))
      graph2.SetPoint(i,50+i*5,results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_mumu_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_outputNRB_Data_emu_METcut'+str(50+i*5)])# alpha = N(ll out)/N(emu out)
      graph2.SetPointEYhigh(i,math.sqrt(Err_allSB_up_MuMu/numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_MuMu)/numpy.square(numpy.square(Num_allSB_EMu))))# err = sqrt (Ell out^2 / Nemu out^2 + Eemu out^2 * Nll out^2/N emu out^4 )
      graph2.SetPointEYlow(i,math.sqrt(Err_allSB_low_MuMu/numpy.square(Num_allSB_EMu)+Err_allSB_low_EMu*numpy.square(Num_allSB_MuMu)/numpy.square(numpy.square(Num_allSB_EMu))))
    if not 'data' in DataorMC and not'closure' in DataorMC:
      Num_allSB_EE       = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_ee_METcut'+str(50+i*5)]
      Err_allSB_up_EE    = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+DataorMC+'_ee_METcut'+str(50+i*5)]
      Err_allSB_low_EE   = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+DataorMC+'_ee_METcut'+str(50+i*5)]
      Num_allSB_MuMu     = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_mumu_METcut'+str(50+i*5)]
      Err_allSB_up_MuMu  = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+DataorMC+'_mumu_METcut'+str(50+i*5)]
      Err_allSB_low_MuMu = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+DataorMC+'_mumu_METcut'+str(50+i*5)]
      Num_allSB_EMu      = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_emu_METcut'+str(50+i*5)]
      Err_allSB_up_EMu   = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+DataorMC+'_emu_METcut'+str(50+i*5)]
      Err_allSB_low_EMu  = results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+DataorMC+'_emu_METcut'+str(50+i*5)]
      #print Num_allSB_EMu,'+',Err_allSB_up_EMu,sidebandtype,'\n'
      graph1.SetPoint(i,50+i*5,results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_ee_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_emu_METcut'+str(50+i*5)])# alpha = N[ll out]/N[emu out]
      graph1.SetPointEYhigh(i,math.sqrt(Err_allSB_up_EE*numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_EE))/numpy.square(Num_allSB_EMu))# err = sqrt (E[ll out]^2 / N[emu out]^2 + E[emu out]^2 * N[ll out]^2/N[emu out]^4 )
      graph1.SetPointEYlow(i,math.sqrt(Err_allSB_low_EE*numpy.square(Num_allSB_EMu)+Err_allSB_low_EMu*numpy.square(Num_allSB_EE))/numpy.square(Num_allSB_EMu))
      graph2.SetPoint(i,50+i*5,results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_mumu_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+DataorMC+'_emu_METcut'+str(50+i*5)])# alpha = N[ll out]/N[emu out]
      graph2.SetPointEYhigh(i,math.sqrt(Err_allSB_up_MuMu*numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_MuMu))/numpy.square(Num_allSB_EMu))# err = sqrt (E[ll out]^2 / N[emu out]^2 + E[emu out]^2 * N[ll out]^2/N[emu out]^4 )
      graph2.SetPointEYlow(i,math.sqrt(Err_allSB_low_MuMu*numpy.square(Num_allSB_EMu)+Err_allSB_low_EMu*numpy.square(Num_allSB_MuMu))/numpy.square(Num_allSB_EMu))
      # print math.sqrt(Err_allSB_up_MuMu*numpy.square(Num_allSB_EMu)+Err_allSB_up_EMu*numpy.square(Num_allSB_MuMu)), sidebandtype
    if 'closure' in DataorMC:
      Num_exp_EE       = results_2D['Num_in_bveto_allMC_ee_METcut125']
      Err_exp_up_EE    = math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_ee_METcut125'])
      Err_exp_low_EE   = math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_ee_METcut125'])
      Num_exp_MuMu     = results_2D['Num_in_bveto_allMC_mumu_METcut125']
      Err_exp_up_MuMu  = math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_mumu_METcut125'])
      Err_exp_low_MuMu = math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_mumu_METcut125'])
      # N pre (or N est) = alpha * N[emu in] = = N[ll out]/N[emu out] *N[emu in]
      Num_pre_EE       = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)] *results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125']
      # E pre  = sqrt (E[emu in]^2 * N[ll in]^2 / N[emu in]^2 + E[ll out]^2 * N[ll out]^2 / N[emu out]^2 + E[emu out]^2 * N[ll out]^2 * N[emu in]^2 / N[emu out]^4 )
      Err_pre_up_EE    = math.sqrt(results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)]  *numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'] )     /numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_in'+'_'+Btagveto+'_low_'+'allMC'+'_emu_METcut125']       *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)])/numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+'allMC'+'_emu_METcut'+str(50+i*5)]*numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'])     *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)])/numpy.square(numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])))
      Err_pre_low_EE    = math.sqrt(results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)]  *numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'] )     /numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_in'+'_'+Btagveto+'_low_'+'allMC'+'_emu_METcut125']       *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)])/numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+'allMC'+'_emu_METcut'+str(50+i*5)]*numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'])     *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'ee'+'_METcut'+str(50+i*5)])/numpy.square(numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])))      # N[pre] (or N est) = alpha * N[emu in] = = N[ll out]/N[emu out] *N[emu in]
      Num_pre_MuMu       = results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)]/results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)] *results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125']
      # E pre  = sqrt (E[emu in]^2 * N[ll in]^2 / N[emu in]^2 + E[ll out]^2 * N[ll out]^2 / N[emu out]^2 + E[emu out]^2 * N[ll out]^2 * N[emu in]^2 / N[emu out]^4 )
      Err_pre_up_MuMu    = math.sqrt(results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)]  *numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'] )     /numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_in'+'_'+Btagveto+'_low_'+'allMC'+'_emu_METcut125']       *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)])/numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_up_'+'allMC'+'_emu_METcut'+str(50+i*5)]*numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'])     *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)])/numpy.square(numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])))
      Err_pre_low_MuMu    = math.sqrt(results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)]  *numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'] )     /numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_in'+'_'+Btagveto+'_low_'+'allMC'+'_emu_METcut125']       *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)])/numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])+
                                   results_2D['Errsquare_'+sidebandtype+'_'+Btagveto+'_low_'+'allMC'+'_emu_METcut'+str(50+i*5)]*numpy.square(results_2D['Num_in_bveto_'+'allMC'+'_emu_METcut125'])     *numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_'+'mumu'+'_METcut'+str(50+i*5)])/numpy.square(numpy.square(results_2D['Num_'+sidebandtype+'_'+Btagveto+'_'+'allMC'+'_emu_METcut'+str(50+i*5)])))      # N[pre] (or N est) = alpha * N[emu in] = = N[ll out]/N[emu out] *N[emu in]     print Num_exp_EE,Err_exp_up_EE
      
      graph1.SetPoint(i,50+i*5,Num_pre_EE/Num_exp_EE -1) #Bias = N[pre] / N[exp] -1
      #Err[Bias] = sqrt (E[pre]^2 / N[exp]^2 + E[exp]^2 * N[pre]^2 /N[exp]^4 )
      graph1.SetPointEYhigh(i,math.sqrt(numpy.square(Err_pre_up_EE)/numpy.square(Num_exp_EE)+numpy.square(Err_exp_up_EE)*numpy.square(Num_pre_EE)/numpy.square(numpy.square(Num_exp_EE))))
      graph1.SetPointEYlow(i,math.sqrt(numpy.square(Err_pre_low_EE)/numpy.square(Num_exp_EE)+numpy.square(Err_exp_low_EE)*numpy.square(Num_pre_EE)/numpy.square(numpy.square(Num_exp_EE))))
      graph2.SetPoint(i,50+i*5,Num_pre_MuMu/Num_exp_MuMu -1) #Bias = N[pre] / N[exp] -1
      #Err[Bias] = sqrt (E[pre]^2 / N[exp]^2 + E[exp]^2 * N[pre]^2 /N[exp]^4 )
      graph2.SetPointEYhigh(i,math.sqrt(numpy.square(Err_pre_up_MuMu)/numpy.square(Num_exp_MuMu)+numpy.square(Err_exp_up_MuMu)*numpy.square(Num_pre_MuMu)/numpy.square(numpy.square(Num_exp_MuMu))))
      graph2.SetPointEYlow(i,math.sqrt(numpy.square(Err_pre_low_MuMu)/numpy.square(Num_exp_MuMu)+numpy.square(Err_exp_low_MuMu)*numpy.square(Num_pre_MuMu)/numpy.square(numpy.square(Num_exp_MuMu))))
def fill_TGraphAsym_kMethod(graph1,graph2,Btagveto,closure,results_2D = {}):    
  for i in range(18):
    Nee_peak = results_2D['Num_in_'+Btagveto+'_outputNRB_Data_ee_METcut'+str(50+i*5)]
    Nmumu_peak = results_2D['Num_in_'+Btagveto+'_outputNRB_Data_mumu_METcut'+str(50+i*5)]
    Eee_peak_up = results_2D['Errsquare_in_'+Btagveto+'_up_outputNRB_Data_ee_METcut'+str(50+i*5)]
    Eee_peak_low = results_2D['Errsquare_in_'+Btagveto+'_low_outputNRB_Data_ee_METcut'+str(50+i*5)]
    Emumu_peak_up = results_2D['Errsquare_in_'+Btagveto+'_up_outputNRB_Data_ee_METcut'+str(50+i*5)]
    Emumu_peak_low = results_2D['Errsquare_in_'+Btagveto+'_low_outputNRB_Data_ee_METcut'+str(50+i*5)]
    Kee = 0.5*math.sqrt( Nee_peak /Nmumu_peak) 
    Kmumu = 0.5*math.sqrt( Nmumu_peak /Nee_peak)
    Eee_up = 0.25*math.sqrt(Eee_peak_up /(Nee_peak*Nmumu_peak)+Emumu_peak_up*Nee_peak / (Nmumu_peak*Nmumu_peak*Nmumu_peak))
    Eee_low= 0.25*math.sqrt(Eee_peak_low /(Nee_peak*Nmumu_peak)+Emumu_peak_low*Nee_peak / (Nmumu_peak*Nmumu_peak*Nmumu_peak))
    Emumu_up= 0.25*math.sqrt(Emumu_peak_up /(Nmumu_peak*Nee_peak)+Eee_peak_up*Nmumu_peak / (Nee_peak*Nee_peak*Nee_peak))
    Emumu_low= 0.25*math.sqrt(Emumu_peak_low /(Nmumu_peak*Nee_peak)+Eee_peak_low*Nmumu_peak / (Nee_peak*Nee_peak*Nee_peak))
    graph1.SetPoint(i,50+i*5,Kee) #Bias = N[pre] / N[exp] -1
    #Err[Bias] = sqrt (E[pre]^2 / N[exp]^2 + E[exp]^2 * N[pre]^2 /N[exp]^4 )
    graph1.SetPointEYhigh(i,Eee_up)
    graph1.SetPointEYlow(i,Eee_low)
    graph2.SetPoint(i,50+i*5,Kmumu) #Bias = N[pre] / N[exp] -1
    #Err[Bias] = sqrt (E[pre]^2 / N[exp]^2 + E[exp]^2 * N[pre]^2 /N[exp]^4 )
    graph2.SetPointEYhigh(i,Emumu_up)
    graph2.SetPointEYlow(i,Emumu_low)
    if closure:
      Num_exp_EE       = results_2D['Num_in_bveto_allMC_ee_METcut125']
      Err_exp_up_EE    = math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_ee_METcut125'])
      Err_exp_low_EE   = math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_ee_METcut125'])
      Num_exp_MuMu     = results_2D['Num_in_bveto_allMC_mumu_METcut125']
      Err_exp_up_MuMu  = math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_mumu_METcut125'])
      Err_exp_low_MuMu = math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_mumu_METcut125'])
      Num_pre_EE       = Kee   * results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125'] 
      Num_pre_MuMu     = Kmumu * results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125'] 
      Err_pre_up_EE    = math.sqrt(results_2D['Errsquare_in_bveto_up_outputNRB_Data_emu_METcut125']*Kee*Kee+Eee_up*Eee_up*numpy.square(results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125']))
      Err_pre_low_EE    = math.sqrt(results_2D['Errsquare_in_bveto_low_outputNRB_Data_emu_METcut125']*Kee*Kee+Eee_low*Eee_low*numpy.square(results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125']))
      Err_pre_up_MuMu   = math.sqrt(results_2D['Errsquare_in_bveto_up_outputNRB_Data_emu_METcut125']*Kmumu*Kmumu+Emumu_up*Emumu_up*numpy.square(results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125']))
      Err_pre_low_MuMu    = math.sqrt(results_2D['Errsquare_in_bveto_low_outputNRB_Data_emu_METcut125']*Kmumu*Kmumu+Emumu_low*Emumu_low*numpy.square(results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125']))
      graph1.SetPoint(i,50+i*5,Num_pre_EE/Num_exp_EE -1) #Bias = N[pre] / N[exp] -1
      
      #Err[Bias] = sqrt (E[pre]^2 / N[exp]^2 + E[exp]^2 * N[pre]^2 /N[exp]^4 )
      graph1.SetPointEYhigh(i,math.sqrt(numpy.square(Err_pre_up_EE)/numpy.square(Num_exp_EE)+numpy.square(Err_exp_up_EE)*numpy.square(Num_pre_EE)/numpy.square(numpy.square(Num_exp_EE))))
      graph1.SetPointEYlow(i,math.sqrt(numpy.square(Err_pre_low_EE)/numpy.square(Num_exp_EE)+numpy.square(Err_exp_low_EE)*numpy.square(Num_pre_EE)/numpy.square(numpy.square(Num_exp_EE))))
      graph2.SetPoint(i,50+i*5,Num_pre_MuMu/Num_exp_MuMu -1) #Bias = N[pre] / N[exp] -1
      #print Num_pre_MuMu,Num_exp_MuMu
      #Err[Bias] = sqrt (E[pre]^2 / N[exp]^2 + E[exp]^2 * N[pre]^2 /N[exp]^4 )
      graph2.SetPointEYhigh(i,math.sqrt(numpy.square(Err_pre_up_MuMu)/numpy.square(Num_exp_MuMu)+numpy.square(Err_exp_up_MuMu)*numpy.square(Num_pre_MuMu)/numpy.square(numpy.square(Num_exp_MuMu))))
      graph2.SetPointEYlow(i,math.sqrt(numpy.square(Err_pre_low_MuMu)/numpy.square(Num_exp_MuMu)+numpy.square(Err_exp_low_MuMu)*numpy.square(Num_pre_MuMu)/numpy.square(numpy.square(Num_exp_MuMu))))
      if i == 4:
        print Eee_up, Emumu_up
def draw_cms_lumi(c1,ytitle):
  t = rt.TLatex()
  t.SetTextAlign(11) # align bottom left corner of text
  t.SetTextColor(rt.kBlack)
  t.SetTextSize(0.04)
  # get top left corner of current pad, and nudge up the y coord a bit
  xcms = rt.c1.GetX1() + rt.c1.GetLeftMargin()
  ycms = rt.c1.GetY2() - rt.c1.GetTopMargin() + 0.01
  xlumi = rt.c1.GetX2() - rt.c1.GetRightMargin()
  xXlabel = rt.c1.GetX2()- rt.c1.GetRightMargin()
  yXlabel =rt.c1.GetY1() + 0.01
  xYlabel =rt.c1.GetX1()+ rt.c1.GetLeftMargin()*0.5
  yYlabel = rt.c1.GetY2()- rt.c1.GetTopMargin()
  cms_label = "preliminary"
  lumi_value = "35.9"
  lumi_unit = "fb"
  energy = 13
  if cms_label is not None:
    t.DrawLatexNDC(xcms,ycms,"#scale[1.25]{#font[61]{CMS}} #scale[1.1]{#font[52]{%s}}" % cms_label)
  if lumi_value:
    t.SetTextSize(0.04)
    t.SetTextAlign(31) # align bottom right
    t.SetTextFont(42) # align bottom right
    t.DrawLatexNDC(xlumi,ycms,"{lumi_str} {lumi_unit}^{{-1}} ({energy} TeV)".format(energy=energy, lumi_str=lumi_value, lumi_unit=lumi_unit))
  t.DrawLatexNDC(xXlabel,yXlabel,"MET cut [GeV]")
  t.SetTextAngle(90)
  t.DrawLatexNDC(xYlabel,yYlabel,ytitle)
  #t.Delete()


def GetHisto (filename, histoname) :
  global xsec, instLumi
  if 'NRBctrl' in histoname:
    h = rt.TH2F()
  else:
    h = rt.TH1F()
  Nevent = rt.TH1F()
  file = rt.TFile.Open("../OUTPUTS/NRB/MERGED/"+filename+'.root')
  file.GetObject("totEventInBaobab_tot",Nevent)
  norm = 1.0
  if not 'Data' in filename:
    norm = instLumi*xsec[filename]/Nevent.Integral();
  else:
    norm = 1
  pointer = file.FindObjectAny(histoname)
  #print 'looking for histo named ',histoname
  #print pointer
  if not pointer == None:
    file.GetObject(histoname,h)
    h_tmp = h.DrawCopy()
    h_tmp.Scale(norm)
    return h_tmp
  else:
    return None

# initiating some lists and dicts
files = ['outputNRB_Data',
'outputNRB_TTJets_DiLept',
'outputNRB_TTJets_SingleLeptFromTbar',
'outputNRB_TTWJetsToLNu',
'outputNRB_TTZToLLNuNu_M-10',
'outputNRB_WWTo2L2Nu',
#'outputNRB_WJetsToLNu',
#'outputNRB_WJetsToLNu_HT-100To200',
#'outputNRB_WJetsToLNu_HT-200To400',
#'outputNRB_WJetsToLNu_HT-400To600',
#'outputNRB_WJetsToLNu_HT-600To800',
#'outputNRB_WJetsToLNu_HT-800To1200',
#'outputNRB_WJetsToLNu_HT-1200To2500',
#'outputNRB_WJetsToLNu_HT-2500ToInf',
'outputNRB_DYJetsToLL_M-50',
'outputNRB_DYJetsToLL_M-10to50',
'outputNRB_ZZToTauTau2Nu',
'outputNRB_ZZToTauTau2Q',
'outputNRB_DYJetsToTauTau_M-50',
'outputNRB_DYJetsToTauTau_M-10to50',
'outputNRB_ST_s-channel_4f_leptonDecays',
'outputNRB_ST_t-channel_antitop_4f_inclusiveDecays',
'outputNRB_ST_t-channel_top_4f_inclusiveDecays',
'outputNRB_ST_tW_antitop_5f_inclusiveDecays',
'outputNRB_ST_tW_top_5f_inclusiveDecays',
'outputNRB_WZTo2L2Q',
'outputNRB_WZTo3LNu',
#'outputNRB_WWW_4F',
'outputNRB_WWZ',
'outputNRB_WZZ',
'outputNRB_ZZZ'
]
#'outputNRB_ZZTo2L2Nu',
#'outputNRB_ZZTo2L2Q',
#'outputNRB_ZZTo4L'] # input Root file names
xsec = {
'outputNRB_TTJets_DiLept':87.31,
'outputNRB_TTJets_SingleLeptFromTbar':182.17,
'outputNRB_TTWJetsToLNu':0.2043,
'outputNRB_TTZToLLNuNu_M-10':0.2529,
'outputNRB_WWTo2L2Nu':12.178,
'outputNRB_WWToLNuQQ':49.997,
'outputNRB_WJetsToLNu':61526.7,
'outputNRB_WJetsToLNu_HT-100To200':1345.,
'outputNRB_WJetsToLNu_HT-200To400':359.7,
'outputNRB_WJetsToLNu_HT-400To600':48.91,
'outputNRB_WJetsToLNu_HT-600To800':12.05,
'outputNRB_WJetsToLNu_HT-800To1200':5.501,
'outputNRB_WJetsToLNu_HT-1200To2500':1.329,
'outputNRB_WJetsToLNu_HT-2500ToInf':0.03216,
'outputNRB_DYJetsToLL_M-50':5765,
'outputNRB_DYJetsToLL_M-10to50':18610.,
'outputNRB_DYJetsToTauTau_M-50':5765,
'outputNRB_DYJetsToTauTau_M-10to50':18610.,
'outputNRB_ST_s-channel_4f_leptonDecays':3.362,
'outputNRB_ST_t-channel_antitop_4f_inclusiveDecays':70.69,
'outputNRB_ST_t-channel_top_4f_inclusiveDecays':70.69,
'outputNRB_ZZToTauTau2Nu':0.564,
'outputNRB_ZZToTauTau2Q':3.22,
'outputNRB_ZZTo4L':1.256,
'outputNRB_ST_tW_antitop_5f_inclusiveDecays':35.6,
'outputNRB_ST_tW_top_5f_inclusiveDecays':35.6,
'outputNRB_WZTo2L2Q':5.595,
'outputNRB_WZTo3LNu':4.42965,
'outputNRB_ZZZ':0.01398 ,
'outputNRB_WZZ':0.05565 ,
'outputNRB_WWZ':0.16510 
#'outputNRB_WWW_4F':0.1651 
}
bins = {'in_bveto':1,'allSB_bveto':2,'upSB_bveto':3,'in_btag':4,'allSB_btag':5,'upSB_btag':6,'in_inclusive':7,'allSB_inclusive':8,'upSB_inclusive':9}
instLumi= 36866.932
histos = ['mt_Inbveto50','mt_Inbveto80','mt_Inbveto125','mt_Outbtag50','mt_Outbtag80','mt_Outbtag125','zmass_bveto125'] #histo names for 1D analysis
channels = ['ee','mumu','emu']
closurelabels= ['Top','WW','allMC','outputNRB_Data']
results = {}
results_2D ={}
for channel in channels:
  for i in range(20):
    for key in bins:
      results_2D['Num_'+key+'_allMC_'+channel+'_METcut'+str(50+i*5)]            = .0
      results_2D['Errsquare_'+key+'_low_allMC_'+channel+'_METcut'+str(50+i*5)]  = .0
      results_2D['Errsquare_'+key+'_up_allMC_'+channel+'_METcut'+str(50+i*5)]   = .0
      results_2D['Num'+key+'_allMC_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Errsquare'+key+'_low_allMC_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key+'_up_allMC_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num'+key+'_allMC_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Errsquare'+key+'_low_allMC_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key+'_up_allMC_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num'+key+'_allMC_'+channel+'_METcut'+str(50+i*5)]       = .0
      results_2D['Errsquare'+key+'_low_allMC_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key+'_up_allMC_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num_'+key          +'_Top_'+channel+'_METcut'+str(50+i*5)]            = .0
      results_2D['Errsquare_'+key+'_low_Top_'+channel+'_METcut'+str(50+i*5)]  = .0
      results_2D['Errsquare_'+key +'_up_Top_'+channel+'_METcut'+str(50+i*5)]   = .0
      results_2D['Num'+key           +'_Top_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Errsquare'+key +'_low_Top_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key  +'_up_Top_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num'+key           +'_Top_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Errsquare'+key +'_low_Top_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key  +'_up_Top_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num'+key           +'_Top_'+channel+'_METcut'+str(50+i*5)]       = .0
      results_2D['Errsquare'+key +'_low_Top_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key  +'_up_Top_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num_'+key          +'_WW_'+channel+'_METcut'+str(50+i*5)]            = .0
      results_2D['Errsquare_'+key+'_low_WW_'+channel+'_METcut'+str(50+i*5)]  = .0
      results_2D['Errsquare_'+key +'_up_WW_'+channel+'_METcut'+str(50+i*5)]   = .0
      results_2D['Num'+key           +'_WW_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Errsquare'+key +'_low_WW_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key  +'_up_WW_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num'+key           +'_WW_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Errsquare'+key +'_low_WW_'+channel+'_METcut'+str(50+i*5)]= .0
      results_2D['Errsquare'+key  +'_up_WW_'+channel+'_METcut'+str(50+i*5)] = .0
      results_2D['Num'+key           +'_WW_'+channel+'_METcut'+str(50+i*5)]       = .0
      results_2D['Errsquare'+key +'_low_WW_'+channel+'_METcut'+str(50+i*5)]= .0
for file in files:
  for histo in histos:
    for channel in channels:
      results['Num_Top_'+histo+'_'+channel] = .0
      results['Errsquare_Top_'+histo+'_'+channel] = .0
      results['Num_WW_'+histo+'_'+channel] = .0
      results['Errsquare_WW_'+histo+'_'+channel] = .0
      results['Num_allMC_'+histo+'_'+channel] = .0
      results['Errsquare_allMC_'+histo+'_'+channel] = .0
#filling the 1D result dict                 
for file in files:
  for histo in histos:
    for channel in channels:
      h = GetHisto(file,histo+'_tot_'+channel)
        # print file+'_'+histo+'_'+channel
      if not h == None :
        Nbins = h.GetNbinsX()
        Err = rt.Double()
        Num = 0
        if 'zmass' in histo:
          Num = h.IntegralAndError(11, 40, Err, "")
        else:
          Num = h.IntegralAndError(0, Nbins+1, Err, "")
        results['Errsquare_'+file+'_'+histo+'_'+channel] = Err*Err
        results['Num_'+file+'_'+histo+'_'+channel] = Num
        if 'TT' in file or 'ST' in file:
          results['Num_Top_'+histo+'_'+channel] += Num
          results['Errsquare_Top_'+histo+'_'+channel] += Err*Err
        if 'WW' in file  or 'WJet' in file:
          results['Num_WW_'+histo+'_'+channel] += Num
          results['Errsquare_WW_'+histo+'_'+channel] += Err*Err
        if not 'Data' in file:
          results['Num_allMC_'+histo+'_'+channel] += Num
          results['Errsquare_allMC_'+histo+'_'+channel] += Err*Err
#filling the 2D result dict
for file in files:
  for channel in channels:
    h2 = GetHisto(file,'mt_shapes_NRBctrl_tot_'+channel)
    if not h2 == None :
      for i in range(20):
        for key in bins:
          if bins[key] < 7:
            # print bins[key]
            results_2D['Num_'+key+'_'+file+'_'+channel+'_METcut'+str(50+i*5)] = h2.GetBinContent(i+2, bins[key])
            results_2D['Errsquare_'+key+'_low_'+file+'_'+channel+'_METcut'+str(50+i*5)] =  h2.GetBinErrorLow(i+1,bins[key]) *h2.GetBinErrorLow(i+1,bins[key])
            results_2D['Errsquare_'+key+'_up_'+file+'_'+channel+'_METcut'+str(50+i*5)] = h2.GetBinErrorUp(i+1,bins[key]) *h2.GetBinErrorUp(i+1,bins[key])
          if bins[key] == 7 :
           #  print bins[key]
            # print key
            results_2D['Num_'+key+'_'+file+'_'+channel+'_METcut'+str(50+i*5)] =     h2.GetBinContent(i+2, 1) +h2.GetBinContent(i+2, 4)
            results_2D['Errsquare_'+key+'_low_'+file+'_'+channel+'_METcut'+str(50+i*5)] = numpy.square(h2.GetBinErrorLow(i+1,1)) +numpy.square(h2.GetBinErrorLow(i+1,4)) 
            results_2D['Errsquare_'+key+'_up_'+file+'_'+channel+'_METcut'+str(50+i*5)] =  numpy.square(h2.GetBinErrorUp(i+1,1)) +numpy.square(h2.GetBinErrorUp(i+1,4)) 
          if bins[key] == 8 :
            # print bins[key]
           #  print key
            results_2D['Num_'+key+'_'+file+'_'+channel+'_METcut'+str(50+i*5)] =     h2.GetBinContent(i+2, 2) +h2.GetBinContent(i+2, 5)
            results_2D['Errsquare_'+key+'_low_'+file+'_'+channel+'_METcut'+str(50+i*5)] = numpy.square(h2.GetBinErrorLow(i+1,2)) +numpy.square(h2.GetBinErrorLow(i+1,5)) 
            results_2D['Errsquare_'+key+'_up_'+file+'_'+channel+'_METcut'+str(50+i*5)] =  numpy.square(h2.GetBinErrorUp(i+1,2)) +numpy.square(h2.GetBinErrorUp(i+1,5)) 
          if bins[key] == 9 :
            results_2D['Num_'+key+'_'+file+'_'+channel+'_METcut'+str(50+i*5)] =          h2.GetBinContent(i+2, 3) +h2.GetBinContent(i+2, 6)
            results_2D['Errsquare_'+key+'_low_'+file+'_'+channel+'_METcut'+str(50+i*5)] = numpy.square(h2.GetBinErrorLow(i+1,3)) +numpy.square(h2.GetBinErrorLow(i+1,6)) 
            results_2D['Errsquare_'+key+'_up_'+file+'_'+channel+'_METcut'+str(50+i*5)] =  numpy.square(h2.GetBinErrorUp(i+1,3)) +numpy.square(h2.GetBinErrorUp(i+1,6)) 
          
          if not 'Data' in file:
            if bins[key] < 7 and not  (key == 'in_bveto' and (channel == 'ee' or channel == 'mumu') and 'WZTo' in file ):
              results_2D['Num_'+key+'_allMC_'+channel+'_METcut'+str(50+i*5)]               += h2.GetBinContent(i+2, bins[key])
              results_2D['Errsquare_'+key+'_low_allMC_'+channel+'_METcut'+str(50+i*5)]     += numpy.square(h2.GetBinErrorLow(i+1,bins[key]))
              results_2D['Errsquare_'+key+'_up_allMC_'+channel+'_METcut'+str(50+i*5)]      += numpy.square(h2.GetBinErrorUp(i+1,bins[key]))
          if 'TT' in file  or 'ST_' in file:
            if bins[key] < 7 :
              results_2D['Num_'+key+'_Top_'+channel+'_METcut'+str(50+i*5)]               += h2.GetBinContent(i+2, bins[key])
              results_2D['Errsquare_'+key+'_low_Top_'+channel+'_METcut'+str(50+i*5)]     += numpy.square(h2.GetBinErrorLow(i+1,bins[key]))
              results_2D['Errsquare_'+key+'_up_Top_'+channel+'_METcut'+str(50+i*5)]      += numpy.square(h2.GetBinErrorUp(i+1,bins[key]))
          if 'WWTo' in file or 'TT' in file  or 'ST_' in file:
            if bins[key] < 7 :
              results_2D['Num_'+key+'_WW_'+channel+'_METcut'+str(50+i*5)]               += h2.GetBinContent(i+2, bins[key])
              if i ==4 and channel == 'emu'  and key == 'upSB_bveto':
                print 'Num_'+key+'_   WW_'+channel+'_METcut'+str(50+i*5),'+=',h2.GetBinContent(i+2, bins[key]),'+-',h2.GetBinErrorLow(i+1,bins[key]),file
              results_2D['Errsquare_'+key+'_low_WW_'+channel+'_METcut'+str(50+i*5)]     += numpy.square(h2.GetBinErrorLow(i+1,bins[key]))
              results_2D['Errsquare_'+key+'_up_WW_'+channel+'_METcut'+str(50+i*5)]      += numpy.square(h2.GetBinErrorUp(i+1,bins[key]))
          if key == 'in_bveto' and (channel == 'ee' or channel == 'mumu'):
            results_2D['Num_'+key+'_allMC_'+channel+'_METcut'+str(50+i*5)] =  results_2D['Num_'+key+'_WW_'+channel+'_METcut'+str(50+i*5)]
            results_2D['Errsquare_'+key+'_low_allMC_'+channel+'_METcut'+str(50+i*5)] = results_2D['Errsquare_'+key+'_low_WW_'+channel+'_METcut'+str(50+i*5)]
            results_2D['Errsquare_'+key+'_up_allMC_'+channel+'_METcut'+str(50+i*5)] = results_2D['Errsquare_'+key+'_up_WW_'+channel+'_METcut'+str(50+i*5)]


print 'Num_upSB_bveto_allMC_ee_METcut70', results_2D['Num_upSB_bveto_allMC_ee_METcut70']
print 'Num_upSB_bveto_   WW_ee_METcut70', results_2D['Num_upSB_bveto_WW_ee_METcut70']
#print results_2D['Num_allSB_allMC_emu_METcut70'],'-',results_2D['Errsquare_allSB_low_allMC_emu_METcut70'],'+',results_2D['Errsquare_allSB_up_allMC_emu_METcut70'],'\n',results_2D['Num_upSB_allMC_emu_METcut70'] ,'-',results_2D['Errsquare_upSB_low_allMC_emu_METcut70'] ,'+',results_2D['Errsquare_upSB_up_allMC_emu_METcut70'],'\n'

ee_gr_data_allSB = rt.TGraphAsymmErrors()
ee_gr_data_upSB = rt.TGraphAsymmErrors()
ee_gr_mc_allSB = rt.TGraphAsymmErrors()
ee_gr_mc_upSB = rt.TGraphAsymmErrors()
mumu_gr_data_allSB = rt.TGraphAsymmErrors()
mumu_gr_data_upSB = rt.TGraphAsymmErrors()
mumu_gr_mc_allSB = rt.TGraphAsymmErrors()
mumu_gr_mc_upSB = rt.TGraphAsymmErrors()
ee_closure_allSB = rt.TGraphAsymmErrors()
mumu_closure_allSB = rt.TGraphAsymmErrors()
ee_closure_upSB = rt.TGraphAsymmErrors()
mumu_closure_upSB = rt.TGraphAsymmErrors()
kmethod_ee = rt.TGraphAsymmErrors()
kmethod_mumu = rt.TGraphAsymmErrors()
kmethod_ee_closure = rt.TGraphAsymmErrors()
kmethod_mumu_closure = rt.TGraphAsymmErrors()

fill_TGraphAsym(ee_gr_data_allSB,mumu_gr_data_allSB,"allSB","btag","data",results_2D)
fill_TGraphAsym(ee_gr_data_upSB,mumu_gr_data_upSB,"upSB","btag","data",results_2D)
fill_TGraphAsym(ee_gr_mc_allSB,mumu_gr_mc_allSB,"allSB","btag","allMC",results_2D)
fill_TGraphAsym(ee_gr_mc_upSB,mumu_gr_mc_upSB,"upSB","btag","allMC",results_2D)
fill_TGraphAsym(ee_closure_allSB,mumu_closure_allSB,"allSB","btag","closure",results_2D)
fill_TGraphAsym(ee_closure_upSB,mumu_closure_upSB,"upSB","btag","closure",results_2D)
fill_TGraphAsym_kMethod(kmethod_ee,kmethod_mumu,"btag",False,results_2D)
fill_TGraphAsym_kMethod(kmethod_ee_closure,kmethod_mumu_closure,"btag",True,results_2D)

c1 = rt.TCanvas("c1","c1",600,600)

pad1 = rt.TPad("pad1","pad1",0.,0.,1.,1.)
pad1.Draw()
pad1.cd()
legend1 = rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym(ee_gr_data_allSB,mumu_gr_data_allSB,"data",legend1)
draw_cms_lumi(pad1,"#alpha value")
c1.Update()
c1.SaveAs("data_allSB.pdf")

pad2 = rt.TPad("pad2","pad2",0.,0.,1.,1.)
pad2.Draw()
pad2.cd()
legend2 = rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym(ee_gr_data_upSB,mumu_gr_data_upSB,"data",legend2)
draw_cms_lumi(pad2,"#alpha value")
c1.Update()
c1.SaveAs("data_upSB.pdf")

pad3 = rt.TPad("pad3","pad3",0.,0.,1.,1.)
pad3.Draw()
pad3.cd()
legend3 = rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym(ee_gr_mc_allSB,mumu_gr_mc_allSB,"mc",legend3)
draw_cms_lumi(pad3,"#alpha value")
c1.Update()
c1.SaveAs("mc_allSB.pdf")


pad4 = rt.TPad("pad4","pad4",0.,0.,1.,1.)
pad4.Draw()
pad4.cd()
legend4= rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym(ee_gr_mc_upSB,mumu_gr_mc_upSB,"mc",legend4)
draw_cms_lumi(pad4,"#alpha value")
c1.Update()
c1.Print("mc_upSB.pdf")

pad5 = rt.TPad("pad5","pad5",0.,0.,1.,1.)
pad5.Draw()
pad5.cd()
legend5= rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym_closure(ee_closure_allSB,mumu_closure_allSB,legend5)
draw_cms_lumi(pad5,"Bias")
c1.Update()
c1.Print("bias_allSB.pdf")

pad6 = rt.TPad("pad6","pad6",0.,0.,1.,1.)
pad6.Draw()
pad6.cd()
legend6= rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym_closure(ee_closure_upSB,mumu_closure_upSB,legend6)
draw_cms_lumi(pad6,"Bias")
c1.Update()
c1.Print("bias_upSB.pdf")

pad7 = rt.TPad("pad7","pad7",0.,0.,1.,1.)
pad7.Draw()
pad7.cd()
legend7= rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym(kmethod_ee,kmethod_mumu,"data kMethod",legend7)
draw_cms_lumi(pad7,"k Value")
c1.Update()
c1.Print("kMethod.pdf")

pad8 = rt.TPad("pad8","pad8",0.,0.,1.,1.)
pad8.Draw()
pad8.cd()
legend8= rt.TLegend(0.8,0.9,0.8,0.9)
draw_TGraphAsym_closure(kmethod_ee_closure,kmethod_mumu_closure,legend8)
draw_cms_lumi(pad8,"Bias")
c1.Update()
c1.Print("kMethod_closure.pdf")



#draw_TGraphAsym(ee_gr_data_upSB,mumu_gr_data_upSB)
#draw_cms_lumi(c1)
#c1.Update()
#c1.SaveAs("data_upSB.pdf")
#
#
#draw_TGraphAsym(ee_gr_mc_allSB,mumu_gr_mc_allSB)
#draw_cms_lumi(c1)
#c1.Print("mc_allSB.pdf")
#
#
#draw_TGraphAsym(ee_gr_mc_upSB,mumu_gr_mc_upSB)
#draw_cms_lumi(c1)
#c1.Print("mc_upSB.pdf")
for sidebandtype in ['upSB','allSB']:
  for Btagveto in ['btag','bveto']:
    fillYields(sidebandtype,Btagveto,results_2D)

print 0.5*math.sqrt(results_2D['Num_in_btag_outputNRB_Data_ee_METcut70']/results_2D['Num_in_btag_outputNRB_Data_mumu_METcut70']),0.5*math.sqrt(results_2D['Num_in_btag_outputNRB_Data_mumu_METcut70']/results_2D['Num_in_btag_outputNRB_Data_ee_METcut70'])


if os.path.exists("AAyields.tex"): 
  os.remove("AAyields.tex")
catalogfile=open("AAyields.tex",'w')
catalogfile.write("\\newpage \n")
catalogfile.write("\\begin{tabular}{l||r|r|r||r||r||r||r} \n")
catalogfile.write("\\hline \\hline\n")
catalogfile.write("\\multicolumn{8}{c}{b-tagged sample} \\\\ \\hline \n")
catalogfile.write("\\hline $ee$  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{ee}^{out}$&  $   \\alpha   $ & $N_{ee}^{SR, est}$& $N_{ee}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_Top_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_Top_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_Top_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_Top_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_Top_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_Top_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_Top_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_Top_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_Top_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_Top_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_Top_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_Top_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_Top_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_Top_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_Top_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_WW_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_WW_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_WW_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_WW_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_WW_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_WW_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_WW_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_WW_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_WW_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_WW_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_WW_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_WW_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_WW_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_WW_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_WW_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_allMC_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_allMC_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_allMC_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_allMC_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_allMC_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_allMC_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_allMC_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_allMC_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_allMC_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_allMC_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_allMC_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_allMC_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_allMC_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_allMC_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_allMC_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline $\\mu\\mu $  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{\\mu\\mu}^{out}$&  $   \\alpha   $ & $N_{\\mu\\mu}^{SR, est}$& $N_{\\mu\\mu}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_Top_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_Top_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_Top_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_Top_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_Top_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_Top_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_Top_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_Top_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_Top_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_Top_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_Top_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_Top_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_Top_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_Top_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_Top_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_WW_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_WW_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_WW_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_WW_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_WW_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_WW_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_WW_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_WW_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_WW_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_WW_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_WW_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_WW_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_WW_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_WW_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_WW_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_allMC_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_allMC_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_allMC_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_allMC_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_allMC_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_allMC_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_allMC_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_allMC_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_allMC_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_allMC_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_allMC_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline \\hline\n")
catalogfile.write("\\multicolumn{8}{c}{b-vetoed sample} \\\\ \\hline \n")
catalogfile.write("\\hline $ee$  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{ee}^{out}$&  $   \\alpha   $ & $N_{ee}^{SR, est}$& $N_{ee}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_bveto_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_bveto_Top_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_Top_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_Top_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_bveto_Top_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_bveto_Top_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_bveto_Top_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_bveto_Top_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_bveto_Top_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_bveto_Top_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_Top_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_Top_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_Top_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_bveto_Top_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_bveto_Top_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_bveto_Top_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_bveto_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_bveto_WW_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_WW_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_WW_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_bveto_WW_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_bveto_WW_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_bveto_WW_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_bveto_WW_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_bveto_WW_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_bveto_WW_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_WW_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_WW_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_WW_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_bveto_WW_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_bveto_WW_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_bveto_WW_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_bveto_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_bveto_allMC_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_allMC_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_allMC_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_bveto_allMC_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_bveto_allMC_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_bveto_allMC_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_bveto_allMC_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_bveto_allMC_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_bveto_allMC_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_allMC_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_allMC_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_allMC_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_bveto_allMC_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_bveto_allMC_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_bveto_allMC_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline $\\mu\\mu $  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{\\mu\\mu}^{out}$&  $   \\alpha   $ & $N_{\\mu\\mu}^{SR, est}$& $N_{\\mu\\mu}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_bveto_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_bveto_Top_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_Top_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_Top_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_bveto_Top_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_bveto_Top_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_bveto_Top_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_bveto_Top_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_bveto_Top_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_bveto_Top_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_Top_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_Top_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_Top_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_bveto_Top_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_bveto_Top_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_bveto_Top_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_bveto_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_bveto_WW_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_WW_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_WW_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_bveto_WW_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_bveto_WW_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_bveto_WW_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_bveto_WW_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_bveto_WW_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_bveto_WW_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_WW_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_WW_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_WW_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_bveto_WW_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_bveto_WW_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_bveto_WW_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_bveto_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_up_allMC_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_bveto_low_allMC_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_bveto_allMC_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_bveto_allMC_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_allMC_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_allMC_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_bveto_allMC_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline \\hline\n")
catalogfile.write("\\end{tabular} \n")
catalogfile.write("\\newpage \n")
catalogfile.write("\n")

catalogfile.write("\\begin{tabular}{l||r|r|r||r||r||r||r} \n")
catalogfile.write("\\hline \\hline\n")
catalogfile.write("\\multicolumn{8}{c}{b-tagged sample} \\\\ \\hline \n")
catalogfile.write("\\hline $ee$  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{ee}^{out}$&  $   \\alpha   $ & $N_{ee}^{SR, est}$& $N_{ee}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_btag_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_btag_Top_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_Top_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_Top_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_btag_Top_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_Top_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_Top_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_btag_Top_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_btag_Top_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_btag_Top_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_in_bveto_Top_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_Top_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_Top_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_btag_Top_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_btag_Top_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_btag_Top_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_btag_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_btag_WW_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_WW_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_WW_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_btag_WW_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_WW_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_WW_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_btag_WW_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_btag_WW_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_btag_WW_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_in_bveto_WW_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_WW_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_WW_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_btag_WW_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_btag_WW_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_btag_WW_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_btag_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_btag_allMC_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_allMC_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_allMC_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_btag_allMC_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_allMC_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_allMC_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_btag_allMC_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_btag_allMC_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_btag_allMC_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_in_bveto_allMC_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_allMC_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_allMC_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_btag_allMC_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_btag_allMC_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_btag_allMC_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline $\\mu\\mu $  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{\\mu\\mu}^{out}$&  $   \\alpha   $ & $N_{\\mu\\mu}^{SR, est}$& $N_{\\mu\\mu}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_btag_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_btag_Top_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_Top_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_Top_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_btag_Top_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_Top_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_Top_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_btag_Top_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_btag_Top_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_btag_Top_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_Top_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_Top_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_Top_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_btag_Top_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_btag_Top_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_btag_Top_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_btag_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_btag_WW_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_WW_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_WW_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_btag_WW_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_WW_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_WW_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_btag_WW_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_btag_WW_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_btag_WW_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_WW_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_WW_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_WW_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_btag_WW_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_btag_WW_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_btag_WW_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_btag_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_btag_allMC_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_up_allMC_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_btag_low_allMC_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_btag_allMC_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_allMC_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_btag_allMC_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_btag_allMC_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_allMC_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_allMC_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_btag_allMC_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_btag_allMC_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline \n")
catalogfile.write("\\multicolumn{8}{c}{b-vetoed sample} \\\\ \\hline \n")
catalogfile.write("\\hline $ee$  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{ee}^{out}$&  $   \\alpha   $ & $N_{ee}^{SR, est}$& $N_{ee}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_bveto_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_bveto_Top_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_Top_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_Top_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_bveto_Top_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_bveto_Top_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_bveto_Top_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_bveto_Top_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_bveto_Top_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_bveto_Top_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_Top_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_Top_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_Top_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_bveto_Top_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_bveto_Top_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_bveto_Top_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_bveto_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_bveto_WW_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_WW_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_WW_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_bveto_WW_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_bveto_WW_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_bveto_WW_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_bveto_WW_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_bveto_WW_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_bveto_WW_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_WW_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_WW_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_WW_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_bveto_WW_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_bveto_WW_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_bveto_WW_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_bveto_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_bveto_allMC_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_allMC_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_allMC_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_bveto_allMC_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_bveto_allMC_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_bveto_allMC_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_bveto_allMC_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_bveto_allMC_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_bveto_allMC_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_allMC_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_allMC_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_allMC_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_bveto_allMC_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_bveto_allMC_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_bveto_allMC_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline $\\mu\\mu $  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{\\mu\\mu}^{out}$&  $   \\alpha   $ & $N_{\\mu\\mu}^{SR, est}$& $N_{\\mu\\mu}^{SR, exp}$ & $f_{est/exp}$ \\\\  \\hline \n")
catalogfile.write("Top & $"+"%.1f" %results_2D['Num_in_bveto_Top_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_Top_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_Top_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_bveto_Top_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_Top_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_Top_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_bveto_Top_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_Top_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_Top_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_bveto_Top_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_bveto_Top_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_bveto_Top_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_bveto_Top_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_bveto_Top_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_bveto_Top_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_Top_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_Top_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_Top_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_bveto_Top_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_bveto_Top_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_bveto_Top_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("Top+WW & $"+"%.1f" %results_2D['Num_in_bveto_WW_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_WW_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_WW_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_bveto_WW_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_WW_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_WW_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_bveto_WW_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_WW_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_WW_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_bveto_WW_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_bveto_WW_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_bveto_WW_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_bveto_WW_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_bveto_WW_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_bveto_WW_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_WW_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_WW_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_WW_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_bveto_WW_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_bveto_WW_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_bveto_WW_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("allMC & $"+"%.1f" %results_2D['Num_in_bveto_allMC_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_allMC_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_allMC_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_upSB_bveto_allMC_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_allMC_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_allMC_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_upSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_up_allMC_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_upSB_bveto_low_allMC_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_upSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_bveto_allMC_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_upSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_upSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_upSB_bveto_allMC_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_upSB_bveto_allMC_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_upSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_upSB_bveto_allMC_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_upSB_bveto_allMC_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_upSB_bveto_allMC_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_upSB_bveto_allMC_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline \n")
catalogfile.write("\\end{tabular} \n")
catalogfile.write("\\newpage \n")
catalogfile.write("\n")

catalogfile.write("\\begin{tabular}{|cc|c|c|}  \n")
catalogfile.write("\\hline \n")
catalogfile.write("\\multicolumn{2}{|c|}{Channel} & b-tag, up side band & b-tag, all side band \\\\ \\hline  \n")
catalogfile.write("\\multirow{2}{0.5in}{$ee$} &Data & "+"$"+"%.3f" %results_2D['alpha_upSB_btag_outputNRB_Data_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_outputNRB_Data_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_outputNRB_Data_ee_METcut70']+"}$"+" & "+"$"+"%.3f" %results_2D['alpha_allSB_btag_outputNRB_Data_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_outputNRB_Data_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_outputNRB_Data_ee_METcut70']+"}$"+" \\\\ & MC & "+"$"+"%.3f" %results_2D['alpha_upSB_btag_allMC_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_allMC_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_allMC_ee_METcut70']+"}$"+" & "+"$"+"%.3f" %results_2D['alpha_allSB_btag_allMC_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_allMC_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_allMC_ee_METcut70']+"}$"+"  \\\\ \\hline \n")
catalogfile.write("\\multirow{2}{0.5in}{$\\mu\\mu$} &Data & "+"$"+"%.3f" %results_2D['alpha_upSB_btag_outputNRB_Data_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_outputNRB_Data_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_outputNRB_Data_mumu_METcut70']+"}$"+" & "+"$"+"%.3f" %results_2D['alpha_allSB_btag_outputNRB_Data_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_outputNRB_Data_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_outputNRB_Data_mumu_METcut70']+"}$"+" \\\\ & MC & "+"$"+"%.3f" %results_2D['alpha_upSB_btag_allMC_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_upSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_upSB_btag_allMC_mumu_METcut70']+"}$"+" & "+"$"+"%.3f" %results_2D['alpha_allSB_btag_allMC_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_allMC_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_allMC_mumu_METcut70']+"}$"+"  \\\\ \\hline \n")
catalogfile.write("\\hline  \n")
catalogfile.write("\\end{tabular}  \n")
catalogfile.write("\\newpage \n")
catalogfile.write("\n")
catalogfile.write("\\begin{tabular}{l||r|r|r|r||r||r||r} \n")
catalogfile.write("\\hline \\hline $ee$  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{ee}^{out}$&  $   \\alpha   $ & $N_{ee}^{SR, est}$& $N_{ee}^{SR, exp}$ & $f_{est/exp}$ \\\\   \n")
catalogfile.write("\\hline NRB & $"+"%.1f" %results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_outputNRB_Data_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_outputNRB_Data_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_outputNRB_Data_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_outputNRB_Data_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_outputNRB_Data_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_outputNRB_Data_ee_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_outputNRB_Data_ee_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_outputNRB_Data_ee_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_outputNRB_Data_ee_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_outputNRB_Data_ee_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_outputNRB_Data_ee_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_outputNRB_Data_ee_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_outputNRB_Data_ee_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_outputNRB_Data_ee_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_outputNRB_Data_ee_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_outputNRB_Data_ee_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_outputNRB_Data_ee_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_outputNRB_Data_ee_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_outputNRB_Data_ee_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_outputNRB_Data_ee_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\hline $\\mu\\mu $  & $N_{e\\mu}^{in}$ & $N_{e\\mu}^{out}$  & $N_{\\mu\\mu}^{out}$&  $   \\alpha   $ & $N_{\\mu\\mu}^{SR, est}$& $N_{\\mu\\mu}^{SR, exp}$ & $f_{est/exp}$ \\\\   \n")
catalogfile.write("\\hline NRB & $"+"%.1f" %results_2D['Num_in_bveto_outputNRB_Data_emu_METcut125']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_up_outputNRB_Data_emu_METcut125'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_in_bveto_low_outputNRB_Data_emu_METcut125'])+"} $ & $ "+"%.1f" %results_2D['Num_allSB_btag_outputNRB_Data_emu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_outputNRB_Data_emu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_outputNRB_Data_emu_METcut70'])+"} $ & $  "+"%.1f" %results_2D['Num_allSB_btag_outputNRB_Data_mumu_METcut70']+"^{+"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_up_outputNRB_Data_mumu_METcut70'])+"}_{-"+"%.1f" %math.sqrt(results_2D['Errsquare_allSB_btag_low_outputNRB_Data_mumu_METcut70'])+"}$ & $"+"%.3f" %results_2D['alpha_allSB_btag_outputNRB_Data_mumu_METcut70']+"^{+"+"%.3f" %results_2D['Erralpha_high_allSB_btag_outputNRB_Data_mumu_METcut70']+"}_{-"+"%.3f" %results_2D['Erralpha_low_allSB_btag_outputNRB_Data_mumu_METcut70']+"}$ & $"+"%.1f" %results_2D['Num_pre_allSB_btag_outputNRB_Data_mumu_METcut70']+"^{+"+"%.1f" %results_2D['ErrNum_pre_high_allSB_btag_outputNRB_Data_mumu_METcut70']+"}_{-"+"%.1f" %results_2D['ErrNum_pre_low_allSB_btag_outputNRB_Data_mumu_METcut70']+"}$ & $ "+"%.1f" %results_2D['Num_exp_allSB_bveto_outputNRB_Data_mumu_METcut125']+"^{+"+"%.1f" %results_2D['ErrNum_exp_high_allSB_bveto_outputNRB_Data_mumu_METcut125']+"}_{-"+"%.1f" %results_2D['ErrNum_exp_low_allSB_bveto_outputNRB_Data_mumu_METcut125']+"}$ & $"+"%.2f" %results_2D['ratio_allSB_btag_outputNRB_Data_mumu_METcut70']+"^{+"+"%.2f" %results_2D['ratio_err_high_allSB_btag_outputNRB_Data_mumu_METcut70']+"}_{-"+"%.2f" %results_2D['ratio_err_low_allSB_btag_outputNRB_Data_mumu_METcut70']+"}$ \\\\ \\hline \n" )
catalogfile.write("\\end{tabular} \n")

catalogfile.close()


