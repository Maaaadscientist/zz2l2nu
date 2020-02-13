#include <istream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <map>

#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TH2F.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TLegend.h>

using namespace std;

void DrawPreliminary(double iLumi, double iEcm, double L, double B, double R, double T, bool isSim=false, bool preliminary=true){

                char LumiText[1024];  sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"fb", iEcm);
                TPaveText* T1 = new TPaveText(1.0-R-0.50, 1.0-T-0.02, 1.02-R, 1.0-0.002, "NDC");
                T1->SetTextFont(43); T1->SetTextSize(23);   T1->SetTextAlign(31);
                T1->SetFillColor(0); T1->SetFillStyle(0);   T1->SetBorderSize(0);
                T1->AddText(LumiText);  T1->Draw();

                TPaveText* T2 = new TPaveText(L+0.005, 1.0-T-0.08, L+0.20, 1.0-T-0.013, "NDC");
                T2->SetTextFont(63); T2->SetTextSize(30);   T2->SetTextAlign(11);
                T2->SetFillColor(0); T2->SetFillStyle(0);   T2->SetBorderSize(0);
                T2->AddText("CMS"); T2->Draw();

                if(preliminary){ //Bellow CMS
                TPaveText* T3 = new TPaveText(L+0.005, 1.0-T-0.183, L+0.20, 1.0-T-0.013, "NDC");
                T3->SetTextFont(53); T3->SetTextSize(23);   T3->SetTextAlign(11);
                T3->SetFillColor(0); T3->SetFillStyle(0);   T3->SetBorderSize(0);
                T3->AddText("Preliminary"); T3->Draw();
                }

                if(isSim){
                TPaveText* T5 = new TPaveText(L+0.085, 1.0-T-0.08, L+0.32, 1.0-T-0.013, "NDC");
                T5->SetTextFont(53); T5->SetTextSize(23);   T5->SetTextAlign(11);
                T5->SetFillColor(0); T5->SetFillStyle(0);   T5->SetBorderSize(0);
                T5->AddText("Simulation"); T5->Draw();
                }
}

void DrawPreliminary(double iLumi, double iEcm, TAttPad* pad=NULL, bool isSim=false,   bool preliminary=true){
              if     (pad )DrawPreliminary(iLumi, iEcm, pad->GetLeftMargin(), pad->GetBottomMargin(), pad->GetRightMargin(), pad->GetTopMargin() , isSim, preliminary);
              else if(gPad)DrawPreliminary(iLumi, iEcm, gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin(), isSim, preliminary);
              else         DrawPreliminary(iLumi, iEcm, 0.15, 0.15, 0.15, 0.15, isSim, preliminary);
}

void NonResonantBckg_TotalUncert_Plotter(){

	//string      MassWindow = "InRegion";
	string      MassWindow = "AllSB";
        string            bTag = "bTag";
        string        NameFile = "normalized.root"; 
        //string         MainDir = "/afs/cern.ch/work/l/lyuan/cmsarea/hzz2l2v/CMSSW_8_0_14/src/UserCode/llvv_fwk/test/hzz2l2v/";
        //string         MainDir = "~georgia/work/public/analysis/preapp/";

        //string process[10] = {"data", "Top", "WW", "W#rightarrow l#nu", "Z#rightarrow #tau#tau_filt15", "ZZ#rightarrow Z#tau#tau_filt15", "WZ", "ZVV", "Z#rightarrow ee-#mu#mu_filt1113", "ZZ_filt1113"}; 
        string process[] = {
        "Data",
        "TTJets_DiLept",
        "TTWJetsToLNu",
        "TTZToLLNuNu_M-10",
        "ST_s-channel_4f_leptonDecays",
        "ST_t-channel_antitop_4f_inclusiveDecays",
        "ST_t-channel_top_4f_inclusiveDecays",
        "ST_tW_antitop_5f_inclusiveDecays",
        "ST_tW_top_5f_inclusiveDecays",
        "WWTo2L2Nu",
        "WJetsToLNu_HT-100To200",
        "WJetsToLNu_HT-200To400",
        "WJetsToLNu_HT-400To600",
        "WJetsToLNu_HT-600To800",
        "WJetsToLNu_HT-800To1200",
        "WJetsToLNu_HT-1200To2500",
        "WJetsToLNu_HT-2500ToInf",
        "WZTo2L2Q",
        "WZTo3LNu",
        "WWZ",
        "WZZ",
        "ZZZ",
        "DYJetsToTauTau_M-50",
        "ZZToTauTau2Nu",
        "ZZToTauTau2Q",
        "DYJetsToLL_M-50",
        "ZZTo2L2Nu",
        "ZZTo2L2Q",
        "ZZTo4L"
        };

        TH1F *emu_SR[31];
        TH1F *ee_SR[31];
        TH1F *mumu_SR[31];
	
	TH2F *emu_mt_2D_NRBctrl[31];
	TH2F *ee_mt_2D_NRBctrl[31];
	TH2F *mumu_mt_2D_NRBctrl[31];


	TH2F *emu_mt_2D_MC;
	TH2F *ee_mt_2D_MC;
	TH2F *mumu_mt_2D_MC;

	TH1F *h_emu_SR_1D;
	TH1F *h_ee_SR_1D;
	TH1F *h_mumu_SR_1D;

	//string InFile = MainDir + NameFile;
	string InFile = NameFile;
	TFile *mf = TFile::Open(InFile.c_str());


        for(int i=0; i<22; i++){

		if(i==1){
			emu_mt_2D_MC = (TH2F* )mf->Get(("mt_shapes_NRBctrl_tot_emu_"+process[i]).c_str());
			ee_mt_2D_MC = (TH2F* )mf->Get(("mt_shapes_NRBctrl_tot_ee_"+process[i]).c_str());
			mumu_mt_2D_MC = (TH2F* )mf->Get(("mt_shapes_NRBctrl_tot_mumu_"+process[i]).c_str());

			h_emu_SR_1D = (TH1F *)mf->Get(("mt_Inbveto125_tot_emu_"+process[i]).c_str());
			h_ee_SR_1D = (TH1F *)mf->Get(("mt_Inbveto125_tot_ee_"+process[i]).c_str());
			h_mumu_SR_1D = (TH1F *)mf->Get(("mt_Inbveto125_tot_mumu_"+process[i]).c_str());
		}
		emu_mt_2D_NRBctrl[i] = (TH2F* )mf->Get(("mt_shapes_NRBctrl_tot_emu_"+process[i]).c_str());
		ee_mt_2D_NRBctrl[i] = (TH2F* )mf->Get(("mt_shapes_NRBctrl_tot_ee_"+process[i]).c_str());
		mumu_mt_2D_NRBctrl[i] = (TH2F* )mf->Get(("mt_shapes_NRBctrl_tot_mumu_"+process[i]).c_str());

		emu_SR[i] = (TH1F *)mf->Get(("mt_Inbveto125_tot_emu_"+process[i]).c_str());
		ee_SR[i] = (TH1F *)mf->Get(("mt_Inbveto125_tot_ee_"+process[i]).c_str());
		mumu_SR[i] = (TH1F *)mf->Get(("mt_Inbveto125_tot_mumu_"+process[i]).c_str());
        }

        for(int i=2; i<22 ; i++){ //attention
		emu_mt_2D_MC->Add(emu_mt_2D_NRBctrl[i]);
		ee_mt_2D_MC->Add(ee_mt_2D_NRBctrl[i]);
		mumu_mt_2D_MC->Add(mumu_mt_2D_NRBctrl[i]);

		if(emu_SR[i])h_emu_SR_1D->Add(emu_SR[i]);
		if(i<17){
			if(ee_SR[i])h_ee_SR_1D->Add(ee_SR[i]);
			if(mumu_SR[i])h_mumu_SR_1D->Add(mumu_SR[i]);
		}
	}
	

	TH2F * h_emu_mt_2D;
	TH2F * h_ee_mt_2D;
	TH2F * h_mumu_mt_2D;

	double N_emu_SR, E_emu_SR;
	double N_ee_SR, E_ee_SR;
	double N_mumu_SR, E_mumu_SR;

	int nbins = h_emu_SR_1D->GetNbinsX();
	N_emu_SR = h_emu_SR_1D->IntegralAndError(0, nbins+1, E_emu_SR, "");
	
	N_ee_SR = h_ee_SR_1D->IntegralAndError(0, nbins+1, E_ee_SR, "");
	N_mumu_SR = h_mumu_SR_1D->IntegralAndError(0, nbins+1, E_mumu_SR, "");
	
	cout<<" ee SR: "<<N_ee_SR<<" err: "<<E_ee_SR<<" mumu SR: "<<N_mumu_SR<<" err: "<<E_mumu_SR<<endl;

	int bin = 5;
	if(MassWindow.find("AllSB")!=string::npos) bin = 5;
	else if(MassWindow.find("upSB")!=string::npos) bin = 6;


	//Loop on PorcessType: MC or Data
	for( int i = 0; i < 1; ++i ){

                string FinalName  = "Plot_13TeV_" + MassWindow + "_TotalUncert_" + bTag + ".pdf";
                //string FinalName  = "Plot_13TeV_" + MassWindow + "_DataMC_RelativeDiff_" + bTag + ".pdf";
                //string FinalName  = "Plot_13TeV_" + MassWindow + "_all_" + bTag + ".pdf";
		string     TagDir = "Method_NonClosure_" + MassWindow + "_" + bTag + "/";
		string  MainDirMC = MainDir + TagDir;
                string OutTxtFile = "TotalUncertanties_" + MassWindow  + "_" + bTag + ".txt";
                ofstream OutFile;
                OutFile.open (OutTxtFile.c_str());

                string     cname = "NRB_NonClosure_" + MassWindow  + "_13TeV";
                TLegend     *leg = new TLegend( 0.56681, 0.667678, 0.85273, 0.863429);
             
                TMultiGraph  *mg = new TMultiGraph();
            
                TCanvas       *c = new TCanvas( cname.c_str(), " ", 800, 800);
                TPad         *p1 = new TPad( "Histo", "Histo", 0.0, 0.0, 1.0, 1.0);
          
                string gr_title = " "; 
                      
                string ee_channel, mumu_channel; 
                double       met_cut[18],         met_err[18];
                double      ee_value[18],      mumu_value[18];
	        double  ee_value_err[18],  mumu_value_err[18];        
                double      ee_alpha_MC[18],      mumu_alpha_MC[18];
                double      ee_alpha_data[18],      mumu_alpha_data[18];
	        double  ee_alpha_data_err[18],  mumu_alpha_data_err[18];        
	        double  ee_alpha_data_maxDiff[18],  mumu_alpha_data_maxDiff[18];        
	        double  ee_alpha_data_maxDiff_rela[18],  mumu_alpha_data_maxDiff_rela[18];        
                double      ee_totalUncert[18],      mumu_totalUncert[18];
	        double  ee_alpha_dataMC_Diff_rela[18],  mumu_alpha_dataMC_Diff_rela[18];        


		h_emu_mt_2D = (TH2F* )emu_mt_2D_MC->Clone("h_emu_mt_2D");
		h_ee_mt_2D = (TH2F* )ee_mt_2D_MC->Clone("h_ee_mt_2D");
		h_mumu_mt_2D = (TH2F* )mumu_mt_2D_MC->Clone("h_mumu_mt_2D");
 
                //Loop on index cut 1->19
		for( int t = 0; t < 19; ++t ){

                    ///////////////////
                    // Running on MC //
                    //////////////////

                    double ee_val, ee_err;
                    double mumu_val, mumu_err;
          
                    double ee_alpha_val, ee_alpha_err;
                    double mumu_alpha_val, mumu_alpha_err;
                    int Line = 0;

                    //Loop on File
                    ee_alpha_val = h_ee_mt_2D->GetBinContent(t+2, bin)/h_emu_mt_2D->GetBinContent(t+2, bin);
		    ee_alpha_err = ee_alpha_val*sqrt(pow(h_ee_mt_2D->GetBinError(t+2, bin)/h_ee_mt_2D->GetBinContent(t+2, bin), 2) + pow(h_emu_mt_2D->GetBinError(t+2, bin)/h_emu_mt_2D->GetBinContent(t+2, bin), 2));


                    mumu_alpha_val = h_mumu_mt_2D->GetBinContent(t+2, bin)/h_emu_mt_2D->GetBinContent(t+2, bin);
		    mumu_alpha_err = mumu_alpha_val*sqrt(pow(h_mumu_mt_2D->GetBinError(t+2, bin)/h_mumu_mt_2D->GetBinContent(t+2, bin), 2) + pow(h_emu_mt_2D->GetBinError(t+2, bin)/h_emu_mt_2D->GetBinContent(t+2, bin), 2));

		    cout<<" ee alpha: "<<ee_alpha_val<<" err: "<<ee_alpha_err<<" mumu alpha: "<<mumu_alpha_val<<" err: "<<mumu_alpha_err<<endl;
		 ee_alpha_MC[t  ] = ee_alpha_val;
		 mumu_alpha_MC[t  ] = mumu_alpha_val;

		 ee_val = N_emu_SR * ee_alpha_val;
	  	 ee_err = ee_val*sqrt( pow(E_emu_SR/N_emu_SR, 2) + pow(ee_alpha_err/ee_alpha_val, 2));

		 mumu_val = N_emu_SR * mumu_alpha_val;
	  	 mumu_err = mumu_val*sqrt( pow(E_emu_SR/N_emu_SR, 2) + pow(mumu_alpha_err/mumu_alpha_val, 2));

                 met_cut[t  ] = (double)( t*5 + 45 );
                 //std::cout << "mumu_alpha: " << mumu_val << "; ee_alpha: " << ee_val << "; met_cut: " << (double)( t*5 + 45 ) << "; t: " << t << std::endl;
                      
                 ee_value[t  ]        = (ee_val - N_ee_SR)/N_ee_SR; 
                 ee_value_err[t  ]    = (ee_val/N_ee_SR)  * sqrt( pow(ee_err/ee_val, 2) + pow(E_ee_SR/N_ee_SR, 2) );
                 mumu_value[t  ]        = (mumu_val - N_mumu_SR)/N_mumu_SR; 
                 mumu_value_err[t  ]    = (mumu_val/N_mumu_SR) * sqrt( pow(mumu_err/mumu_val, 2) + pow(E_mumu_SR/N_mumu_SR, 2) );
                 met_err[t  ]         = 0.;

		 //OutFile << " MET cut: "<<met_cut[t  ]<<" ** ee ** NonClosure: "<<ee_value[t  ]<<"  Err: "<<ee_value_err[t  ]<< endl;
		 //OutFile << " MET cut: "<<met_cut[t  ]<<" -- mumu -- NonClosure: "<<mumu_value[t  ]<<"  Err: "<<mumu_value_err[t  ]<< endl;
		 //OutFile << endl;


                    ///////////////////
                    // Running on Data //
                    //////////////////
                    //
                    ee_alpha_val = ee_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin)/emu_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin);
		    ee_alpha_err = ee_alpha_val*sqrt(pow(ee_mt_2D_NRBctrl[0]->GetBinError(t+2, bin)/ee_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin), 2) + pow(emu_mt_2D_NRBctrl[0]->GetBinError(t+2, bin)/emu_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin), 2));


                    mumu_alpha_val = mumu_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin)/emu_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin);
		    mumu_alpha_err = mumu_alpha_val*sqrt(pow(mumu_mt_2D_NRBctrl[0]->GetBinError(t+2, bin)/mumu_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin), 2) + pow(emu_mt_2D_NRBctrl[0]->GetBinError(t+2, bin)/emu_mt_2D_NRBctrl[0]->GetBinContent(t+2, bin), 2));

		    ee_alpha_data[t  ] = ee_alpha_val;
		    ee_alpha_data_err[t  ] = ee_alpha_err;

		    mumu_alpha_data[t  ] = mumu_alpha_val;
		    mumu_alpha_data_err[t  ] = mumu_alpha_err;
		    cout<<" Data: ee alpha: "<<ee_alpha_val<<" err: "<<ee_alpha_err<<" mumu alpha: "<<mumu_alpha_val<<" err: "<<mumu_alpha_err<<endl;

		    ee_alpha_data_maxDiff[t  ] = 0;
		    mumu_alpha_data_maxDiff[t  ] = 0;

            }//End Loop Met Cut
  
           for( unsigned int n = 0; n < 18; n++){
             for( unsigned int m = 0; m < 18 ; m++){
		  double ee_tmpDiff = fabs(ee_alpha_data[n] - ee_alpha_data[m]);
		  if(ee_tmpDiff > ee_alpha_data_maxDiff[n])ee_alpha_data_maxDiff[n] = ee_tmpDiff;

		  double mumu_tmpDiff = fabs(mumu_alpha_data[n] - mumu_alpha_data[m]);
		  if(mumu_tmpDiff > mumu_alpha_data_maxDiff[n])mumu_alpha_data_maxDiff[n] = mumu_tmpDiff;
      cout<<ee_alpha_data[n]<<"-"<<ee_alpha_data[m] <<" and " <<mumu_alpha_data[n]<<"-"<<mumu_alpha_data[m]<<endl;
       cout<<" data ee alpha difference: "<<ee_alpha_data_maxDiff[n]<<" mumu difference: "<<mumu_alpha_data_maxDiff[n]<<endl;
	     }

           }

           for( unsigned int n = 0; n < 18; n++){
		ee_alpha_data_maxDiff_rela[n] = ee_alpha_data_maxDiff[n]/ee_alpha_data[n];
		mumu_alpha_data_maxDiff_rela[n] = mumu_alpha_data_maxDiff[n]/mumu_alpha_data[n];

		ee_alpha_dataMC_Diff_rela[n] = fabs(ee_alpha_data[n] - ee_alpha_MC[n])/ee_alpha_data[n];
		mumu_alpha_dataMC_Diff_rela[n] = fabs(mumu_alpha_data[n] - mumu_alpha_MC[n])/mumu_alpha_data[n];

		ee_totalUncert[n] = sqrt( pow(ee_value[n], 2) + pow(ee_alpha_data_err[n]/ee_alpha_data[n], 2) + pow(ee_alpha_data_maxDiff[n]/ee_alpha_data[n], 2));
		mumu_totalUncert[n] = sqrt( pow(mumu_value[n], 2) + pow(mumu_alpha_data_err[n]/mumu_alpha_data[n], 2) + pow(mumu_alpha_data_maxDiff[n]/mumu_alpha_data[n], 2));
//		ee_totalUncert[n] = sqrt( pow(ee_value[n], 2) + pow(ee_alpha_data_err[n]/ee_alpha_data[n], 2) + pow(ee_alpha_data_maxDiff[n]/ee_alpha_data[n], 2) + pow(ee_alpha_dataMC_Diff_rela[n], 2));
//		mumu_totalUncert[n] = sqrt( pow(mumu_value[n], 2) + pow(mumu_alpha_data_err[n]/mumu_alpha_data[n], 2) + pow(mumu_alpha_data_maxDiff[n]/mumu_alpha_data[n], 2) + pow(mumu_alpha_dataMC_Diff_rela[n], 2));

		OutFile << " MET cut: "<<met_cut[n]<<" **  ee  ** "<<" NonClosure: "<<ee_value[n]<<endl;
		OutFile << " MET cut: "<<met_cut[n]<<" **  ee  ** "<<" StatUncert: "<<fabs(ee_alpha_data_err[n]/ee_alpha_data[n])<<endl;
		OutFile << " MET cut: "<<met_cut[n]<<" **  ee  ** "<<" MaxDiff: "<<fabs(ee_alpha_data_maxDiff[n]/ee_alpha_data[n])<<endl;
		OutFile << " MET cut: "<<met_cut[n]<<" **  ee  ** "<<" ***** Total Uncertainty: "<<ee_totalUncert[n]<<endl;
		OutFile << endl;

		OutFile << " MET cut: "<<met_cut[n]<<" -- mumu -- "<<" NonClosure: "<<mumu_value[n]<<endl;
		OutFile << " MET cut: "<<met_cut[n]<<" -- mumu -- "<<" StatUncert: "<<fabs(mumu_alpha_data_err[n]/mumu_alpha_data[n])<<endl;
		OutFile << " MET cut: "<<met_cut[n]<<" -- mumu -- "<<" MaxDiff: "<<fabs(mumu_alpha_data_maxDiff[n]/mumu_alpha_data[n])<<endl;
		OutFile << " MET cut: "<<met_cut[n]<<" -- mumu -- "<<" NonClosure: "<<mumu_value[n]<<"Total Uncertainty: "<<mumu_totalUncert[n]<< endl;
		OutFile << endl;
		OutFile << endl;

	   }

           //Choose the lowest value of alpha both for data and MC
           double min_err_ee_mc      = ee_value_err[0];
           double min_err_mumu_mc    = mumu_value_err[0];
           int index_ee_mc     = 0;
           int index_mumu_mc   = 0;
        
           for( unsigned int n = 0; n < 18; n++){
             if( min_err_ee_mc > ee_value_err[n] && ee_value_err[n] > 0 ){
                min_err_ee_mc = ee_value_err[n];
                index_ee_mc   = n;
             }
           }

           for( unsigned int p = 0; p < 18; p++){
             if( min_err_mumu_mc > mumu_value_err[p] && mumu_value_err[p] > 0 ){
                min_err_mumu_mc = mumu_value_err[p];
                index_mumu_mc   = p;
             }
           }
  
           OutFile << "Min Err, at cut " << index_ee_mc  << " => ee Err: "  << min_err_ee_mc     << "; Alfa Value: " << ee_value[index_ee_mc] << "\n";
           OutFile << "Min Err, at cut " << index_mumu_mc << " => mumu Err: "  << min_err_mumu_mc   << "; Alfa Value: " << mumu_value[index_mumu_mc] << "\n"; 

//           TGraphErrors  *ee_gr         = new TGraphErrors(  18,  met_cut,         ee_value,    met_err,         ee_value_err);
//           TGraphErrors  *mumu_gr       = new TGraphErrors(  18,  met_cut,       mumu_value,    met_err,       mumu_value_err);

           TGraph  *ee_gr         = new TGraph(  18,  met_cut,         ee_totalUncert);
           TGraph  *mumu_gr       = new TGraph(  18,  met_cut,       mumu_totalUncert);
//           TGraph  *ee_gr         = new TGraph(  18,  met_cut,         ee_alpha_dataMC_Diff_rela);
//           TGraph  *mumu_gr       = new TGraph(  18,  met_cut,         mumu_alpha_dataMC_Diff_rela);

   
           string hname = " ";
           string ee_leg_name = "ee channel"; 
           string mumu_leg_name = "mumu channel";

           ee_gr->SetTitle(gr_title.c_str()); 
           ee_gr->SetMarkerColor(602);
           ee_gr->SetMarkerStyle(22);
           ee_gr->SetMarkerSize(1.5);                 

           mumu_gr->SetTitle(gr_title.c_str());
           mumu_gr->SetMarkerColor(634);
           mumu_gr->SetMarkerStyle(22);
           mumu_gr->SetMarkerSize(1.5);

           mg->Add( ee_gr );
           mg->Add( mumu_gr );

           leg->AddEntry(   ee_gr,   ee_leg_name.c_str(),  "P");
           leg->AddEntry( mumu_gr, mumu_leg_name.c_str(),  "P");

           p1->SetFillColor(0);
           p1->SetBorderMode(0);
           p1->SetBorderSize(2);
           p1->SetTickx(1);
           p1->SetTicky(1);
           p1->SetLeftMargin(0.10);
           p1->SetRightMargin(0.05);
           p1->SetTopMargin(0.05);
           p1->SetBottomMargin(0.10);
           p1->SetFrameFillStyle(0);
           p1->SetFrameBorderMode(0);
           p1->SetFrameFillStyle(0);
           p1->SetFrameBorderMode(0);
           p1->Draw();
           p1->cd();

           mg->Draw("APE");
           mg->SetTitle(gr_title.c_str());
           mg->GetXaxis()->SetTitle("MET cut [GeV]");
           mg->GetYaxis()->SetTitle("Total Uncertainty");
           //mg->GetYaxis()->SetTitle("Relative Data-MC #alpha");
           mg->GetYaxis()->SetRangeUser(-0.0, 0.2);
           mg->GetYaxis()->SetTitleOffset(1.4);
                  
           leg->Draw();
           leg->SetFillColor(0);
           leg->SetBorderSize(0);

           if( i == 0 ) { DrawPreliminary(35.9, 13, p1, true, false); }
           else if( i == 1 ){ DrawPreliminary(35.9, 13, p1); }

           c->SaveAs(FinalName.c_str());

           OutFile.close();	

	}//End loop on ProcessType
}                                                                                                     
