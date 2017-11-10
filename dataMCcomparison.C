float instLumi;

struct MCentry{
  TString nameSample;
  TString legendEntry;
  TString fileSuffix;
  float crossSection;
  int color;
  TFile *sampleFile;
  MCentry(TString theNameSample, TString theLegendEntry, TString theFileSuffix, float theCrossSection, int theColor)
    : nameSample(theNameSample), legendEntry(theLegendEntry), fileSuffix(theFileSuffix), crossSection(theCrossSection), color(theColor)
  {}
};

void drawTheHisto(TFile *dataFile, std::vector<MCentry> allMCsamples, TString theHistoName){
  TH1F *MZ_data = (TH1F*) dataFile->Get(theHistoName);
  TH1F *totEventInBaobab_tot_data = (TH1F*) dataFile->Get("totEventInBaobab_tot");
  cout << "the tot events data =" << totEventInBaobab_tot_data->Integral() << endl;


  TCanvas *c0 = new TCanvas("c0","coucou",600,800);
  TPad *pad =new TPad("haut","haut",0,0.25,1,1);
  pad->SetNumber(1);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetTicky();
  cout << pad->GetBottomMargin() << endl;
  pad->SetBottomMargin(0.006);
  pad->Draw();

  TPad *pad2 =new TPad("bas","bas",0,0,1,0.25);
  pad2->SetNumber(2);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();

  c0->cd(1);
  TLegend *t = new TLegend(0.79,0.66,0.89,0.89);
  t->SetLineColor(0);

  MZ_data->SetMarkerColor(kBlack);
  MZ_data->SetLineColor(kBlack);
  MZ_data->Draw("E1");

  TH1F* MChistos[10];
  TH1F* sumMC;
  int iteHisto=0;
  TString lastLegend = "";
  THStack *stackMCsamples = new THStack("stackMCsamples","Stacked MC");
  for (MCentry theMCentry: allMCsamples){
  //  cout << "doing " << theMCentry.nameSample << endl;
    MChistos[iteHisto] = (TH1F*) (theMCentry.sampleFile)->Get(theHistoName);
    if (MChistos[iteHisto] == 0) continue;
  //  cout << "found" << endl;
    TH1F *totEventInBaobab = (TH1F*) (theMCentry.sampleFile)->Get("totEventInBaobab_tot");
    float norm = instLumi*theMCentry.crossSection/totEventInBaobab->Integral()*0.9;//FIXME here the 0.9 corresponds to the double muon trigger efficiency
  //  cout << "scale is " << norm << endl;
    MChistos[iteHisto]->Scale(norm);
    if (iteHisto==0) sumMC = (TH1F*) MChistos[iteHisto]->Clone("sumHisto");
    else sumMC->Add(MChistos[iteHisto]);
    MChistos[iteHisto]->SetLineColor(theMCentry.color);
    MChistos[iteHisto]->SetFillColor(theMCentry.color);
    stackMCsamples->Add(MChistos[iteHisto]);
    if (lastLegend !=theMCentry.legendEntry){
      t->AddEntry(MChistos[iteHisto], theMCentry.legendEntry, "F");
      lastLegend = theMCentry.legendEntry;
    }
    //MChistos[iteHisto]->Draw("HIST:same");
    delete totEventInBaobab;
    iteHisto++;
  }


  MZ_data->Draw("E1:same");
  stackMCsamples->Draw("HIST:same");
  MZ_data->Draw("E1:same");
  t->Draw();

  c0->cd(2);

  TH1F *ratio = (TH1F*) MZ_data->Clone("ratio");
  ratio->Sumw2();
  ratio->Divide(MZ_data, sumMC, 1,1);
  ratio->SetMaximum(1.3);
  ratio->SetMinimum(0.7);
  ratio->SetTitle("");
  ratio->GetYaxis()->SetTitle("");
  ratio->GetYaxis()->SetLabelSize(0.08);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetLabelOffset(0.02);
  ratio->Draw("E1");

  c0->Print("plots/"+theHistoName+".png");
  c0->cd(1);
  c0->cd(1)->SetLogy();
  MZ_data->SetMinimum(0.01);
  MZ_data->Draw("E1:same");

  c0->Print("plots/"+theHistoName+"_log.png");
//  for (int i=0 ; i < iteHisto ; i++){  delete MChistos[i];}
//  delete stackMCsamples;
}


void dataMCcomparison(){
  TString curentDirectory="merged_2ndNov";
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);

  std::vector<MCentry> allMCsamples;
  allMCsamples.push_back(MCentry("DY",        "DY",    "DYJetsToLL",    5765,    833));
  allMCsamples.push_back(MCentry("TTbar",     "Top",   "TTJets",       87.31,      8));
  allMCsamples.push_back(MCentry("WZTo2L2Q",  "WZ",    "WZTo2L2Q",     5.595,    594));
  allMCsamples.push_back(MCentry("WZTo3LNu",  "WZ",    "WZTo3LNu",   4.42965,    594));
  allMCsamples.push_back(MCentry("WWTo2L2Nu", "WW",    "WWTo2L2Nu",   12.178,    590));
  allMCsamples.push_back(MCentry("ZZTo2L2Q",  "ZZ",    "ZZTo2L2Q",      3.22,    595));
  allMCsamples.push_back(MCentry("ZZTo2L2Nu", "ZZ",    "ZZTo2L2Nu",    0.564,    595));
  allMCsamples.push_back(MCentry("ZZTo4L",    "ZZ",    "ZZTo4L",       1.256,    595));

  instLumi=16916.98*9.03906/9.26099;//FIXME 16916.16 = int. lumi of the double muon sample;
  //FIXME 9.03906~nb of data event in the baobab actually used; 9.26099=nb of event in the baobab from the crab report
  //FIXME the 2 last numbers should not be here because ALL the DATA should be included ! 

  TFile *dataFile = new TFile(curentDirectory+"/output_DoubleMuon.root");

  /*for (int i=0 ; i<allMCsamples.size(); i++){
    allMCsamples.at(i).sampleFile  = new TFile("merged_newWeights/output_"+allMCsamples.at(i).fileSuffix+".root");
  }*/
  for (MCentry &theEntry: allMCsamples){
    theEntry.sampleFile = new TFile(curentDirectory+"/output_"+theEntry.fileSuffix+".root");
  }


TIter listPlots(dataFile->GetListOfKeys());
 TKey *keyPlot;
 while ((keyPlot = (TKey*)listPlots())) {
     TString typeObject = keyPlot->GetClassName();
     TString nomObject = keyPlot->GetTitle();
     if (nomObject.Contains("totEventInBaobab")) continue;
     cout << "nom=" << typeObject << " title=" << nomObject << endl;
     drawTheHisto(dataFile, allMCsamples, nomObject);
  }

//drawTheHisto(dataFile, allMCsamples, "M_Z_tot_mumu");
//  drawTheHisto(dataFile, allMCsamples, "eventflow_tot");
//  drawTheHisto(dataFile, allMCsamples, "MET_beforeMETcut");
//drawTheHisto(dataFile, DYfile, "M_Z_tot_mumu_eq0jets");
//drawTheHisto(dataFile, DYfile, "jetCategory_tot");




}
