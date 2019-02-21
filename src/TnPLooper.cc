#define TnPLooper_cxx

#include <LooperMain.h>

#include <bitset> 
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>

struct Electron
{
  float Pt,Eta,Phi,E;
  float EtaSc, PfIsoRho;
  float Ch; 
  int Id;
  int Seq;// avoiding a pair containing 2 same lep
  TLorentzVector lvector;     
};

struct Muon
{
  float Pt,Eta,Phi,E;
  float PfIso;
  float Ch;
  int Id, IdTight,HltMatch;
  // int TriggerInf[20];
  int Seq;// avoiding a pair containing 2 same lep
  TLorentzVector lvector;
};
namespace TnPobject
{
  bool TriggerMatchResult (ULong64_t MuHltMatch){
    bool MatchResult;
    bool IsoMu24 = MuHltMatch &(1<<2);
    bool IsotkMu24 = MuHltMatch &(1<<3);
    MatchResult = IsoMu24 || IsotkMu24;
    return MatchResult;
  }
  bool selectElectrons(std::vector<Electron> & tagElectrons, std::vector<Electron> & probeElectrons,std::vector<float> *ElCh, std::vector<float> *ElPt, std::vector<float> *ElEta, std::vector<float> *ElPhi, std::vector<float> *ElE, std::vector<unsigned int> *ElId, std::vector<float> *ElEtaSc, std::vector<float> *ElPfIsoRho)
  {
    for(int i = 0 ; i<ElPt->size() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false;
      Electron currentLepton; 
      currentLepton.lvector.SetPtEtaPhiE(ElPt->at(i),ElEta->at(i),ElPhi->at(i),ElE->at(i));
      currentLepton.Pt = ElPt->at(i);
      currentLepton.Phi = ElPhi->at(i);
      currentLepton.Eta = ElEta->at(i);
      currentLepton.PfIsoRho = ElPfIsoRho->at(i);
      currentLepton.EtaSc = ElEtaSc->at(i);
      currentLepton.Ch = ElCh->at(i);
      currentLepton.Seq = i;

      passId = (ElId->at(i) & (1<<17))?1:0;//need to check this later
      currentLepton.Id = passId;

      int eta = fabs(ElEtaSc->at(i));//I took the supercluster eta since it's really the geometry which is taken here.
      passEta = (eta<=2.5 && (eta>=1.5660 || eta<=1.4442));
      currentLepton.E = fabs(ElEta->at(i));
      if(eta>=1.5660 && ElPfIsoRho->at(i)<0.0646) passIso = true;
      if(eta<=1.4442 && ElPfIsoRho->at(i)<0.0354) passIso = true; //Numbers are taken from llvv_fwk and have not been checked.
      passPt = (currentLepton.lvector.Pt() >=10);
      if(passIso && passId ) tagElectrons.push_back(currentLepton);
      if(passEta && passPt) probeElectrons.push_back(currentLepton);
    }
    return true;
  }

  bool selectMuons(std::vector<Muon> & tagMuons, std::vector<Muon> & probeMuons,std::vector<float> *MuCh, std::vector<float> *MuPt, std::vector<float> *MuEta, std::vector<float> *MuPhi, std::vector<float> *MuE, std::vector<unsigned int> *MuId, std::vector<unsigned int> *MuIdTight, std::vector<float> *MuPfIso,std::vector<unsigned int> *MuHltMatch)
  {
    for(int i = 0 ; i<MuPt->size() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false ;
      Muon currentLepton; 
      currentLepton.lvector.SetPtEtaPhiE(MuPt->at(i),MuEta->at(i),MuPhi->at(i),MuE->at(i));

      currentLepton.Pt = MuPt->at(i);
      currentLepton.Phi = MuPhi->at(i);
      currentLepton.Eta = MuEta->at(i);
      currentLepton.E = MuE->at(i);
      currentLepton.PfIso = MuPfIso->at(i);
      currentLepton.Id = MuId->at(i);
      currentLepton.IdTight = MuIdTight->at(i)& (1<<0);
      //std::cout<< "The hltmatch number is "<<std::bitset<sizeof(int)*8>(MuHltMatch->at(i))<<std::endl;
      currentLepton.Ch = MuCh->at(i);
      currentLepton.Seq = i;
      currentLepton.HltMatch =1;// TriggerMatchResult(MuHltMatch->at(i));
      // for(int k=0;k<sizeof(currentLepton.TriggerInf) / sizeof(currentLepton.TriggerInf[0]);k++){currentLepton.TriggerInf[k]=((MuHltMatch->at(i))>>k)&1;}
      passId = MuIdTight->at(i) & (1<<0); //Look at the first vertex, hence the bit 0.
      float eta = fabs(MuEta->at(i));
      currentLepton.E = eta;
      passEta = (eta<=2.4);
      passIso = (MuPfIso->at(i)<0.15); //Numbers are taken from llvv_fwk and have not been checked.
      passPt = (currentLepton.lvector.Pt() >=10);
      if(passIso && passId && currentLepton.HltMatch) tagMuons.push_back(currentLepton);
      if(passEta && passPt) probeMuons.push_back(currentLepton);
    }
    return true;
  }

}

void LooperMain::Loop_TnP()
{
  if (fChain == 0) return;

  //###############################################################
  //################## DECLARATION OF HISTOGRAMS ##################
  //###############################################################
  TFile *outFile = new TFile(outputFile_,"RECREATE");
  TDirectoryFile *dir = new TDirectoryFile("tpTree","tpTree","",outFile);
  gDirectory->cd("tpTree");
  TTree *tpTree = new TTree("fitter_tree","fitter_tree");
  std::vector<Electron>  tagElectrons ;
  std::vector<Electron>  probeElectrons;
  std::vector<Muon>  tagMuons ;
  std::vector<Muon>  probeMuons ;



  // TString variables[]={"mass","tag_pt","pt","tag_eta","eta","tag_phi","phi","tag_abseta","abseta","tag_charge","charge","tkIso","Tight2012","tag_HltMatch","HltMatch","IsoMu22","IsoTkMu22","IsoMu24","IsoTkMu24","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_leg17","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_leg8","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_filter","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_leg17","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_leg8","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_filter"};
  TString variables[]={"mass","tag_pt","pt","tag_eta","eta","tag_phi","phi","tag_abseta","abseta","tag_charge","charge","tkIso","Tight2012","tag_HltMatch","HltMatch"} ;
  Int_t varnum = sizeof(variables) / sizeof(variables[0]);
  // enum varnames
  // {mass,tag_pt,pt,tag_eta,eta,tag_phi,phi,tag_abseta,abseta,tag_charge,charge,tkIso,Tight2012,tag_HltMatch,HltMatch,IsoMu22,IsoTkMu22,IsoMu24,IsoTkMu24,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_leg17,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_leg8,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_filter,HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_leg17,HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_leg8,HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_filter };

  enum varnames
  {mass,tag_pt,pt,tag_eta,eta,tag_phi,phi,tag_abseta,abseta,tag_charge,charge,tkIso,Tight2012,tag_HltMatch,HltMatch};
  int name = Tight2012;
  Float_t fill_float[Tight2012];

  Int_t fill_int[3];

  for(int i=0;i<varnum;i++){

    if (i<name) tpTree->Branch(variables[i],&fill_float[i],variables[i]+"/F");

    else tpTree->Branch(variables[i],&fill_int[i-name],variables[i]+"/I");

  }

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  cout << "nb of entries in the input file =" << nentries << endl;

  //###############################################################
  //##################     EVENT LOOP STARTS     ##################
  //###############################################################

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry % 10000 ==0) cout << jentry << " of " << nentries << endl;
    if (0){
      TnPobject::selectElectrons(tagElectrons, probeElectrons, ElCh,ElPt, ElEta, ElPhi, ElE, ElId, ElEtaSc,ElPfIsoRho);
      for(int i=0;i<tagElectrons.size();i++)
        for(int j=0;j<probeElectrons.size();j++)
        {
          if(tagElectrons[i].Ch== probeElectrons[j].Ch) continue;  
          if(tagElectrons[i].Seq == probeElectrons[j].Seq) continue;// avoiding a pair containing 2 same lepton
          float mll;
          mll= (tagElectrons[i].lvector+probeElectrons[j].lvector).M();
          if(mll<60||mll>120) continue;
          fill_float[mass        ] = mll;
          fill_float[tag_pt      ] = tagElectrons  [i].Pt   ;
          fill_float[pt          ] = probeElectrons[j].Pt   ;
          fill_float[tag_eta     ] = tagElectrons  [i].Eta  ;
          fill_float[eta         ] = probeElectrons[j].Eta  ;
          fill_float[tag_phi     ] = tagElectrons  [i].Phi  ;
          fill_float[phi         ] = probeElectrons[j].Phi  ;
          fill_float[tag_abseta  ] = tagElectrons  [i].E    ;
          fill_float[abseta      ] = probeElectrons[j].E    ;
          fill_float[tag_charge  ] = tagElectrons  [i].Ch   ;
          fill_float[charge      ] = probeElectrons[j].Ch   ;
          fill_float[tkIso       ] = probeElectrons[j].PfIsoRho;
          fill_int  [0] = probeElectrons[j].Id   ;
          tpTree->Fill();
        }
    }
    if(1){
      TnPobject::selectMuons(tagMuons, probeMuons, MuCh,MuPt, MuEta, MuPhi, MuE, MuId, MuIdTight, MuPfIso,MuHltMatch);
      for(int i=0;i<tagMuons.size();i++)
        for(int j=0;j<probeMuons.size();j++)
        {
          if(tagMuons[i].Seq == probeMuons[j].Seq) continue;// avoiding a pair containing 2 same lep
          if(tagMuons[i].Ch == probeMuons[j].Ch) continue;
          float mll;
          mll= (tagMuons[i].lvector+probeMuons[j].lvector).M();
          if(mll<60||mll>120) continue;

          fill_float[mass        ] = mll;
          fill_float[tag_pt      ] = tagMuons  [i].Pt   ;
          fill_float[pt          ] = probeMuons[j].Pt   ;
          fill_float[tag_eta     ] = tagMuons  [i].Eta  ;
          fill_float[eta         ] = probeMuons[j].Eta  ;
          fill_float[tag_phi     ] = tagMuons  [i].Phi  ;
          fill_float[phi         ] = probeMuons[j].Phi  ;
          fill_float[tag_abseta  ] = tagMuons  [i].E    ;
          fill_float[abseta      ] = probeMuons[j].E    ;
          fill_float[tag_charge  ] = tagMuons  [i].Ch   ;
          fill_float[charge      ] = probeMuons[j].Ch   ;
          fill_float[tkIso       ] = probeMuons[j].PfIso;
          fill_int  [0] = probeMuons[j].IdTight           ;
          //    fill_int  [1] = tagMuons[i].HltMatch           ;
          //   fill_int  [2] = probeMuons[j].HltMatch         ;
          // for(int k=0;k<sizeof(probeMuons[j].TriggerInf)/sizeof(probeMuons[j].TriggerInf[0]);k++){fill_int[k+3]=probeMuons[j].TriggerInf[k];}



          tpTree->Fill();
        }
    }


    tagMuons.clear();
    probeMuons.clear();
    tagElectrons.clear();
    probeElectrons.clear();


  }

  //###############################################################
  //##################        END OF LOOP        ##################
  //###############################################################

  tpTree->Write();
  outFile->Close();

}
