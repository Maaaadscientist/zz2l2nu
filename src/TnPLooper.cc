#define TnPLooper_cxx

#include <LooperMain.h>

#include <bitset> 

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TTreeReaderArray.h>


namespace tnp {

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

};  // namespace tnp

namespace TnPobject
{
  bool TriggerMatchResult (ULong64_t MuHltMatch){
    bool MatchResult;
    bool IsoMu24 = MuHltMatch &(1<<2);
    bool IsotkMu24 = MuHltMatch &(1<<3);
    MatchResult = IsoMu24 || IsotkMu24;
    return MatchResult;
  }
  bool selectElectrons(std::vector<tnp::Electron> & tagElectrons, std::vector<tnp::Electron> & probeElectrons,TTreeReaderArray<float> const &ElCh, TTreeReaderArray<float> const &ElPt, TTreeReaderArray<float> const &ElEta, TTreeReaderArray<float> const &ElPhi, TTreeReaderArray<float> const &ElE, TTreeReaderArray<unsigned int> const &ElId, TTreeReaderArray<float> const &ElEtaSc, TTreeReaderArray<float> const &ElPfIsoRho)
  {
    for(int i = 0 ; i<ElPt.GetSize() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false;
      tnp::Electron currentLepton;
      currentLepton.lvector.SetPtEtaPhiE(ElPt[i],ElEta[i],ElPhi[i],ElE[i]);
      currentLepton.Pt = ElPt[i];
      currentLepton.Phi = ElPhi[i];
      currentLepton.Eta = ElEta[i];
      currentLepton.PfIsoRho = ElPfIsoRho[i];
      currentLepton.EtaSc = ElEtaSc[i];
      currentLepton.Ch = ElCh[i];
      currentLepton.Seq = i;

      passId = (ElId[i] & (1<<17))?1:0;//need to check this later
      currentLepton.Id = passId;

      int eta = fabs(ElEtaSc[i]);//I took the supercluster eta since it's really the geometry which is taken here.
      passEta = (eta<=2.5 && (eta>=1.5660 || eta<=1.4442));
      currentLepton.E = fabs(ElEta[i]);
      if(eta>=1.5660 && ElPfIsoRho[i]<0.0646) passIso = true;
      if(eta<=1.4442 && ElPfIsoRho[i]<0.0354) passIso = true; //Numbers are taken from llvv_fwk and have not been checked.
      passPt = (currentLepton.lvector.Pt() >=10);
      if(passIso && passId ) tagElectrons.push_back(currentLepton);
      if(passEta && passPt) probeElectrons.push_back(currentLepton);
    }
    return true;
  }

  bool selectMuons(std::vector<tnp::Muon> & tagMuons, std::vector<tnp::Muon> & probeMuons,TTreeReaderArray<float> const &MuCh, TTreeReaderArray<float> const &MuPt, TTreeReaderArray<float> const &MuEta, TTreeReaderArray<float> const &MuPhi, TTreeReaderArray<float> const &MuE, TTreeReaderArray<unsigned int> const &MuId, TTreeReaderArray<unsigned int> const &MuIdTight, TTreeReaderArray<float> const &MuPfIso,TTreeReaderArray<unsigned int> const &MuHltMatch)
  {
    for(int i = 0 ; i<MuPt.GetSize() ; i++){
      bool passEta = false, passIso = false, passId = false, passPt = false ;
      tnp::Muon currentLepton;
      currentLepton.lvector.SetPtEtaPhiE(MuPt[i],MuEta[i],MuPhi[i],MuE[i]);

      currentLepton.Pt = MuPt[i];
      currentLepton.Phi = MuPhi[i];
      currentLepton.Eta = MuEta[i];
      currentLepton.E = MuE[i];
      currentLepton.PfIso = MuPfIso[i];
      currentLepton.Id = MuId[i];
      currentLepton.IdTight = MuIdTight[i]& (1<<0);
      //std::cout<< "The hltmatch number is "<<std::bitset<sizeof(int)*8>(MuHltMatch[i])<<std::endl;
      currentLepton.Ch = MuCh[i];
      currentLepton.Seq = i;
      currentLepton.HltMatch =1;// TriggerMatchResult(MuHltMatch[i]);
      // for(int k=0;k<sizeof(currentLepton.TriggerInf) / sizeof(currentLepton.TriggerInf[0]);k++){currentLepton.TriggerInf[k]=((MuHltMatch[i])>>k)&1;}
      passId = MuIdTight[i] & (1<<0); //Look at the first vertex, hence the bit 0.
      float eta = fabs(MuEta[i]);
      currentLepton.E = eta;
      passEta = (eta<=2.4);
      passIso = (MuPfIso[i]<0.15); //Numbers are taken from llvv_fwk and have not been checked.
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
  std::vector<tnp::Electron>  tagElectrons ;
  std::vector<tnp::Electron>  probeElectrons;
  std::vector<tnp::Muon>  tagMuons ;
  std::vector<tnp::Muon>  probeMuons ;



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

  cout << "nb of entries in the input file =" << nentries << endl;

  //###############################################################
  //##################     EVENT LOOP STARTS     ##################
  //###############################################################

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    fReader.SetEntry(jentry);

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
