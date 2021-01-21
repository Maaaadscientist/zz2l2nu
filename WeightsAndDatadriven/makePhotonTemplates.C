void makePhotonTemplates(bool isPhoton = true) {
  TFile *f1 = isPhoton? TFile::Open("../templates_photon_2018_v2.root") : TFile::Open("../templates_SR_2018_Asimov_v3.root"); // Update paths if needed
  TFile *fTransfer = TFile::Open("../data/InstrMetReweighting/meanWeights_2018.root");
  std::map<std::string,TH1*> transferFunction;
  transferFunction["eq0jets"] = (TH1*) fTransfer->Get("mean_weights_tot_eq0jets");
  transferFunction["eq1jets"] = (TH1*) fTransfer->Get("mean_weights_tot_eq1jets");
  transferFunction["geq2jets"] = (TH1*) fTransfer->Get("mean_weights_tot_geq2jets");
  TH1* histModel = (TH1*) f1->Get("geq2jets/data_obs/nominal");
  std::cout << "transferFunction has " << transferFunction["geq2jets"]->GetNbinsX() << " bins. histModel has " << histModel->GetNbinsX() << " bins." << std::endl;
  int nbins = 14;
  for (int i = 1 ; i <= nbins ; i++) {
    TString fileName, fileNumber;
    fileNumber.Form("%d",i);
    if (isPhoton) fileName = "template_photon_2018_v3_bin_"+fileNumber+".root";
    else fileName = "temporary_SR_bin_"+fileNumber+".root";
    TFile *fout = new TFile(fileName,"recreate");
    fout->cd();
    const auto *keys = f1->GetListOfKeys();
    for (const auto &&object : *keys) {
      std::string name = dynamic_cast<TKey*>(object)->GetName();
      const auto *key = dynamic_cast<TKey*>(object);
      TDirectory *dir = nullptr;
      f1->GetObject(key->GetName(), dir);
      TDirectory* subD = fout->mkdir(key->GetName());
      if(object->IsFolder()) {
        if (dir != nullptr) {
          const auto *keys2 = dir->GetListOfKeys();
          for (const auto &&object2 : *keys2) {
            if (!isPhoton and i != 1) continue;
            std::string name2 = dynamic_cast<TKey*>(object2)->GetName();
            const auto *key2 = dynamic_cast<TKey*>(object2);
            TDirectory *dir2 = nullptr;
            f1->GetObject((TString)key->GetName()+"/"+(TString)key2->GetName(), dir2);
            TDirectory* subD2 = subD->mkdir(key2->GetName());
            if(object2->IsFolder()) {
              if (dir2 != nullptr) {
                const auto *keys3 = dir2->GetListOfKeys();
                subD2->cd();
                for (const auto &&object3 : *keys3) {
                  std::string name3 = dynamic_cast<TKey*>(object3)->GetName();
                  const auto *key3 = dynamic_cast<TKey*>(object3);
                  if (!((std::string)key3->GetClassName() == "TH1D")) continue;
                  TH1D *h = nullptr;
                  f1->GetObject((TString)key->GetName()+"/"+(TString)key2->GetName()+"/"+(TString)key3->GetName(), h);
                  if (isPhoton) {
                    for (int iBin = 0 ; iBin <= h->GetNbinsX()+1 ; iBin++) {
                      if (iBin != i) {
                        h->SetBinContent(iBin,0);
                        h->SetBinError(iBin,0);
                      }
                    }
                  }
                  h->Write();
                  h->Delete();
                }
              }
            }
          }
          if (name == "emu") continue;
          TH1D* histInstrMET = (TH1D*)histModel->Clone("nominal");
          histInstrMET->Reset();
          if (isPhoton) {
            histInstrMET->SetBinContent(i,1);
          }
          else {
            if (transferFunction[name]->GetBinContent(i) == 0) {
              histInstrMET->SetBinContent(i,1);
            }
            else {
              histInstrMET->SetBinContent(i,transferFunction[name]->GetBinContent(i));
            }
          }
          histInstrMET->SetBinError(i,0);
          histInstrMET->SetName("nominal");
          TDirectory* subInstrMET = subD->mkdir("InstrMET_bin"+fileNumber);
          subInstrMET->cd();
          histInstrMET->Write();
          histInstrMET->Delete();
          histModel->SetDirectory(0);
        }
      }
    }
  }
}
