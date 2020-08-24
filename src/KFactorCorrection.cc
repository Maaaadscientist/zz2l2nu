#include <KFactorCorrection.h>

#include <string>

#include <FileInPath.h>
#include <HZZException.h>
#include <Logger.h>

#include <TFile.h>
#include <TLorentzVector.h>


KFactorCorrection::KFactorCorrection(Dataset &dataset, Options const &)
    : genPartPt_{dataset.Reader(), "GenPart_pt"},
      genPartEta_{dataset.Reader(), "GenPart_eta"},
      genPartPhi_{dataset.Reader(), "GenPart_phi"},
      genPartMass_{dataset.Reader(), "GenPart_mass"},
      genPartStatus_{dataset.Reader(), "GenPart_status"},
      genPartStatusFlags_{dataset.Reader(), "GenPart_statusFlags"} {

  auto const settingsNode = dataset.Info().Parameters()["k_factor"];

  if (settingsNode and not settingsNode.IsNull()) {
    enabled_ = true;
    auto const typeLabel = settingsNode.as<std::string>();

    if (typeLabel != "ggF") {
      HZZException exception;
      exception << "Unknown type \"" << typeLabel << "\" for k factor.";
      throw exception;
    }
  } else
    enabled_ = false;

  if (enabled_) {
    LOG_DEBUG << "Will apply k factors.";
    TFile kFactorFile(FileInPath::Resolve(
      "corrections/Kfactor_Collected_"
      "ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root").c_str());
    kfactorGraph_.reset(dynamic_cast<TGraph *>(
      kFactorFile.Get("kfactor_Nominal")));
  }
  else
    LOG_DEBUG << "Will not apply k factors.";
}


double KFactorCorrection::HiggsMass() const {
  if (enabled_) {
    TLorentzVector higgs;
    int numberOfLepton = 0;

    for (int i = 0; i < int(genPartPt_.GetSize()); i++) {
      // Status: 1=stable
      // flags bits are: 0 : isPrompt, 8 : fromHardProcess
      if (genPartStatus_[i] != 1 || (genPartStatusFlags_[i] & 1 << 0) == 0 ||
          (genPartStatusFlags_[i] & 1 << 8) == 0 ) 
        continue;

      TLorentzVector lepton;
      lepton.SetPtEtaPhiM(genPartPt_[i], genPartEta_[i], genPartPhi_[i],
        genPartMass_[i]);
      higgs += lepton;
      numberOfLepton++;
    }

    if (numberOfLepton != 4) {
      HZZException exception;
      exception << "Found " << numberOfLepton
          << " generator-level leptons while 4 is expected. Higgs boson is not "
          "reconstructed correctly.";
      throw exception;
    }
    return higgs.M();
  }
  else
    return 0.;
}


double KFactorCorrection::NominalWeight() const {
  if (enabled_)
    return kfactorGraph_->Eval(HiggsMass());
  else
    return 1.;
}
