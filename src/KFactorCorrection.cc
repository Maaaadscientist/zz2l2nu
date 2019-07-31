#include <KFactorCorrection.h>

#include <stdexcept>
#include <sstream>
#include <string>

#include <FileInPath.h>
#include <Logger.h>

#include <TFile.h>
#include <TLorentzVector.h>


KFactorCorrection::KFactorCorrection(Dataset &dataset, Options const &options)
    : genPartPt_{dataset.Reader(), "GenPart_pt"},
      genPartEta_{dataset.Reader(), "GenPart_eta"},
      genPartPhi_{dataset.Reader(), "GenPart_phi"},
      genPartMass_{dataset.Reader(), "GenPart_mass"},
      genPartMomId_{dataset.Reader(), "GenPart_genPartIdxMother"} {

  auto const settingsNode = dataset.Info().Parameters()["k_factor"];

  if (settingsNode and not settingsNode.IsNull()) {
    enabled_ = true;
    auto const typeLabel = settingsNode.as<std::string>();

    if (typeLabel != "ggF") {
      std::ostringstream message;
      message << "Unknown type \"" << typeLabel << "\" for k factor.";
      throw std::runtime_error(message.str());
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

    for (int i = 0; i < genPartPt_.GetSize(); i++) {
      if (genPartMomId_[i] != 23)
        continue;

      TLorentzVector lepton;
      lepton.SetPtEtaPhiM(genPartPt_[i], genPartEta_[i], genPartPhi_[i],
        genPartMass_[i]);
      higgs += lepton;
      numberOfLepton++;
    }

    if (numberOfLepton != 4)
      LOG_WARN << "WARNING: Found " << numberOfLepton
          << " generator-level leptons while 4 is expected. Higgs boson is not "
          "reconstructed correctly.";

    return higgs.M();
  }
  else
    return 0.;
}


double KFactorCorrection::operator()() const {
  if (enabled_)
    return kfactorGraph_->Eval(HiggsMass());
  else
    return 1.;
}

