#include <KFactorCorrection.h>

#include <stdexcept>
#include <sstream>
#include <string>

#include <FileInPath.h>
#include <Logger.h>

#include <TFile.h>
#include <TLorentzVector.h>


KFactorCorrection::KFactorCorrection(Dataset &dataset, Options const &options)
    : gLepBarePt_{dataset.Reader(), "GenPart_pt"},
      gLepBareEta_{dataset.Reader(), "GenPart_eta"},
      gLepBarePhi_{dataset.Reader(), "GenPart_phi"},
      gLepBareMass_{dataset.Reader(), "GenPart_mass"},
      gLepBareMomId_{dataset.Reader(), "GenPart_genPartIdxMother"} {

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

    for (int i = 0; i < gLepBarePt_.GetSize(); i++) {
      if (gLepBareMomId_[i] != 23)
        continue;

      TLorentzVector lepton;
      lepton.SetPtEtaPhiM(gLepBarePt_[i], gLepBareEta_[i], gLepBarePhi_[i],
        gLepBareMass_[i]);
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

