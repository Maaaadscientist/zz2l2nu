#include <PhotonBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


PhotonBuilder::PhotonBuilder(Dataset &dataset, Options const &options)
    : CollectionBuilder{dataset.Reader()},
      minPt_{55.},
      isSim_{dataset.Info().IsSimulation()},
      srcPt_{dataset.Reader(), "Photon_pt"},
      srcEta_{dataset.Reader(), "Photon_eta"},
      srcPhi_{dataset.Reader(), "Photon_phi"},
      srcIsEtaScEb_{dataset.Reader(), "Photon_isScEtaEB"} {

  std::string year = Options::NodeAs<std::string>(options.GetConfig(), {"period"});
  if (year == "2016") {
    idBranchName_ = "Photon_cutBased";
    isIdBitmap_ = false;
  }
  else {
    idBranchName_ = "Photon_cutBasedBitmap";
    isIdBitmap_ = true;
  }
  srcId_.reset(new TTreeReaderArray<int>(dataset.Reader(), idBranchName_.c_str()));

  if (isSim_) {
    srcGenPt_.reset(
      new TTreeReaderArray<float>(dataset.Reader(), "GenPart_pt"));
    srcGenEta_.reset(
      new TTreeReaderArray<float>(dataset.Reader(), "GenPart_eta"));
    srcGenPhi_.reset(
      new TTreeReaderArray<float>(dataset.Reader(), "GenPart_phi"));
    srcPhotonGenPartIndex_.reset(
      new TTreeReaderArray<int>(dataset.Reader(), "Photon_genPartIdx"));
    srcFlavour_.reset(
      new TTreeReaderArray<UChar_t>(dataset.Reader(), "Photon_genPartFlav"));
  }
}


std::vector<Photon> const &PhotonBuilder::Get() const {
  Update();
  return photons_;
}


void PhotonBuilder::Build() const {
  photons_.clear();

  for (unsigned i = 0; i < srcPt_.GetSize(); ++i) {
    Photon photon;
    // Gen particle matching
    if (isSim_) {
      if (srcPhotonGenPartIndex_->At(i) >= 0){
        photon.genP4.SetPtEtaPhiM(srcGenPt_->At(srcPhotonGenPartIndex_->At(i)),
            srcGenEta_->At(srcPhotonGenPartIndex_->At(i)),
            srcGenPhi_->At(srcPhotonGenPartIndex_->At(i)), 0);
        // Photon flavour. See NanoAOD documentation for the conventions.
        if (srcFlavour_->At(i) == 1)
          photon.flavour = Photon::Origin::PromptPhoton;
        else if (srcFlavour_->At(i) == 11)
          photon.flavour = Photon::Origin::PromptElectron;
        else if (srcFlavour_->At(i) == 0)
          photon.flavour = Photon::Origin::Unmatched;
      }
      else {  // if no gen-level photon matched
        photon.genP4.SetPtEtaPhiM(0, 0, 0, 0);
        photon.flavour = Photon::Origin::Unmatched;
      }

      if (IsDuplicate(photon.genP4, 0.1))
        continue;
    }

    // Tight ID
    bool passId = false;
    if (isIdBitmap_) {
      passId = (srcId_->At(i) & (1 << 2));
    }
    else {
      passId = (srcId_->At(i) >= 3);
    }

    
    if (srcPt_[i] < minPt_ or not passId)
      continue;
    
    // Only consider photons in the barrel
    if (!srcIsEtaScEb_[i])
      continue;

    // Conversion safe electron veto and pixel seed veto could be added:
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2

    photon.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);

    // Perform angular cleaning
    if (IsDuplicate(photon.p4, 0.1))
      continue;

    photons_.emplace_back(photon);
  }

  // Make sure the collection is ordered in pt
  std::sort(photons_.begin(), photons_.end(), PtOrdered);
}

