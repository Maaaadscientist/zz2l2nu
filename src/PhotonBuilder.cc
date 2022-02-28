#include <PhotonBuilder.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <Utils.h>


PhotonBuilder::PhotonBuilder(Dataset &dataset)
    : CollectionBuilder{dataset.Reader()},
      minPt_{20.},
      isSim_{dataset.Info().IsSimulation()},
      srcPt_{dataset.Reader(), "Photon_pt"},
      srcEta_{dataset.Reader(), "Photon_eta"},
      srcPhi_{dataset.Reader(), "Photon_phi"},
      srcIsEtaScEb_{dataset.Reader(), "Photon_isScEtaEB"},
      srcPixelSeed_{dataset.Reader(), "Photon_pixelSeed"},
      srcElecronVeto_{dataset.Reader(), "Photon_electronVeto"},
      srcR9_{dataset.Reader(), "Photon_r9"},
      srcSieie_{dataset.Reader(), "Photon_sieie"}
{
  srcId_.reset(
      new TTreeReaderArray<int>(dataset.Reader(), "Photon_cutBased"));

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
    passId = (srcId_->At(i) >= 3);
    //std::cout<< "photon No."<<i<<" : pt="<<srcPt_[i]<< " eta="<<srcEta_[i] << " phi="<<srcPhi_[i]<<" passID="<<(srcId_->At(i) >= 3)<<std::endl;
    
    if (srcPt_[i] < minPt_ or not passId)
      continue;
    
    photon.p4.SetPtEtaPhiM(srcPt_[i], srcEta_[i], srcPhi_[i], 0.);

    // Only consider photons in the barrel except for Njet >= 2
    photon.isEB = srcIsEtaScEb_[i];

    // Conversion safe electron veto and pixel seed vetos, stored for event veto
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
    photon.passElecVeto = (!srcPixelSeed_[i] && srcElecronVeto_[i]);

    // Only consider events with photon sieie > 0.001
    photon.sieie = srcSieie_[i];

    // Extra variables
    photon.r9 = srcR9_[i];

    // Perform angular cleaning
    //std::cout<<"photon No."<<i<<" isEB="<<srcIsEtaScEb_[i] << " passEleVeto="<<(!srcPixelSeed_[i] && srcElecronVeto_[i])<< " sieie" << srcSieie_[i] << " r9="<<srcR9_[i]<< std::endl;
    if (IsDuplicate(photon.p4, 0.1)){
      //std::cout<<"deltaR cleaning < 0.1"<<std::endl;
      continue;
    }
    if (photon.p4.Eta()<= 1.58 && photon.p4.Eta()>= 1.48 && photon.p4.Phi()>= -0.78 && photon.p4.Phi() <= -0.55  ){
      continue;
    }
    if (abs(photon.p4.Eta())> 1.58 && ( abs(photon.p4.Phi()) > 3.14159* 11/12 || abs(photon.p4.Phi()) < 3.14159/12)) continue;
    photons_.emplace_back(photon);
  }

  // Make sure the collection is ordered in pt
  std::sort(photons_.begin(), photons_.end(), PtOrdered);
}

