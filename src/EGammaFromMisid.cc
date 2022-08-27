#include <EGammaFromMisid.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <TLorentzVector.h>
#include <TVector2.h>

#include <Utils.h>


namespace po = boost::program_options;


EGammaFromMisid::EGammaFromMisid(Options const &options, Dataset &dataset)
    : EventTrees{options, dataset},
      storeMoreVariables_{options.Exists("more-vars")},
      // triggerFilter_{dataset, options, &runSampler_},
      srcEvent_{dataset.Reader(), "event"},
      photonBuilder_{dataset},
      // photonFilter_{dataset, options},
      photonPrescales_{dataset, options},
      photonWeight_{dataset, options, &photonBuilder_},
      srcNumPVGood_{dataset.Reader(), "PV_npvsGood"} {

  photonBuilder_.EnableCleaning({&muonBuilder_, &electronBuilder_});

  weightCollector_.Add(&photonWeight_);

  CreateWeightBranches();

  AddBranch("event_cat", &eventCat_);
  AddBranch("tot_mass", &totMass_);
  AddBranch("num_pv_good", &numPVGood_);
  AddBranch("probe_pt", &probePt_);
  AddBranch("probe_eta", &probeEta_);

  if (storeMoreVariables_) {
    AddBranch("event", &event_);
  }
}


po::options_description EGammaFromMisid::OptionsDescription() {
  auto optionsDescription = AnalysisCommon::OptionsDescription();
  optionsDescription.add_options()
    ("more-vars", "Store additional variables");
  return optionsDescription;
}


bool EGammaFromMisid::ProcessEvent() {
  if (not ApplyCommonFilters())
    return false;

  // int ee_count = 0;
  Electron const *e0, *e1;
  Electron const *e;
  Photon const *photon;
  bool e0_is_probe = false, e1_is_probe = false;

  auto const &electrons = electronBuilder_.GetTight();
  auto const &photons = photonBuilder_.Get(); // baseline photons

  auto const &looseElectrons = electronBuilder_.GetLoose();
  auto const &looseMuons = muonBuilder_.GetLoose();

  if (looseElectrons.size() != electrons.size())
    return false;

  if (looseMuons.size() > 0)
    return false;

  if (isotrkBuilder_.Get().size() > 0)
    return false;

  EventCat eventCat;
  if (electrons.size() == 2 && photons.size() == 0) {
    eventCat = EventCat::kEE;
  } else if (electrons.size() == 1 && photons.size() == 1) {
    eventCat = EventCat::kEGamma;
  } else {
    return false;
  }

  TLorentzVector p4tot;
  switch (eventCat) {
    case EventCat::kEE:
      e0 = &electrons[0];
      e1 = &electrons[1];

      e0_is_probe = CheckProbe(e0);
      e1_is_probe = CheckProbe(e1);

      if (!(e0_is_probe || e1_is_probe)) {
        return false;
      }

      p4tot = e0->p4 + e1->p4;
      break;

    case EventCat::kEGamma:
      e = &electrons[0];
      photon = &photons[0];
      if (!CheckProbe(photon)) {
        return false;
      }
      p4tot = e->p4 + photon->p4;
      break;
  }

  totPt_ = p4tot.Pt();
  totEta_ = p4tot.Eta();
  totPhi_ = p4tot.Phi();
  totMass_ = p4tot.M();

  if (std::abs(p4tot.M() - kNominalMZ_) > 10)
    return false;

  eventCat_ = int(eventCat);
  numPVGood_ = *srcNumPVGood_;

  switch (eventCat) {
    case EventCat::kEE:
      if (e0_is_probe) {
        probePt_ = e0->p4.Pt();
        probeEta_ = e0->p4.Eta();

        if (storeMoreVariables_)
          FillMoreVariables();

        FillTree();
      }
      if (e1_is_probe) {
        probePt_ = e1->p4.Pt();
        probeEta_ = e1->p4.Eta();

        if (storeMoreVariables_)
          FillMoreVariables();

        FillTree();
      }
      break;

    case EventCat::kEGamma:
      probePt_ = photon->p4.Pt();
      probeEta_ = photon->p4.Eta();

      if (storeMoreVariables_)
        FillMoreVariables();

      FillTree();
      break;
  }

  return true;
}

bool EGammaFromMisid::CheckProbe(std::variant<Electron const *, Photon const *> particle) {
  if (std::holds_alternative<Photon const *>(particle)) {
      auto photon = std::get<Photon const *>(particle);
      if (!photon->passElecVeto) {
        return false;
      }
      if (photon->sieie <= 0.001) {
        return false;
      }
      return photon->p4.Pt() > minPtLL_;
  } else {
    auto e = std::get<Electron const *>(particle);
    return e->p4.Pt() > minPtLL_;
  }
}

void EGammaFromMisid::FillMoreVariables() {

  event_ = *srcEvent_;
}
