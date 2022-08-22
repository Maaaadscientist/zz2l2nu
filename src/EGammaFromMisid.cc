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
  bool e1_is_probe = false, e2_is_probe = false;

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
    return false;
  } else {
    return false;
  }

  auto l1 = &electrons[0], l2 = &electrons[1];

  e1_is_probe = CheckProbe(l1);
  e2_is_probe = CheckProbe(l2);

  if (!(e1_is_probe || e2_is_probe)) {
    return false;
  }

  eventCat_ = int(eventCat);
  TLorentzVector const p4tot = l1->p4 + l2->p4;
  totPt_ = p4tot.Pt();
  totEta_ = p4tot.Eta();
  totPhi_ = p4tot.Phi();
  totMass_ = p4tot.M();

  if (std::abs(p4tot.M() - kNominalMZ_) > 10)
    return false;

  numPVGood_ = *srcNumPVGood_;

  if (e1_is_probe) {
    probePt_ = l1->p4.Pt();
    probeEta_ = l1->p4.Eta();

    if (storeMoreVariables_)
      FillMoreVariables();

    FillTree();
  }
  if (e2_is_probe) {
    probePt_ = l2->p4.Pt();
    probeEta_ = l2->p4.Eta();

    if (storeMoreVariables_)
      FillMoreVariables();

    FillTree();
  }

  return true;
}

bool EGammaFromMisid::CheckProbe(std::variant<Electron const *, Photon const *> particle) {
  auto e = std::get<Electron const *>(particle);
  return e->p4.Pt() > minPtLL_;
  // TODO: implement the photon check
}

void EGammaFromMisid::FillMoreVariables() {

  event_ = *srcEvent_;
}
