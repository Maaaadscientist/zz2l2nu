#include "Utils.h"
#define PI 3.141592654
#include <utility>

namespace utils
{
  double deltaPhi (TLorentzVector v1, TLorentzVector v2)
  {
    return deltaPhi (v1.Phi(), v2.Phi());
  }

  double deltaPhi (float phi1, float phi2)
  {
    double dPhi = phi1 - phi2; 
    if (dPhi > PI) dPhi -= 2*PI;
    if (dPhi <= -PI) dPhi += 2*PI; //So that dPhi is always between -pi and +pi.
    return dPhi;
  }

  double deltaR(TLorentzVector v1, TLorentzVector v2)
  {
    double dEta = v1.Eta()-v2.Eta();
    double dPhi = deltaPhi(v1,v2);
    return sqrt(pow(dEta,2) + pow(dPhi,2));
  }

  double deltaR(float eta1, float phi1, float eta2, float phi2)
  {
    double dEta = eta1 - eta2; 
    double dPhi = deltaPhi(phi1, phi2);
    return sqrt(pow(dEta,2) + pow(dPhi,2));
  }

  double getPhotonEnergy (double pT, double eta)
  {
    double theta = 2*atan(exp(-eta));
    return fabs(pT/sin(theta));
  }

  bool passVBFcuts(std::vector<TLorentzVectorWithIndex> selJets, TLorentzVector boson)
  {
    if(selJets.size()>=2){
      float etamin=0., etamax=0;
      if(selJets[0].Eta()>selJets[1].Eta()) {etamax = selJets[0].Eta(); etamin = selJets[1].Eta();}
      else {etamax = selJets[1].Eta(); etamin = selJets[0].Eta();}
      bool centralJetVeto = true;
      if(selJets.size()>2){
        for(unsigned int i = 2 ; i < selJets.size() ; i++){
          if(selJets[i].Eta()>etamin && selJets[i].Eta()<etamax) centralJetVeto = false;
        }
      }
      bool centralBoson = (boson.Eta()>etamin && boson.Eta()<etamax);
      bool passDeltaEta = (fabs(selJets[0].Eta() - selJets[1].Eta())>4);
      bool passMjj = ((selJets[0]+selJets[1]).M()>500);
      if(centralJetVeto && centralBoson && passDeltaEta && passMjj) return true;
    }
    return false;
  }

  bool passMetFilter(ULong64_t TrigMET, std::vector<std::pair<int, int> > & listMETFilter, bool isMC){
    //from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#ICHEP_2016
    std::vector<int> trigMET_position_all{14,  // primary vertex filter 
                                          8,  // beam halo filter 
                                          //see below: 3,  // HBHE noise filter 
                                          //see below: 4,  // HBHEiso noise filter 
                                          12,  // ECAL TP filter 
                                          15,  // ee badSC noise filter (not suggested for MC, but should be 100% according to David. To check) 
                                          24,  // badMuon (actually they went badPFMuon) 
                                          25   // badCharged hadron 
                                            };
    std::vector<int> trigMET_position_dataOnly{0, //duplicateMuons
                                               1, //badMuon
                                               2, //noBadMuon (i.e. no duplicate and no bad muons)
                                               3, //for the data, the condition on HBHE noise filter is inversed!
                                               4  //for the data, the condition on HBHEiso noise filter is inversed!
                                                };
    std::vector<int> trigMET_position_MCOnly{3, // HBHE noise filter (is also applied on data, but with the opposite value for the flag!)
                                             4  // HBHEiso noise filter (is also applied on data, but with the opposite value for the flag!)
                                              };

    bool passAllFilters = true;
    for(unsigned int i = 0; i < trigMET_position_all.size(); i++){
      if((TrigMET & (1<<trigMET_position_all[i]))) listMETFilter.push_back(std::make_pair(trigMET_position_all[i], 1));
      else{
        listMETFilter.push_back(std::make_pair(trigMET_position_all[i], 0));
        passAllFilters=false;
      }
    }
    if(!isMC){
      //For the 5 flags on data above, the logic is inversed: we ask the bin to be 0!
      for(unsigned int i = 0; i < trigMET_position_dataOnly.size(); i++){
        if(!(TrigMET & (1<<trigMET_position_dataOnly[i]))) listMETFilter.push_back(std::make_pair(trigMET_position_dataOnly[i], 1));
        else{
          listMETFilter.push_back(std::make_pair(trigMET_position_dataOnly[i], 0));
          passAllFilters=false;
        }
      }
    }
    else{
      for(unsigned int i = 0; i < trigMET_position_MCOnly.size(); i++){
        if((TrigMET & (1<<trigMET_position_MCOnly[i]))) listMETFilter.push_back(std::make_pair(trigMET_position_MCOnly[i], 1));
        else{
          listMETFilter.push_back(std::make_pair(trigMET_position_MCOnly[i], 0));
          passAllFilters=false;
        }
      }

    }
    return passAllFilters;
  }

}
