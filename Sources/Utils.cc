#include "../Includes/Utils.h"
#include "../Includes/LooperMain.h"
#define PI 3.141592654
#include <utility>

namespace utils
{
  double deltaPhi (TLorentzVector v1, TLorentzVector v2)
  {
    double dPhi = v1.Phi()-v2.Phi();
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

  double getPhotonEnergy (double pT, double eta)
  {
    double theta = 2*atan(exp(-eta));
    return fabs(pT/sin(theta));
  }

  bool passVBFcuts(std::vector<TLorentzVector> selJets, TLorentzVector boson)
  {
    if(selJets.size()>=2){
      float etamin=0., etamax=0;
      if(selJets[0].Eta()>selJets[1].Eta()) {etamax = selJets[0].Eta(); etamin = selJets[1].Eta();}
      else {etamax = selJets[1].Eta(); etamin = selJets[0].Eta();}
      bool centralJetVeto = true;
      if(selJets.size()>2){
        for(int i = 2 ; i < selJets.size() ; i++){
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

  int trigDoubleMu[4] = {8,9,10,11};
  int trigSingleMu[4] = {10,11,15,16};
  int trigDoubleE[2] = {12,13};
  int trigHighPtE[1] = {17};//Located in DoubleElectron
  int trigSingleE[2] = {11,12};
  int trigEMu[2] = {0,3}; //all DZ paths are still missing, as well as Mu12Ele23
  //int trigSinglePhoton[12] = {1,21,0,20,32,31,30,29,28,27,26,25}; //ordered by decreasing Pt, see below. //FIXME FIXME FIXME FIXME FIXME
  int trigSinglePhoton[10] = {1,0,30,29,28,27,26,25,24,23}; //ordered by decreasing Pt, see below.


  //The pt threshold corresponding to 
  //int lowThresholdPt[12] = {300,300,250,250,165,120,90,75,50,36,30,22}; //ordered by Pt to avoid double counting. //FIXME FIXME FIXME FIXME FIXME
  //int highThresholdPt[12] = {9999,9999,300,300,250,165,120,90,75,50,36,30}; //ordered by Pt to avoid double counting. //FIXME FIXME FIXME FIXME FIXME
  int lowThresholdPt[10] = {300,250,165,120,90,75,50,36,30,22}; //ordered by Pt to avoid double counting.
  int highThresholdPt[10] = {9999,300,250,165,120,90,75,50,36,30}; //ordered by Pt to avoid double counting.

  int passTrigger(int trig, ULong64_t TrigHltDiMu, ULong64_t TrigHltMu, ULong64_t TrigHltDiEl, ULong64_t TrigHltEl, ULong64_t TrigHltElMu, ULong64_t TrigHltPhot, vector<unsigned int> *TrigHltDiMu_prescale, vector<unsigned int> *TrigHltMu_prescale, vector<unsigned int> *TrigHltDiEl_prescale, vector<unsigned int> *TrigHltEl_prescale, vector<unsigned int> *TrigHltElMu_prescale, vector<unsigned int> *TrigHltPhot_prescale, double selectedPhotonPt){
    std::vector<std::vector<int> > trigList(Ntrig);
    trigList[DoubleMu].insert(trigList[DoubleMu].end(),trigDoubleMu,trigDoubleMu+(sizeof(trigDoubleMu)/sizeof(trigDoubleMu[0])));
    trigList[SingleMu].insert(trigList[SingleMu].end(),trigSingleMu,trigSingleMu+(sizeof(trigSingleMu)/sizeof(trigSingleMu[0])));
    trigList[DoubleE].insert(trigList[DoubleE].end(),trigDoubleE,trigDoubleE+(sizeof(trigDoubleE)/sizeof(trigDoubleE[0])));
    trigList[HighPtE].insert(trigList[HighPtE].end(),trigHighPtE,trigHighPtE+(sizeof(trigHighPtE)/sizeof(trigHighPtE[0])));
    trigList[SingleE].insert(trigList[SingleE].end(),trigSingleE,trigSingleE+(sizeof(trigSingleE)/sizeof(trigSingleE[0])));
    trigList[EMu].insert(trigList[EMu].end(),trigEMu,trigEMu+(sizeof(trigEMu)/sizeof(trigEMu[0])));
    trigList[SinglePhoton].insert(trigList[SinglePhoton].end(),trigSinglePhoton,trigSinglePhoton+(sizeof(trigSinglePhoton)/sizeof(trigSinglePhoton[0])));
    switch(trig){
      case DoubleMu:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return TrigHltDiMu_prescale->at(trigList.at(DoubleMu).at(i));
        break;
      case SingleMu:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return TrigHltMu_prescale->at(trigList.at(SingleMu).at(i));
        break;
      case DoubleE: //Includes also HighPtE
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(DoubleE).at(i));
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(HighPtE).at(i)); //Accepted also for HighPtE
        break;
      case SingleE:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return TrigHltEl_prescale->at(trigList.at(SingleE).at(i));
        break;
      case EMu:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1<<trigList.at(EMu).at(i))) return TrigHltElMu_prescale->at(trigList.at(EMu).at(i));
        break;
      case SinglePhoton:
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1<<trigList.at(EMu).at(i))) return 0;
        for(unsigned int i = 0 ; i < trigList[SinglePhoton].size() ; i++){
          if(TrigHltPhot & (1<<trigList.at(SinglePhoton).at(i))){
            if(selectedPhotonPt){
              if((selectedPhotonPt >= lowThresholdPt[i]) && (selectedPhotonPt <= highThresholdPt[i]+10)) return TrigHltPhot_prescale->at(trigList.at(SinglePhoton).at(i));
            }
            else return true;
          }
        }
        break;
      case Ntrig://In this case (used for MC), take if any trigger passed
        for(unsigned int i = 0 ; i < trigList[SinglePhoton].size() ; i++){
          if(TrigHltPhot & (1<<trigList.at(SinglePhoton).at(i))){
            if(selectedPhotonPt){
              if((selectedPhotonPt >= lowThresholdPt[i]) && (selectedPhotonPt <= highThresholdPt[i]+10)) return TrigHltPhot_prescale->at(trigList.at(SinglePhoton).at(i));
            }
            else return true;
          }
        }
        for(unsigned int i = 0 ; i < trigList[EMu].size() ; i++)  if(TrigHltElMu & (1<<trigList.at(EMu).at(i))) return TrigHltElMu_prescale->at(trigList.at(EMu).at(i));
        for(unsigned int i = 0 ; i < trigList[SingleE].size() ; i++)  if(TrigHltEl & (1<<trigList.at(SingleE).at(i))) return TrigHltEl_prescale->at(trigList.at(SingleE).at(i));
        for(unsigned int i = 0 ; i < trigList[HighPtE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(HighPtE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(HighPtE).at(i));
        for(unsigned int i = 0 ; i < trigList[DoubleE].size() ; i++)  if(TrigHltDiEl & (1<<trigList.at(DoubleE).at(i))) return TrigHltDiEl_prescale->at(trigList.at(DoubleE).at(i));
        for(unsigned int i = 0 ; i < trigList[SingleMu].size() ; i++)  if(TrigHltMu & (1<<trigList.at(SingleMu).at(i))) return TrigHltMu_prescale->at(trigList.at(SingleMu).at(i));
        for(unsigned int i = 0 ; i < trigList[DoubleMu].size() ; i++)  if(TrigHltDiMu & (1<<trigList.at(DoubleMu).at(i))) return TrigHltDiMu_prescale->at(trigList.at(DoubleMu).at(i));
        break;
      default:
        return 0;
    }
    return 0; //If nothing found.
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
