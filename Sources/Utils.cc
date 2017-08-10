#include "../Includes/Utils.h"
#define PI 3.141592654

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

}
