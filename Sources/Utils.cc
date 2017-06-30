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

}
