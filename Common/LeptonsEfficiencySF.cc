#include "LeptonsEfficiencySF.h"



namespace trigAndIDsfs
{
  std::pair<float,float> leptonEffSF(int particlePDGid, float pt, float eta, int cutType, int cutVersion){
    std::pair<float,float> eff; eff = std::make_pair(-1,-1);
    switch(particlePDGid){
      case 11 :
        switch (cutType){
          case llvvElecIdIso::ElecIdIso::Tight :
            switch (cutVersion){
              case CutVersion::CutSet::Moriond17Cut :
                if( eta >= -2.5 && eta < -2.0){
                  if( pt < 20.0){ eff.first=0.807; eff.second=0.018;
                  } else if( pt < 35.0){ eff.first=0.882; eff.second=0.010;
                  } else if( pt < 50.0){ eff.first=0.919; eff.second=0.009;
                  } else if( pt < 90.0){ eff.first=0.940; eff.second=0.007;
                  } else if( pt < 150.0){ eff.first=1.051; eff.second=0.022;
                  } else { eff.first=1.051; eff.second=0.106;
                  }
                }else if( eta >= -2.0 && eta < -1.57){
                  if( pt < 20.0){ eff.first=0.829; eff.second=0.018;
                  } else if( pt < 35.0){ eff.first=0.927; eff.second=0.018;
                  } else if( pt < 50.0){ eff.first=0.967; eff.second=0.007;
                  } else if( pt < 90.0){ eff.first=0.981; eff.second=0.006;
                  } else if( pt < 150.0){ eff.first=1.006; eff.second=0.022;
                  } else { eff.first=0.973; eff.second=0.030;
                  }
                }else if( eta >= -1.57 && eta < -1.444){
                  if( pt < 20.0){ eff.first=1.033; eff.second=0.106;
                  } else if( pt < 35.0){ eff.first=1.008; eff.second=0.110;
                  } else if( pt < 50.0){ eff.first=0.988; eff.second=0.017;
                  } else if( pt < 90.0){ eff.first=0.995; eff.second=0.024;
                  } else if( pt < 150.0){ eff.first=1.104; eff.second=0.050;
                  } else { eff.first=1.038; eff.second=0.075;
                }
                }else if( eta >= -1.444 && eta < -0.8){
                  if( pt < 20.0){ eff.first=1.008; eff.second=0.027;
                  } else if( pt < 35.0){ eff.first=0.972; eff.second=0.013;
                  } else if( pt < 50.0){ eff.first=0.975; eff.second=0.007;
                  } else if( pt < 90.0){ eff.first=0.972; eff.second=0.019;
                  } else if( pt < 150.0){ eff.first=0.989; eff.second=0.009;
                  } else { eff.first=0.982; eff.second=0.019;
                  }
                }else if( eta >= -0.8 && eta < 0.0){
                  if( pt < 20.0){ eff.first=0.941; eff.second=0.026;
                  } else if( pt < 35.0){ eff.first=0.953; eff.second=0.015;
                  } else if( pt < 50.0){ eff.first=0.953; eff.second=0.005;
                  } else if( pt < 90.0){ eff.first=0.953; eff.second=0.017;
                  } else if( pt < 150.0){ eff.first=0.975; eff.second=0.013;
                  } else { eff.first=0.982; eff.second=0.013;
                  }
                }else if( eta >= 0.0 && eta < 0.8){
                  if( pt < 20.0){ eff.first=0.946; eff.second=0.026;
                  } else if( pt < 35.0){ eff.first=0.982; eff.second=0.015;
                  } else if( pt < 50.0){ eff.first=0.980; eff.second=0.005;
                  } else if( pt < 90.0){ eff.first=0.978; eff.second=0.017;
                  } else if( pt < 150.0){ eff.first=1.012; eff.second=0.013;
                  } else { eff.first=1.021; eff.second=0.013;
                }
                }else if( eta >= 0.8 && eta < 1.444){
                  if( pt < 20.0){ eff.first=0.990; eff.second=0.027;
                  } else if( pt < 35.0){ eff.first=0.975; eff.second=0.013;
                  } else if( pt < 50.0){ eff.first=0.975; eff.second=0.007;
                  } else if( pt < 90.0){ eff.first=0.979; eff.second=0.019;
                  } else if( pt < 150.0){ eff.first=1.011; eff.second=0.009;
                  } else { eff.first=1.000; eff.second=0.019;
                  }
                }else if( eta >= 1.444 && eta < 1.57){
                  if( pt < 20.0){ eff.first=1.034; eff.second=0.106;
                  } else if( pt < 35.0){ eff.first=0.975; eff.second=0.110;
                  } else if( pt < 50.0){ eff.first=0.966; eff.second=0.017;
                  } else if( pt < 90.0){ eff.first=0.980; eff.second=0.024;
                  } else if( pt < 150.0){ eff.first=1.007; eff.second=0.049;
                  } else { eff.first=0.884; eff.second=0.076;
                  }
                }else if( eta >= 1.57 && eta < 2.0){
                  if( pt < 20.0){ eff.first=0.827; eff.second=0.018;
                  } else if( pt < 35.0){ eff.first=0.909; eff.second=0.018;
                  } else if( pt < 50.0){ eff.first=0.958; eff.second=0.007;
                  } else if( pt < 90.0){ eff.first=0.969; eff.second=0.006;
                  } else if( pt < 150.0){ eff.first=0.988; eff.second=0.022;
                  } else { eff.first=0.979; eff.second=0.030;
                  }
                }else {
                  if( pt < 20.0){ eff.first=0.797; eff.second=0.018;
                  } else if( pt < 35.0){ eff.first=0.863; eff.second=0.010;
                  } else if( pt < 50.0){ eff.first=0.908; eff.second=0.009;
                  } else if( pt < 90.0){ eff.first=0.938; eff.second=0.007;
                  } else if( pt < 150.0){ eff.first=1.021; eff.second=0.022;
                  } else { eff.first=1.048; eff.second=0.106;
                  }
                }
                break;
              }
            break;
          }
        break;
      }
      return eff;
    }
}

