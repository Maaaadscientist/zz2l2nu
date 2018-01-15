#include "LeptonsEfficiencySF.h"



namespace trigAndIDsfs
{
  std::pair<float,float> leptonEffSF(int particlePDGid, float pt, float eta, int cutType, int cutVersion){
    std::pair<float,float> eff; eff = std::make_pair(-1,-1);
    switch(particlePDGid){
      case 11 :
        switch (cutType){
          case llvvElecRecoIdIso::ElecRecoIdIso::Tight :
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
            case llvvElecRecoIdIso::ElecRecoIdIso::Reco :
              switch (cutVersion){
                case CutVersion::CutSet::Moriond17Cut :
                  if (eta >= -2.5 && eta < -2.45){
                    eff.first = 1.318 ; eff.second = 0.018;
                  }else if (eta >= -2.450 && eta < -2.400){
                    eff.first = 1.114 ; eff.second = 0.011;
                  }else if (eta >= -2.400 && eta < -2.300){
                    eff.first = 1.025 ; eff.second = 0.008;
                  }else if (eta >= -2.300 && eta < -2.200){
                    eff.first = 1.014 ; eff.second = 0.007;
                  }else if (eta >= -2.200 && eta < -2.000){
                    eff.first = 1.007 ; eff.second = 0.004;
                  }else if (eta >= -2.000 && eta < -1.800){
                    eff.first = 0.995 ; eff.second = 0.006;
                  }else if (eta >= -1.800 && eta < -1.630){
                    eff.first = 0.995 ; eff.second = 0.005;
                  }else if (eta >= -1.630 && eta < -1.566){
                    eff.first = 0.992 ; eff.second = 0.006;
                  }else if (eta >= -1.566 && eta < -1.444){
                    eff.first = 0.963 ; eff.second = 0.026;
                  }else if (eta >= -1.444 && eta < -1.200){
                    eff.first = 0.990 ; eff.second = 0.004;
                  }else if (eta >= -1.200 && eta < -1.000){
                    eff.first = 0.986 ; eff.second = 0.005;
                  }else if (eta >= -1.000 && eta < -0.600){
                    eff.first = 0.982 ; eff.second = 0.003;
                  }else if (eta >= -0.600 && eta < -0.400){
                    eff.first = 0.985 ; eff.second = 0.006;
                  }else if (eta >= -0.400 && eta < -0.200){
                    eff.first = 0.982 ; eff.second = 0.006;
                  }else if (eta >= -0.200 && eta < 0.000){
                    eff.first = 0.980 ; eff.second = 0.005;
                  }else if (eta >= 0.000 && eta < 0.200){
                    eff.first = 0.985 ; eff.second = 0.005;
                  }else if (eta >= 0.200 && eta < 0.400){
                    eff.first = 0.989 ; eff.second = 0.006;
                  }else if (eta >= 0.400 && eta < 0.600){
                    eff.first = 0.988 ; eff.second = 0.006;
                  }else if (eta >= 0.600 && eta < 1.000){
                    eff.first = 0.988 ; eff.second = 0.003;
                  }else if (eta >= 1.000 && eta < 1.200){
                    eff.first = 0.988 ; eff.second = 0.005;
                  }else if (eta >= 1.200 && eta < 1.444){
                    eff.first = 0.988 ; eff.second = 0.004;
                  }else if (eta >= 1.444 && eta < 1.566){
                    eff.first = 0.968 ; eff.second = 0.026;
                  }else if (eta >= 1.566 && eta < 1.630){
                    eff.first = 0.990 ; eff.second = 0.006;
                  }else if (eta >= 1.630 && eta < 1.800){
                    eff.first = 0.993 ; eff.second = 0.005;
                  }else if (eta >= 1.800 && eta < 2.000){
                    eff.first = 0.992 ; eff.second = 0.006;
                  }else if (eta >= 2.000 && eta < 2.200){
                    eff.first = 0.998 ; eff.second = 0.004;
                  }else if (eta >= 2.200 && eta < 2.300){
                    eff.first = 1.001 ; eff.second = 0.007;
                  }else if (eta >= 2.300 && eta < 2.400){
                    eff.first = 0.990 ; eff.second = 0.008;
                  }else if (eta >= 2.400 && eta < 2.450){
                    eff.first = 0.971 ; eff.second = 0.011;
                  }else if (eta >= 2.450 && eta < 2.500){
                    eff.first = 0.907 ; eff.second = 0.018;
                  }
                  break;
              }
                break;
          }
        break;
        case 13:
          switch (cutType){
            case llvvMuonIdIso::MuonIdIso::Tight :
              switch (cutVersion){
                case CutVersion::CutSet::Moriond17CutRunGH :
                  if( std::abs(eta) >= 0.0 && std::abs(eta) < 0.9){
                    if( pt < 25.0){ eff.first=0.9932; eff.second=0.9932;
                    } else if( pt < 30.0){ eff.first=0.9870; eff.second=0.9870;
                    } else if( pt < 40.0){ eff.first=0.9876; eff.second=0.9876;
                    } else if( pt < 50.0){ eff.first=0.9898; eff.second=0.9898;
                    } else if( pt < 60.0){ eff.first=0.9847; eff.second=0.9847;
                    } else { eff.first=0.9914; eff.second=0.9914;
                    }
                  } else if( std::abs(eta) >= 0.9 && std::abs(eta) < 1.2){
                    if( pt < 25.0){ eff.first=0.9856; eff.second=0.9856;
                    } else if( pt < 30.0){ eff.first=0.9847; eff.second=0.9847;
                    } else if( pt < 40.0){ eff.first=0.9839; eff.second=0.9839;
                    } else if( pt < 50.0){ eff.first=0.9833; eff.second=0.9833;
                    } else if( pt < 60.0){ eff.first=0.9806; eff.second=0.9806;
                    } else { eff.first=0.9839; eff.second=0.9839;
                    }
                  } else if( std::abs(eta) >= 1.2 && std::abs(eta) < 2.1){
                    if( pt < 25.0){ eff.first=0.9909; eff.second=0.9909;
                    } else if( pt < 30.0){ eff.first=0.9909; eff.second=0.9909;
                    } else if( pt < 40.0){ eff.first=0.9921; eff.second=0.9921;
                    } else if( pt < 50.0){ eff.first=0.9938; eff.second=0.9938;
                    } else if( pt < 60.0){ eff.first=0.9857; eff.second=0.9857;
                    } else { eff.first=0.9886; eff.second=0.9886;
                    }
                  } else if( std::abs(eta) >= 2.1 && std::abs(eta) < 2.4){
                    if( pt < 25.0){ eff.first=0.9815; eff.second=0.9815;
                    } else if( pt < 30.0){ eff.first=0.9791; eff.second=0.9791;
                    } else if( pt < 40.0){ eff.first=0.9715; eff.second=0.9715;
                    } else if( pt < 50.0){ eff.first=0.9748; eff.second=0.9748;
                    } else if( pt < 60.0){ eff.first=0.9677; eff.second=0.9677;
                    } else { eff.first=0.9632; eff.second=0.9632;
                    }
                  }
                  break;
              }
              break;
            case llvvMuonIdIso::MuonIdIso::TightIso :
              switch (cutVersion){
                case CutVersion::CutSet::Moriond17CutRunGH :
                  if( std::abs(eta) >= 0.0 && std::abs(eta) < 0.9){
                    if( pt < 25.0){ eff.first=0.9811; eff.second=0.9811;
                    } else if( pt < 30.0){ eff.first=0.9928; eff.second=0.9928;
                    } else if( pt < 40.0){ eff.first=0.9935; eff.second=0.9935;
                    } else if( pt < 50.0){ eff.first=0.9952; eff.second=0.9952;
                    } else if( pt < 60.0){ eff.first=0.9967; eff.second=0.9967;
                    } else { eff.first=0.9991; eff.second=0.9991;
                    }
                  } else if( std::abs(eta) >= 0.9 && std::abs(eta) < 1.2){
                    if( pt < 25.0){ eff.first=0.9977; eff.second=0.9977;
                    } else if( pt < 30.0){ eff.first=0.9997; eff.second=0.9997;
                    } else if( pt < 40.0){ eff.first=0.9996; eff.second=0.9996;
                    } else if( pt < 50.0){ eff.first=0.9983; eff.second=0.9983;
                    } else if( pt < 60.0){ eff.first=0.9989; eff.second=0.9989;
                    } else { eff.first=0.9989; eff.second=0.9989;
                    }
                  } else if( std::abs(eta) >= 1.2 && std::abs(eta) < 2.1){
                    if( pt < 25.0){ eff.first=0.9937; eff.second=0.9937;
                    } else if( pt < 30.0){ eff.first=0.9980; eff.second=0.9980;
                    } else if( pt < 40.0){ eff.first=0.9990; eff.second=0.9990;
                    } else if( pt < 50.0){ eff.first=0.9986; eff.second=0.9986;
                    } else if( pt < 60.0){ eff.first=0.9988; eff.second=0.9988;
                    } else { eff.first=0.9995; eff.second=0.9995;
                    }
                  } else if( std::abs(eta) >= 2.1 && std::abs(eta) < 2.4){
                    if( pt < 25.0){ eff.first=0.9942; eff.second=0.9942;
                    } else if( pt < 30.0){ eff.first=0.9991; eff.second=0.9991;
                    } else if( pt < 40.0){ eff.first=1.0001; eff.second=1.0001;
                    } else if( pt < 50.0){ eff.first=1.0002; eff.second=1.0002;
                    } else if( pt < 60.0){ eff.first=1.0000; eff.second=1.0000;
                    } else { eff.first=1.0015; eff.second=1.0015;
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

