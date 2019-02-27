#include <BTagger.h>

namespace btagger
{

  BTagCalibrationReader loadCalibrationReader()
  {
    std::string calib_file = (std::string)getenv("HZZ2L2NU_BASE") + "/data/efficiencyTables/CSVv2_Moriond17_B_H.csv";
    BTagCalibration calib("", calib_file);
    BTagEntry::OperatingPoint wp = BTagEntry::OP_LOOSE;
    BTagCalibrationReader _btag_calibration_reader;
    _btag_calibration_reader = BTagCalibrationReader(wp, "central", {"up", "down"});
    _btag_calibration_reader.load(calib, BTagEntry::FLAV_B, "mujets");
    _btag_calibration_reader.load(calib, BTagEntry::FLAV_C, "mujets");
    _btag_calibration_reader.load(calib, BTagEntry::FLAV_UDSG, "incl");
    return _btag_calibration_reader;
  }

  utils::tables loadEffTables()
  {
    utils::tables theMapOfTables;
    std::vector<std::string> btagWP = {"loose","medium","tight"}; // Actually, we could load the maps for our WP only. But the gain in time is very small.
    std::vector<std::string> btagFlavor = {"b","c","udsg"};
    std::string btagEffTablePath = (std::string)getenv("HZZ2L2NU_BASE") + "/data/efficiencyTables/";
    for(auto i : btagWP){
      for(auto j : btagFlavor){
        theMapOfTables.insert( std::pair<std::string, utils::table> (i+"-"+j,utils::table(btagEffTablePath+"btag-"+i+"-"+j+".txt")) );
      }
    }
    return theMapOfTables;
  }

  
  double apply_sf(std::vector<TLorentzVectorWithIndex> selCentralJets, std::vector<double> btags, std::vector<float> *JetAk04HadFlav, const utils::tables &tab, BTagCalibrationReader _btag_calibration_reader, TString systName)
  {
    std::string workingPoint = "loose"; // If you change this, make sure it's consistent in the whole code.
    double bjetCut = 0.5426; // Same here.
    double btagWeight = 1.;
    for(unsigned int i = 0 ; i < selCentralJets.size() ; i ++){
      std::string tag = "";
      BTagEntry::JetFlavor flavor;
      if(fabs(JetAk04HadFlav->at(selCentralJets.at(i).GetIndex()))==5) {flavor = BTagEntry::FLAV_B; tag = "b";}
      else if(fabs(JetAk04HadFlav->at(selCentralJets.at(i).GetIndex()))==4) {flavor = BTagEntry::FLAV_C; tag = "c";}
      else {flavor = BTagEntry::FLAV_UDSG; tag = "udsg";}
      double eff = tab.at(workingPoint+"-"+tag).getEfficiency(selCentralJets[i].Pt(), selCentralJets[i].Eta());
      bool tagged = btags.at(i) > bjetCut;
      double sf = 1.; 
      if (systName == "btag_up"){
        sf = _btag_calibration_reader.eval_auto_bounds("up", flavor, std::abs(selCentralJets[i].Eta()), selCentralJets[i].Pt());
      }
      else if (systName == "btag_down"){
        sf = _btag_calibration_reader.eval_auto_bounds("down", flavor, std::abs(selCentralJets[i].Eta()), selCentralJets[i].Pt());
      }
      else{
        sf = _btag_calibration_reader.eval_auto_bounds("central", flavor, std::abs(selCentralJets[i].Eta()), selCentralJets[i].Pt());
      }
      btagWeight *= tagged ? sf : (1 - sf * eff) / (1 - eff);
      //std::cout << "Jet nb " << i << " : tag = " << tag << " ; (pT,eta) = (" << selCentralJets[i].Pt() << "," << selCentralJets[i].Eta() << ") ; btag = " << btags.at(i) << " ; eff = " << eff << " ; sf = " << sf << std::endl;
    }
    //std::cout << "btagWeight: " << btagWeight << std::endl;
    return btagWeight;
  }

}
