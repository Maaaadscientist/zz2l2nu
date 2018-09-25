#include "BTagger.h"

namespace btagger
{

  BTagCalibrationReader loadCalibrationReader()
  {
    std::string calib_file = (std::string)getenv("CMSSW_BASE")+"/src/shears/HZZ2l2nu/data/efficiencyTables/CSVv2_Moriond17_B_H.csv";
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
    std::vector<std::string> btagWP = {"loose","medium","tight"};
    std::vector<std::string> btagFlavor = {"b","c","udsg"};
    std::string btagEffTablePath = (std::string)getenv("CMSSW_BASE")+(std::string)"/src/shears/HZZ2l2nu/data/efficiencyTables/";
    for(auto i : btagWP){
      for(auto j : btagFlavor){
        theMapOfTables.insert( std::pair<std::string, utils::table> (i+"-"+j,utils::table(btagEffTablePath+"btag-"+i+"-"+j+".txt")) );
      }
    }
    return theMapOfTables;
  }

  bool fill_eff(std::vector<TLorentzVectorWithIndex> selCentralJets, std::vector<double> btags, std::vector<float> *JetAk04HadFlav, double weight, SmartSelectionMonitor_hzz & mon)
  {
    for(unsigned int i = 0 ; i < selCentralJets.size() ; i ++){
      std::string tag = "";
      if(fabs(JetAk04HadFlav->at(selCentralJets.at(i).GetIndex()))==5) tag = "bjet";
      else if(fabs(JetAk04HadFlav->at(selCentralJets.at(i).GetIndex()))==4) tag = "cjet";
      else tag = "udsgjet";
      bool tagged_loose = btags.at(i) > 0.5426;
      bool tagged_medium = btags.at(i) > 0.8484;
      bool tagged_tight = btags.at(i) > 0.9535;
      mon.fillHisto("btagEff","den_"+tag,selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
      if(tagged_loose) mon.fillHisto("btagEff","num_"+tag+"_tagged_loose",selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
      if(tagged_medium) mon.fillHisto("btagEff","num_"+tag+"_tagged_medium",selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
      if(tagged_tight) mon.fillHisto("btagEff","num_"+tag+"_tagged_tight",selCentralJets.at(i).Pt(),selCentralJets.at(i).Eta(),weight);
    }
    return true;
  }
  
  double apply_sf(std::vector<TLorentzVectorWithIndex> selCentralJets, std::vector<double> btags, std::vector<float> *JetAk04HadFlav, const utils::tables &tab, BTagCalibrationReader _btag_calibration_reader)
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
      double sf = _btag_calibration_reader.eval_auto_bounds("central", flavor, std::abs(selCentralJets[i].Eta()), selCentralJets[i].Pt());
      btagWeight *= tagged ? sf : (1 - sf * eff) / (1 - eff);
      //std::cout << "Jet nb " << i << " : tag = " << tag << " ; (pT,eta) = (" << selCentralJets[i].Pt() << "," << selCentralJets[i].Eta() << ") ; btag = " << btags.at(i) << " ; eff = " << eff << " ; sf = " << sf << std::endl;
    }
    //std::cout << "btagWeight: " << btagWeight << std::endl;
    return btagWeight;
  }

}
