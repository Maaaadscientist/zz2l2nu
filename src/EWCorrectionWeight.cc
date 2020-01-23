#include <EWCorrectionWeight.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <Logger.h>
#include <FileInPath.h>


namespace fs = std::filesystem;
using namespace std;


EWCorrectionWeight::EWCorrectionWeight(Dataset &dataset, Options const &options)
    : cache_{dataset.Reader()},
      genPartPt_{dataset.Reader(), "GenPart_pt"},
      genPartEta_{dataset.Reader(), "GenPart_eta"},
      genPartPhi_{dataset.Reader(), "GenPart_phi"},
      genPartMass_{dataset.Reader(), "GenPart_mass"},
      genPartPdgId_{dataset.Reader(), "GenPart_pdgId"},
      genPartIdxMother_{dataset.Reader(), "GenPart_genPartIdxMother"},
      generatorX1_{dataset.Reader(), "Generator_x1"},
      generatorX2_{dataset.Reader(), "Generator_x2"},
      generatorId1_{dataset.Reader(), "Generator_id1"},
      generatorId2_{dataset.Reader(), "Generator_id2"} {

  auto const settingsNode = dataset.Info().Parameters()["ew_correction"];
  std::string typeLabel;

  if (settingsNode and not settingsNode.IsNull()) {
    typeLabel = settingsNode.as<std::string>();

    if (typeLabel == "ZZ")
      correctionType_ = Type::ZZ;
    else if (typeLabel == "WZ")
      correctionType_ = Type::WZ;
    else {
      std::ostringstream message;
      message << "Unknown type \"" << typeLabel << "\" for EW correction.";
      throw std::runtime_error(message.str());
    }
  } else
    correctionType_ = Type::None;

  if (correctionType_ != Type::None) {
    LOG_DEBUG << "Will apply EW corrections of type \"" << typeLabel << "\".";
    readFile_and_loadEwkTable();
  } else
    LOG_DEBUG << "Will not apply EW corrections.";

  auto const systLabel = options.GetAs<std::string>("syst");
  if (systLabel == "ewk_up")
    systDirection_ = +1;
  else if (systLabel == "ewk_down")
    systDirection_ = -1;
  else
    systDirection_ = 0;
}


double EWCorrectionWeight::NominalWeight() const {
  if (correctionType_ == Type::None)
    return 1.;
  if (cache_.IsUpdated())
    Update();
  return weightNominal_;
}

int EWCorrectionWeight::NumVariations() const {
  if (correctionType_ == Type::None)
    return 0;
  else
    return 2;
}


double EWCorrectionWeight::operator()() const {
  if (correctionType_ == Type::None)
    return 1.;
  if (cache_.IsUpdated())
    Update();
  return weightNominal_ + systDirection_ * weightError_;
}


double EWCorrectionWeight::RelWeight(int variation) const {
  if (cache_.IsUpdated())
    Update();
  if (variation == 0)
    return 1. + weightError_ / weightNominal_;
  else
    return 1. - weightError_ / weightNominal_;
}


std::string_view EWCorrectionWeight::VariationName(int variation) const {
  switch (variation) {
    case 0:
      return "ewk_up";
    case 1:
      return "ewk_down";
    default:
      return "";
  }
}

  
void EWCorrectionWeight::readFile_and_loadEwkTable(){
  std::vector<float> Table_line;
  ewTable_.clear();

  std::string name;
  
  if (correctionType_ == Type::ZZ)
    name = "ZZ_EwkCorrections.dat";
  else if (correctionType_ == Type::WZ)
    name = "WZ_EwkCorrections.dat";

  std::ifstream myReadFile{FileInPath::Resolve("corrections", name)};

  if (not myReadFile.is_open()) {
    std::ostringstream message;
    message << "File \"" << name << "\" with EW corrections is not found.";
    throw std::runtime_error(message.str());
  }

  int Start=0;
  while (!myReadFile.eof()){
    Start++;
    std::string output;
    myReadFile >> output;
    if(Start%5!=0) Table_line.push_back(atof(output.c_str()));
    if(Start%5==0){
      Table_line.push_back(atof(output.c_str()));
      ewTable_.push_back(Table_line);
      Table_line.clear();
    }
  }
  myReadFile.close();
}


std::vector<float> EWCorrectionWeight::findCorrection(float sqrt_s_hat, float t_hat) const {
  //find the range of sqrt s hat (each 200 lines it changes)
  unsigned int j = 0;
  float best = 0.8E+04; //highest value of sqrt s hat in the table
  if( sqrt_s_hat > best) j = 39800; //in the very rare case where we have bigger s than our table (table is for 8TeV and we run at 13TeV)
  else{
    for(unsigned int i = 0 ; i < 40000 ; i = i+200){
      if(fabs(sqrt_s_hat - ewTable_[i][0]) < best){
        best = fabs(sqrt_s_hat - ewTable_[i][0]);
        j = i;
      }
      else break ;
    }
  }
  best = ewTable_[j+199][1];
  if(t_hat > best) j = j+199; //in the very rare case where we have bigger t than our table
  else{
    best = 0.1E+09;
    for(unsigned int k = j ; k < j + 200 ; k++){
      if(fabs(t_hat - ewTable_[k][1]) < best){
        best = fabs(t_hat - ewTable_[k][1]);
        j = k;
      }
      else break ;
    }
  }
  std::vector<float> EWK_w2_vec;
  EWK_w2_vec.push_back(ewTable_[j][2]); //ewk corrections for quark u/c
  EWK_w2_vec.push_back(ewTable_[j][3]); //ewk corrections for quark d/s
  EWK_w2_vec.push_back(ewTable_[j][4]); //ewk corrections for quark b
  return EWK_w2_vec ;
}


std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> EWCorrectionWeight::reconstructGenLevelBosons() const {
  std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons; //Convention: For Z, first is lepton and second is antilepton. For W, first is charged lepton and second is neutrino. Warning: does not work for ZZ->4l or for WW->2l2nu.
  //std::cout << "====================================================================================================================================" << std::endl;
  //std::cout << "New event." << std::endl;
  for (int i = 0; i < int(genPartPt_.GetSize()); i++) {
    //std::cout << "BareLepton with ID = " << genPartPdgId_[i] << " and status = " << GenPart_status[0] << " and MomId = " << genPartPdgId_[genPartIdxMother_[i]] <<  " and pT = " << genPartPt_[i] << std::endl;
    if(fabs(genPartPdgId_[genPartIdxMother_[i]]) == 23 && (genPartPdgId_[i] == 11 || genPartPdgId_[i] == 13 || genPartPdgId_[i] == 15)) genLevelLeptons["leptonsFromZ"].first.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(fabs(genPartPdgId_[genPartIdxMother_[i]]) == 23 && (genPartPdgId_[i] == -11 || genPartPdgId_[i] == -13 || genPartPdgId_[i] == -15)) genLevelLeptons["leptonsFromZ"].second.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(fabs(genPartPdgId_[genPartIdxMother_[i]]) == 23 && (genPartPdgId_[i] == 12 || genPartPdgId_[i] == 14 || genPartPdgId_[i] == 16)) genLevelLeptons["neutrinosFromZ"].first.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(fabs(genPartPdgId_[genPartIdxMother_[i]]) == 23 && (genPartPdgId_[i] == -12 || genPartPdgId_[i] == -14 || genPartPdgId_[i] == -16)) genLevelLeptons["neutrinosFromZ"].second.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(genPartPdgId_[genPartIdxMother_[i]] == 24 && (fabs(genPartPdgId_[i]) == 11 || fabs(genPartPdgId_[i]) == 13 || fabs(genPartPdgId_[i]) == 15)) genLevelLeptons["leptonsFromWp"].first.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(genPartPdgId_[genPartIdxMother_[i]] == 24 && (fabs(genPartPdgId_[i]) == 12 || fabs(genPartPdgId_[i]) == 14 || fabs(genPartPdgId_[i]) == 16)) genLevelLeptons["leptonsFromWp"].second.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(genPartPdgId_[genPartIdxMother_[i]] == -24 && (fabs(genPartPdgId_[i]) == 11 || fabs(genPartPdgId_[i]) == 13 || fabs(genPartPdgId_[i]) == 15)) genLevelLeptons["leptonsFromWm"].first.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
    if(genPartPdgId_[genPartIdxMother_[i]] == -24 && (fabs(genPartPdgId_[i]) == 12 || fabs(genPartPdgId_[i]) == 14 || fabs(genPartPdgId_[i]) == 16)) genLevelLeptons["leptonsFromWm"].second.SetPtEtaPhiM(genPartPt_[i],genPartEta_[i],genPartPhi_[i],genPartMass_[i]);
  }
  return genLevelLeptons;
}


double EWCorrectionWeight::getEwkCorrections(std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons, double & error) const {
  double kFactor = 1.;
  enum {ZZ, WZp, WZm};
  int event_type = -1;

  if (correctionType_ == Type::ZZ)
    event_type = ZZ;
  else if (correctionType_ == Type::WZ) {
    event_type = WZp;
  }
  else
    return 1.;

  if(event_type==ZZ && (genLevelLeptons.find("leptonsFromZ")==genLevelLeptons.end() || genLevelLeptons.find("neutrinosFromZ")==genLevelLeptons.end() )) return 1.;
  if(event_type==WZp && genLevelLeptons.find("leptonsFromZ")==genLevelLeptons.end()) return 1.;
  if(event_type==WZp && genLevelLeptons.find("leptonsFromWp")==genLevelLeptons.end()) event_type = WZm;
  if(event_type==WZm && genLevelLeptons.find("leptonsFromWm")==genLevelLeptons.end()) return 1.;
  //if(event_type==ZZ) std::cout << "Event is of type ZZ." << std::endl;
  //if(event_type==WZp) std::cout << "Event is of type W+Z." << std::endl;
  //if(event_type==WZm) std::cout << "Event is of type W-Z." << std::endl;

  TLorentzVector l1, l2, l3, l4, V1, V2, VV;
  l1 = genLevelLeptons.at("leptonsFromZ").first;
  l2 = genLevelLeptons.at("leptonsFromZ").second;
  if(event_type==ZZ){
    l3 = genLevelLeptons.at("neutrinosFromZ").first;
    l4 = genLevelLeptons.at("neutrinosFromZ").second;
  }
  else if(event_type==WZp){
    l3 = genLevelLeptons.at("leptonsFromWp").first;
    l4 = genLevelLeptons.at("leptonsFromWp").second;
  }
  else if(event_type==WZm){
    l3 = genLevelLeptons.at("leptonsFromWm").first;
    l4 = genLevelLeptons.at("leptonsFromWm").second;
  }
  V1 = l1+l2;
  V2 = l3+l4;
  VV = V1+V2;
  // From there, same as in the old framework.

  double s_hat = pow(VV.M(),2);

  TLorentzVector V1_b = V1;
  TLorentzVector p1_b, p2_b;
  double energy = 6500. ; //13 TeV in total
  double x1 = *generatorX1_;
  double x2 = *generatorX2_;
  p1_b.SetXYZT(0.,0.,x1*energy,x1*energy); //x1 = fraction of momentum taken by the particle initiating the hard process
  p2_b.SetXYZT(0.,0.,-x2*energy,x2*energy);
  V1_b.Boost( -VV.BoostVector()); //Inverse Lorentz transformation, to get to the center-of-mass frame
  p1_b.Boost( -VV.BoostVector());
  p2_b.Boost( -VV.BoostVector());

  //Unitary vectors
  TLorentzVector V1_b_u = V1_b*(1/V1_b.P()); //Normalized to 1
  TLorentzVector p1_b_u = p1_b*(1/p1_b.P());
  TLorentzVector p2_b_u = p2_b*(1/p2_b.P());

  //Effective beam axis
  TLorentzVector diff_p = p1_b_u - p2_b_u;
  TLorentzVector eff_beam_axis = diff_p*(1./diff_p.P());
  double cos_theta = eff_beam_axis.X()*V1_b_u.X() + eff_beam_axis.Y()*V1_b_u.Y() + eff_beam_axis.Z()*V1_b_u.Z();

  double m_z = 91.1876; //Z bosons assumed to be on-shell
  double m_w = 80.385;
  double t_hat = 0.;

  if(event_type==ZZ) t_hat = m_z*m_z - 0.5*s_hat + cos_theta * sqrt( 0.25*s_hat*s_hat - m_z*m_z*s_hat );
  if((event_type==WZp) || (event_type==WZm)) {
    double b = 1./2./sqrt(s_hat) * sqrt(pow(s_hat-m_z*m_z-m_w*m_w,2) - 4*m_w*m_w*m_z*m_z);
    double a = sqrt(b*b + m_z*m_z);
    t_hat = m_z*m_z - sqrt(s_hat) * (a - b * cos_theta); //awful calculation, needed to put ourselves to the center-of-mass frame with the 2 particles having a different mass !
  }
  //std::cout << "Computing corrections. The value of sqrt(s_hat) is " << sqrt(s_hat) << " and t_hat is " << t_hat << std::endl;

  int quark_type = 0; //Flavour of incident quark
  if(fabs(*generatorId1_) != 21){
    if((event_type == ZZ) && (fabs((*generatorId2_) != 21) && (fabs(*generatorId2_) != fabs(*generatorId1_)))) {/*std::cout << "Different flavours!" << std::endl;*/ return 1.;} //No correction applied if 2 different flavours
    else quark_type = fabs(*generatorId1_);
  }
  else{
    if(fabs(*generatorId2_) == 21) {/*std::cout << "gg case, impossible to compute corrections!" << std::endl;*/ return 1.;} //No correction can be applied in the gg->ZZ case
    else quark_type = fabs(*generatorId2_);
  }
  std::vector<float> Correction_vec = findCorrection(sqrt(s_hat), t_hat ); //Extract the corrections for the values of s and t computed
  //std::cout << "Correction_vec = (" << Correction_vec[0] << "," << Correction_vec[1] << "," << Correction_vec[2] << ")" << std::endl;
  //std::cout << "quark_type = " << quark_type << std::endl;

  if(quark_type==1) kFactor = 1. + Correction_vec[1]; //d
  if(quark_type==2) kFactor = 1. + Correction_vec[0]; //u
  if(quark_type==3) kFactor = 1. + Correction_vec[1]; //s as d
  if(quark_type==4) kFactor = 1. + Correction_vec[0]; //c as u
  if(quark_type==5) kFactor = 1. + Correction_vec[2]; //b  //Notice that the quark types are irrelevant for the case of WZ (same numbers in the last 3 columns).

  if(sqrt(s_hat)< 2*m_z && event_type == ZZ) {/*std::cout << "Event is off-shell!" << std::endl;*/ kFactor = 1.;} //Off-shell cases, not corrected to avoid non-defined values for t.
  if(sqrt(s_hat)< m_z + m_w && (event_type == WZp || event_type == WZm)) {/*std::cout << "Event is off-shell!" << std::endl;*/ kFactor = 1.;}

  //std::cout << "For only virtual corrections, kFactor = " << kFactor << std::endl;

  //Uncertainty
  double kFactor_QCD = 1.;
  if(event_type == ZZ) kFactor_QCD = 15.99/9.89; //From arXiv1105.0020 //FIXME Check if this number is still up-to-date
  if(event_type == WZp) kFactor_QCD = 28.55/15.51; //for W+Z
  if(event_type == WZm) kFactor_QCD = 18.19/9.53; //for W-Z
  //Rho variable
  double rho = (l1+l2+l3+l4).Pt() / (l1.Pt() + l2.Pt() + l3.Pt() + l4.Pt());
  //std::cout << "rho = " << rho << std::endl;
  if(rho<0.3) {error = fabs((kFactor-1)*(kFactor_QCD -1));/*std::cout << "rho is small, uncertainty only in EWK X QCD." << std::endl;*/} //If rho is small: only corrections in QCD X EWK
  else {error = fabs(1-kFactor);/*std::cout << "rho is big, uncertainty is 1 on the k-factor." << std::endl;*/} //If rho is large: 100% because of large colinear gluon radiations

  //At this point, we have the relative error on the delta_ewk ( = k_ewk -1 )
  //Let's - instead - return the absolute error on k: we do delta_ewk* the_relative_errir_on_it. This gives absolute error on delta, and so on k
  error = fabs(error*kFactor);

  //std::cout << "Uncertainty = " << error << std::endl;

  //WZ: gamma-induced contribution
  if(event_type == WZp){
    kFactor *= (1 + 0.00559445 - 5.17082e-6 * sqrt(s_hat) + 3.63331e-8 * s_hat); //FIXME Check the function with latest version of LUXqed.
    double gamma_induced_uncertainty = 0.;//0.00286804 - 8.4624e-6 * sqrt(s_hat) + 3.90611e-8 * s_hat; //FIXME These were neglected in the previous incarnation of the code. But when I compute these uncertainties, they are not particularly negligible. This needs to be checked (I don't know if the formula above is strictly correct).
    //std::cout << "Uncertainty on the gamma-induced component = " << gamma_induced_uncertainty << std::endl;
    error = sqrt(pow(error,2) + pow(gamma_induced_uncertainty,2));
  }
  if(event_type == WZm){
    kFactor *= (1 + 0.00174737 + 1.70668e-5 * sqrt(s_hat) + 2.26398e-8 * s_hat);
    double gamma_induced_uncertainty = 0.;//0.00417376 - 1.51319e-5 * sqrt(s_hat) + 5.68576e-8 * s_hat;
    //std::cout << "Uncertainty on the gamma-induced component = " << gamma_induced_uncertainty << std::endl;
    error = sqrt(pow(error,2) + pow(gamma_induced_uncertainty,2));
  }
  //std::cout << "Final uncertainty on electroweak corrections = " << error << std::endl;

  //std::cout << "Total kFactor = " << kFactor << std::endl;
  return kFactor;
}


void EWCorrectionWeight::Update() const {
  if (correctionType_ == Type::None)
    return;

  auto const genLevelLeptons = reconstructGenLevelBosons();
  weightNominal_ = getEwkCorrections(genLevelLeptons, weightError_);
}

