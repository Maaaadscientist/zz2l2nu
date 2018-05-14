#include "./EwkCorrections.h"

namespace EwkCorrections
{
  //Reads correction table
  std::vector<std::vector<float>> readFile_and_loadEwkTable(TString catalogInputFile){
    std::ifstream myReadFile;
    std::vector<float> Table_line;
    std::vector<std::vector<float>> Table_EWK;
    TString name;
    TString cmssw_path;
    cmssw_path = getenv("CMSSW_BASE");
    TString path = cmssw_path+"/src/shears/HZZ2l2nu/data/";

    if(catalogInputFile.Contains("ZZ")) name = path+"corrections/ZZ_EwkCorrections.dat";
    if(catalogInputFile.Contains("WZ")) name = path+"corrections/WZ_EwkCorrections.dat";
    myReadFile.open(name);
    if(!myReadFile.is_open()) std::cout<<"WARNING: "+name+" NOT FOUND"<<std::endl;
    int Start=0;
    while (!myReadFile.eof()){
      Start++;
      std::string output;
      myReadFile >> output;
      if(Start%5!=0) Table_line.push_back(atof(output.c_str()));
      if(Start%5==0){
        Table_line.push_back(atof(output.c_str()));
        Table_EWK.push_back(Table_line);
        Table_line.clear();
      }
    }
    myReadFile.close();
    return Table_EWK;
  }

  //Finds the right correction in the file
  std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat){
    //find the range of sqrt s hat (each 200 lines it changes)
    unsigned int j = 0;
    float best = 0.8E+04; //highest value of sqrt s hat in the table
    if( sqrt_s_hat > best) j = 39800; //in the very rare case where we have bigger s than our table (table is for 8TeV and we run at 13TeV)
    else{
      for(unsigned int i = 0 ; i < 40000 ; i = i+200){
        if(fabs(sqrt_s_hat - Table_EWK[i][0]) < best){
          best = fabs(sqrt_s_hat - Table_EWK[i][0]);
          j = i;
        }
        else break ;
      }
    }
    best = Table_EWK[j+199][1];
    if(t_hat > best) j = j+199; //in the very rare case where we have bigger t than our table
    else{
      best = 0.1E+09;
      for(unsigned int k = j ; k < j + 200 ; k++){
        if(fabs(t_hat - Table_EWK[k][1]) < best){
          best = fabs(t_hat - Table_EWK[k][1]);
          j = k;
        }
        else break ;
      }
    }
    std::vector<float> EWK_w2_vec;
    EWK_w2_vec.push_back(Table_EWK[j][2]); //ewk corrections for quark u/c
    EWK_w2_vec.push_back(Table_EWK[j][3]); //ewk corrections for quark d/s
    EWK_w2_vec.push_back(Table_EWK[j][4]); //ewk corrections for quark b
    return EWK_w2_vec ;
  }

  std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> reconstructGenLevelBosons(std::vector<float> *GLepBarePt, std::vector<float> *GLepBareEta, std::vector<float> *GLepBarePhi, std::vector<float> *GLepBareE, std::vector<int> *GLepBareId, std::vector<int> *GLepBareSt, std::vector<int> *GLepBareMomId){
    std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons; //Convention: For Z, first is lepton and second is antilepton. For W, first is charged lepton and second is neutrino. Warning: does not work for ZZ->4l or for WW->2l2nu.
    std::cout << "Hello" << std::endl;
    for(int i = 0 ; i < GLepBarePt->size() ; i++){
      std::cout << "BareLepton with ID = " << GLepBareId->at(i) << " and status = " << GLepBareSt->at(0) << " and MomId = " << GLepBareMomId->at(i) <<  " and pT = " << GLepBarePt->at(i) << std::endl;
      if(fabs(GLepBareMomId->at(i)) == 23 && (GLepBareId->at(i) == 11 || GLepBareId->at(i) == 13 || GLepBareId->at(i) == 15)) genLevelLeptons["leptonsFromZ"].first.SetPtEtaPhiE(GLepBarePt->at(i),GLepBareEta->at(i),GLepBarePhi->at(i),GLepBareE->at(i));
      if(fabs(GLepBareMomId->at(i)) == 23 && (GLepBareId->at(i) == -11 || GLepBareId->at(i) == -13 || GLepBareId->at(i) == -15)) genLevelLeptons["leptonsFromZ"].second.SetPtEtaPhiE(GLepBarePt->at(i),GLepBareEta->at(i),GLepBarePhi->at(i),GLepBareE->at(i));
      if(fabs(GLepBareMomId->at(i)) == 23 && (GLepBareId->at(i) == 12 || GLepBareId->at(i) == 14 || GLepBareId->at(i) == 16)) genLevelLeptons["neutrinosFromZ"].first.SetPtEtaPhiE(GLepBarePt->at(i),GLepBareEta->at(i),GLepBarePhi->at(i),GLepBareE->at(i));
      if(fabs(GLepBareMomId->at(i)) == 23 && (GLepBareId->at(i) == -12 || GLepBareId->at(i) == -14 || GLepBareId->at(i) == -16)) genLevelLeptons["neutrinosFromZ"].second.SetPtEtaPhiE(GLepBarePt->at(i),GLepBareEta->at(i),GLepBarePhi->at(i),GLepBareE->at(i));
      if(fabs(GLepBareMomId->at(i)) == 24 && (fabs(GLepBareId->at(i)) == 11 || fabs(GLepBareId->at(i)) == 13 || fabs(GLepBareId->at(i)) == 15)) genLevelLeptons["leptonsFromW"].first.SetPtEtaPhiE(GLepBarePt->at(i),GLepBareEta->at(i),GLepBarePhi->at(i),GLepBareE->at(i)); //For now, no distinction between W+ and W-. It will have to come later.
      if(fabs(GLepBareMomId->at(i)) == 24 && (fabs(GLepBareId->at(i)) == 12 || fabs(GLepBareId->at(i)) == 14 || fabs(GLepBareId->at(i)) == 16)) genLevelLeptons["leptonsFromW"].second.SetPtEtaPhiE(GLepBarePt->at(i),GLepBareEta->at(i),GLepBarePhi->at(i),GLepBareE->at(i)); //For now, no distinction between W+ and W-. It will have to come later. //Don't take tau neutrinos, it means that it comes with a tau and we don't want that.
    }
    return genLevelLeptons;
  }

  //The main function, returns the kfactor
  double getEwkCorrections(TString catalogInputFile, std::map<std::string,std::pair<TLorentzVector,TLorentzVector>> genLevelLeptons, const std::vector<std::vector<float>> & Table, double & ewkCorrections_error, std::vector<float> *GPdfx1, std::vector<float> *GPdfx2, std::vector<int> *GPdfId1, std::vector<int> *GPdfId2){
    double kFactor = 1.;
    enum {ZZ, WZ};
    int event_type = -1;
    if(catalogInputFile.Contains("ZZ")) event_type = ZZ;
    else if (catalogInputFile.Contains("WZ")) event_type = WZ;
    else return 1.;
    if(event_type==ZZ && (genLevelLeptons.find("leptonsFromZ")==genLevelLeptons.end() || genLevelLeptons.find("neutrinosFromZ")==genLevelLeptons.end() )) return 1.;
    if(event_type==WZ && (genLevelLeptons.find("leptonsFromZ")==genLevelLeptons.end() || genLevelLeptons.find("leptonsFromW")==genLevelLeptons.end() )) return 1.;
    if(event_type==ZZ) std::cout << "Event is of type ZZ." << std::endl;
    if(event_type==WZ) std::cout << "Event is of type WZ." << std::endl;

    TLorentzVector l1, l2, l3, l4, V1, V2, VV;
    l1 = genLevelLeptons.at("leptonsFromZ").first;
    l2 = genLevelLeptons.at("leptonsFromZ").second;
    if(event_type==ZZ){
      l3 = genLevelLeptons.at("neutrinosFromZ").first;
      l4 = genLevelLeptons.at("neutrinosFromZ").second;
    }
    else if(event_type==WZ){
      l3 = genLevelLeptons.at("leptonsFromW").first;
      l4 = genLevelLeptons.at("leptonsFromW").second;
    }
    V1 = l1+l2;
    V2 = l3+l4;
    VV = V1+V2;
    // From there, same as in the old framework.

    double s_hat = pow(VV.M(),2);
    std::cout << "Just to check: size of x1 is " << GPdfx1->size() << " and x1 = " << GPdfx1->at(0) << std::endl;
    std::cout << "Just to check: size of x2 is " << GPdfx2->size() << " and x2 = " << GPdfx2->at(0) << std::endl;
    std::cout << "Just to check: size of Id1 is " << GPdfId1->size() << " and Id1 = " << GPdfId1->at(0) << std::endl;
    std::cout << "Just to check: size of Id2 is " << GPdfId2->size() << " and Id2 = " << GPdfId2->at(0) << std::endl;

    TLorentzVector V1_b = V1;
    TLorentzVector p1_b, p2_b;
    double energy = 6500. ; //13 TeV in total
    double x1 = GPdfx1->at(0);
    double x2 = GPdfx2->at(0);
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
    if(event_type==WZ) {
      double b = 1./2./sqrt(s_hat) * sqrt(pow(s_hat-m_z*m_z-m_w*m_w,2) - 4*m_w*m_w*m_z*m_z);
      double a = sqrt(b*b + m_z*m_z);
      t_hat = m_z*m_z - sqrt(s_hat) * (a - b * cos_theta); //awful calculation, needed to put ourselves to the center-of-mass frame with the 2 particles having a different mass !
    }
    std::cout << "Computing corrections. The value of sqrt(s_hat) is " << sqrt(s_hat) << " and t_hat is " << t_hat << std::endl;

    int quark_type = 0; //Flavour of incident quark
    if(fabs(GPdfId1->at(0)) != 21){
      if((event_type == ZZ) && (fabs((GPdfId2->at(0)) != 21) && (fabs(GPdfId2->at(0)) != fabs(GPdfId1->at(0))))) {std::cout << "Different flavours!" << std::endl; return 1.;} //No correction applied if 2 different flavours
      else quark_type = fabs(GPdfId1->at(0));
    }
    else{
      if(fabs(GPdfId2->at(0)) == 21) {std::cout << "gg case, impossible to compute corrections!" << std::endl; return 1.;} //No correction can be applied in the gg->ZZ case
      else quark_type = fabs(GPdfId2->at(0));
    }
    std::vector<float> Correction_vec = findCorrection( Table, sqrt(s_hat), t_hat ); //Extract the corrections for the values of s and t computed
    std::cout << "Correction_vec = (" << Correction_vec[0] << "," << Correction_vec[1] << "," << Correction_vec[2] << ")" << std::endl;
    std::cout << "quark_type = " << quark_type << std::endl;

    if(quark_type==1) kFactor = 1. + Correction_vec[1]; //d
    if(quark_type==2) kFactor = 1. + Correction_vec[0]; //u
    if(quark_type==3) kFactor = 1. + Correction_vec[1]; //s as d
    if(quark_type==4) kFactor = 1. + Correction_vec[0]; //c as u
    if(quark_type==5) kFactor = 1. + Correction_vec[2]; //b  //Notice that the quark types are irrelevant for the case of WZ (same numbers in the last 3 columns).

    if(sqrt(s_hat)< 2*m_z && event_type == ZZ) {std::cout << "Event is off-shell!" << std::endl; kFactor = 1.;} //Off-shell cases, not corrected to avoid non-defined values for t.
    if(sqrt(s_hat)< m_z + m_w && event_type == WZ) {std::cout << "Event is off-shell!" << std::endl; kFactor = 1.;}

    //All part about uncertainties comes here

    //All part about WZ gamma-induced contribution here

    std::cout << "kFactor = " << kFactor << std::endl;
    return kFactor;
  }

}
