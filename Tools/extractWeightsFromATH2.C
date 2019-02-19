#include <iostream>
#include <fstream>

//////////////////////////////////
//// Config global parameters ////
//////////////////////////////////

TString base_path = std::string(getenv("CMSSW_BASE")) + "/src/shears/HZZ2l2nu/"; //The location of the script and the weight file
TString weightFile = base_path + "/WeightsAndDatadriven/egammaEffi.txt_EGM2D_tight.root"; //The name (location) of the weight file
TString nameOfHisto = "EGamma_SF2D"; //The name of the histo with the weights
TString output = "Tools/generated_cpp_code.C"; //Name of the output
TString bin_X_name = "eta"; //The name you want for your X axis
TString bin_Y_name = "pt"; //The name you want for your Y axis

//////////////////////////////////
//////////////  End  /////////////
//////////////////////////////////


void writeCondition(ofstream& outputFile, TString bin_name, double binLow, double binUp, int bin, int lastBin = -1, int numberOfTabulations =0) 
{
  TString condition;
  if(bin == 1) condition = "if";
  else condition = "}else if";
  TString tabs = "";
  for(int i = 0; i < numberOfTabulations; i++) tabs += "\t";
  
  
  if(bin != lastBin) outputFile << tabs << condition << "(" << bin_name <<" >= " << binLow << " && " << bin_name << " < " << binUp << "){";
  else outputFile << "\t}else{";
}

void writeEfficiency(ofstream& outputFile, double weight)
{
  outputFile << " eff.first=" << weight << "; eff.second=0.00; //Error is not in the file...\n";
}

void getWeight()
{
  TFile *f_weight = TFile::Open(weightFile);
  TH2D *h_weight = (TH2D*) f_weight->Get(nameOfHisto);
  double weight = 0.;
  double binLow_X = 0., binUp_X = 0., binLow_Y = 0., binUp_Y =0.;
  ofstream outputFile;
  outputFile.open(base_path+output);
  for(int bin_X = 1; bin_X <= h_weight->GetNbinsX(); bin_X++){
    binLow_X = h_weight->GetXaxis()->GetBinLowEdge(bin_X);
    binUp_X = binLow_X + h_weight->GetXaxis()->GetBinWidth(bin_X);
    writeCondition(outputFile, bin_X_name, binLow_X, binUp_X, bin_X);
    outputFile << "\n";
    for(int bin_Y = 1; bin_Y <= h_weight->GetNbinsY(); bin_Y++){
      binLow_Y = h_weight->GetYaxis()->GetBinLowEdge(bin_Y);
      binUp_Y = binLow_Y + h_weight->GetYaxis()->GetBinWidth(bin_Y);
      writeCondition(outputFile, bin_Y_name, binLow_Y, binUp_Y, bin_Y, h_weight->GetNbinsY(), 1);
      weight=h_weight->GetBinContent(bin_X, bin_Y);
//      std::cout << "weight for bin " << h_weight->FindBin(bin_X, bin_Y) << " is " << weight <<std::endl;
      writeEfficiency(outputFile, weight);
    }
    outputFile << "\t}\n";//close parenthesis
  }
  outputFile << "}";//close parenthesis
  f_weight->Close();
  outputFile.close();
}

void extractWeightsFromATH2()
{
  gROOT->SetBatch();
  getWeight();
}
