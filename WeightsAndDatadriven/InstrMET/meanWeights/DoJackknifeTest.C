{
TFile *f = TFile::Open("jackknife_validation_2017_geq2jets_withEta.root");

//std::vector<float> mT_thresholds_ = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000, 9999};
std::vector<float> mT_thresholds_ = {100, 200, 300, 350, 400, 450, 500, 550, 600, 700, 850, 1000, 1250, 1500, 3000, 9999};

std::vector<TH1F*> distributions; // index = bin

std::cout << std::endl;
std::cout << std::endl;
std::cout << "{\\scriptsize" << std::endl;
std::cout << "\\begin{center}" << std::endl;
std::cout << "\\begin{tabular}{c|c|c|c}" << std::endl;

std::cout << "$m_T$ range (GeV) & $\\lambda$ & $\\sigma_{\\text{jackknife replicas}}$ & $r$ \\\\ \\hline" << std::endl;

for (int i = 0 ; i < int(mT_thresholds_.size()-1) ; i++) {
  distributions.push_back( (TH1F*) f->Get("lambdas_jackknife_"+(TString)std::to_string(i)));
  float m_min = mT_thresholds_[i];
  float m_max = mT_thresholds_[i+1];
  float lambda = distributions[i]->GetMean();
  if (lambda == 0) continue;
  float sigma = distributions[i]->GetStdDev();
  float delta_sigma = sigma/sqrt(2.*99.);
  float r = (sqrt(100.)*sigma - sqrt(lambda))/sqrt(lambda);
  float delta_r = sqrt(100./lambda)*delta_sigma;
  std::cout << "{[}" << m_min << " --- " << m_max << "] & $" << lambda << "$ & $" << sigma << " \\pm " << delta_sigma << "$ & $" << r << " \\pm " << delta_r << "$ \\\\" << std::endl;
}

std::cout << "\\end{tabular}" << std::endl;
std::cout << "\\end{center}" << std::endl;
std::cout << "}" << std::endl;
std::cout << std::endl;
std::cout << std::endl;


}
