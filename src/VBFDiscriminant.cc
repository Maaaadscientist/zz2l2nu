#include <FileInPath.h>
#include <HZZException.h>
#include <Logger.h>
#include <MelaHandler.h>
#include <VBFDiscriminant.h>

#include <filesystem>


namespace fs = std::filesystem;
namespace PDG = PDGHelpers;


VBFDiscriminant::VBFDiscriminant(Options const &options)
    : melaHandle_{MelaHandler::smartMela_.Get()} {

  auto mepConfig = options.GetConfig()["vbf_discriminant"];
  if (not mepConfig)
    throw HZZException{" vbf discriminant configuration is empty."};

  // Loaing Matrix Element probablities flags from configuration file
  auto const &mepFlag = mepConfig["mep_flags"];
  mepFlag_[MEP::kSIGa1] = mepFlag["p_sig_a1"].as<std::string>();
  mepFlag_[MEP::kSIGa2] = mepFlag["p_sig_a2"].as<std::string>();
  mepFlag_[MEP::kSIGa3] = mepFlag["p_sig_a3"].as<std::string>();
  mepFlag_[MEP::kSIGl1] = mepFlag["p_sig_l1"].as<std::string>();
  mepFlag_[MEP::kALT] = mepFlag["p_alt"].as<std::string>();

  // Preparing MELA hypotheses, clusters and computers
  auto const &melaOptionsList = mepConfig["meps"];
  std::string optStr;
  for (auto const &opt : melaOptionsList) {
    optStr = opt.as<std::string>();
    auto hypo = std::make_unique<MELAHypothesis>(melaHandle_, optStr);
    auto computer = std::make_unique<MELAComputation>(hypo.get());
    GMECHelperFunctions::addToMELACluster(computer.get(), clusters_);
    MELAOptionParser melaOption(optStr);
    if (melaOption.isAliased())
      aliasedHypos_.emplace_back(hypo.get());

    hypotheses_.emplace_back(std::move(hypo));
    computers_.emplace_back(std::move(computer));
  }

  for (auto &computer : computers_)
    computer->addContingencies(aliasedHypos_);

  if (not clusters_.empty()) {
    LOG_DEBUG << "ME clusters:" << std::endl;
    for (auto const cluster : clusters_) {
      LOG_DEBUG << "\t- Cluster " << cluster->getName() <<
        " has " << cluster->getComputations()->size() <<
        " computations registered." << std::endl;
    }
  }

  // Loading g and c constants
  TFile cConstantFile(FileInPath::Resolve(
        "VBFDiscriminant/SmoothKDConstant_m4l_DjjVBF_13TeV.root").c_str());
  TFile ga2ConstantFile(FileInPath::Resolve(
        "VBFDiscriminant/gConstant_VBF_g2.root").c_str());
  TFile ga3ConstantFile(FileInPath::Resolve(
        "VBFDiscriminant/gConstant_VBF_g4.root").c_str());
  TFile gl1ConstantFile(FileInPath::Resolve(
        "VBFDiscriminant/gConstant_VBF_L1.root").c_str());
  cConstant_.reset(dynamic_cast<TSpline3 *>(
              cConstantFile.Get("sp_gr_varReco_Constant_Smooth")));
  gConstant_[MEP::kSIGa2].reset(dynamic_cast<TSpline3 *>(
              ga2ConstantFile.Get("sp_tgfinal_VBF_SM_over_tgfinal_VBF_g2")));
  gConstant_[MEP::kSIGa3].reset(dynamic_cast<TSpline3 *>(
              ga3ConstantFile.Get("sp_tgfinal_VBF_SM_over_tgfinal_VBF_g4")));
  gConstant_[MEP::kSIGl1].reset(dynamic_cast<TSpline3 *>(
              gl1ConstantFile.Get("sp_tgfinal_VBF_SM_over_tgfinal_VBF_L1")));
}


VBFDiscriminant::~VBFDiscriminant() {
  for (auto cluster : clusters_)
    delete cluster;
}


std::array<double, 4> const &VBFDiscriminant::Get(
    TLorentzVector const &p4LL,
    TLorentzVector const &p4Miss,
    std::vector<Jet> const &jets) {
  Reset();
  // Check if jets size is less than 2 then
  // return invalid values for Djj VBF discriminants i.e. -1
  if (jets.size() < 2)
    // Default values of Djj VBF discriminants are -1
    return dJJVBF_;

  BuildMelaCandidate(p4LL, p4Miss, jets);
  ComputeClusters();

  // Storing P_sigs and P_alt
  double pSiga1 = 0, pSiga2 = 0, pSiga3 = 0, pSigl1 = 0, pAlt = 0;
  for (auto &computer : computers_) {
    auto const &name = computer->getName();
    if (name == mepFlag_[MEP::kSIGa1])
      pSiga1 = computer->getVal(MELAHypothesis::METype::UseME);
    else if (name == mepFlag_[MEP::kSIGa2])
      pSiga2 = computer->getVal(MELAHypothesis::METype::UseME);
    else if (name == mepFlag_[MEP::kSIGa3])
      pSiga3 = computer->getVal(MELAHypothesis::METype::UseME);
    else if (name == mepFlag_[MEP::kSIGl1])
      pSigl1 = computer->getVal(MELAHypothesis::METype::UseME);
    else if (name == mepFlag_[MEP::kALT])
      pAlt = computer->getVal(MELAHypothesis::METype::UseME);
  }

  auto const mZZ = p4ZZApprox_.M();
  auto const c = cConstant_->Eval(mZZ);
  auto const a2C = c * std::pow(gConstant_[MEP::kSIGa2]->Eval(mZZ), -2);
  auto const a3C = c * std::pow(gConstant_[MEP::kSIGa3]->Eval(mZZ), -2);
  auto const l1C = c * std::pow(gConstant_[MEP::kSIGl1]->Eval(mZZ), -2) * 1e+8;
  // Calculating the DjjVBF discriminants
  dJJVBF_[DjjVBF::a1] = pSiga1 / (pSiga1 + c * pAlt);
  dJJVBF_[DjjVBF::a2] = pSiga2 / (pSiga2 + a2C * pAlt);
  dJJVBF_[DjjVBF::a3] = pSiga3 / (pSiga3 + a3C * pAlt);
  dJJVBF_[DjjVBF::L1] = pSigl1 / (pSigl1 + l1C * pAlt);

  return dJJVBF_;
}


void VBFDiscriminant::Reset() {
  melaHandle_->resetInputEvent();
  dJJVBF_.fill(-1);
  for (auto cluster : clusters_)
    cluster->reset();
}


void VBFDiscriminant::BuildMelaCandidate(TLorentzVector const &p4LL,
    TLorentzVector const &p4Miss, std::vector<Jet> const &jets) {
  SimpleParticleCollection_t daughters, associated;

  // Construncing approximate ZZ candidate
  p4ZZApprox_.SetPtEtaPhiM(p4Miss.Pt(), p4LL.Eta(), p4Miss.Phi(), PDG::Zmass);
  p4ZZApprox_ += p4LL;

  // Adding reconstructed ZZ candidate as a duaghter particle
  daughters.push_back(SimpleParticle_t(25, p4ZZApprox_));

  // Adding jets as associated particles
  for (auto const& jet : jets)
    associated.push_back(SimpleParticle_t(0, jet.p4));

  melaHandle_->setCandidateDecayMode(TVar::CandidateDecay_Stable);
  melaHandle_->setInputEvent(&daughters, &associated, nullptr, false);
}


void VBFDiscriminant::ComputeClusters() const {
  for (auto cluster : clusters_) {
    cluster->computeAll();
    cluster->update();
  }
}

