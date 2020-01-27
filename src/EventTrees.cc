#include <EventTrees.h>


EventTrees::EventTrees(Options const &options, Dataset &dataset,
                       std::string const treeName)
    : AnalysisCommon{options, dataset},
      storeWeightSyst_{options.GetAs<std::string>("syst") == "weights"},
      outputFile_{options.GetAs<std::string>("output").c_str(), "recreate"} {
  tree_ = new TTree(treeName.c_str(), "");
  tree_->SetDirectory(&outputFile_);

  if (isSim_) {
    tree_->Branch("weight", &weight_);

    if (storeWeightSyst_) {
      int const numVariations = weightCollector_.NumVariations();
      systWeights_.resize(numVariations);
      for (int i = 0; i < numVariations; ++i) {
        auto const name = "weight_"
            + std::string{weightCollector_.VariationName(i)};
        tree_->Branch(name.c_str(), &systWeights_[i]);
      }
    }
  }
}


void EventTrees::PostProcessing() {
  outputFile_.Write();
  outputFile_.Close();
}


void EventTrees::FillTree() {
  if (isSim_) {
    if (not storeWeightSyst_)
      weight_ = weightCollector_() * intLumi_;
    else {
      weight_ = weightCollector_.NominalWeight() * intLumi_;
      for (int i = 0; i < int(systWeights_.size()); ++i)
        systWeights_[i] = weightCollector_.RelWeight(i) * weight_;
    }
  }
  tree_->Fill();
}

