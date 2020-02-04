#ifndef GENPHOTONBUILDER_H_
#define GENPHOTONBUILDER_H_

#include <vector>

#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>


/// Lazily builds a collection of generator-level photons
class GenPhotonBuilder : public CollectionBuilder<GenPhoton> {
 public:
  GenPhotonBuilder(Dataset &dataset, Options const &);

  /// Returns collection of generator-level jets
  std::vector<GenPhoton> const &Get() const override;

 private:
  /// Constructs generator-level jets in the current event
  void Build() const override;

  /// Collection of generator-level jets
  mutable std::vector<GenPhoton> photons_;
  
  mutable TTreeReaderArray<int> srcPhotonGenPartIndex_;
  mutable TTreeReaderArray<UChar_t> srcFlavour_;
  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_;
};

#endif  // GENPHOTONBUILDER_H_

