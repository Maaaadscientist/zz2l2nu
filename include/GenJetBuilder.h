#ifndef GENJETBUILDER_H_
#define GENJETBUILDER_H_

#include <vector>

#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <Options.h>
#include <PhysicsObjects.h>


/// Lazily builds a collection of generator-level jets
class GenJetBuilder : public CollectionBuilder<GenJet> {
 public:
  GenJetBuilder(Dataset &dataset, Options const &);

  /// Returns collection of generator-level jets
  std::vector<GenJet> const &Get() const override;

 private:
  /// Constructs generator-level jets in the current event
  void Build() const override;

  /// Collection of generator-level jets
  mutable std::vector<GenJet> jets_;
  
  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_;
};

#endif  // GENJETBUILDER_H_

