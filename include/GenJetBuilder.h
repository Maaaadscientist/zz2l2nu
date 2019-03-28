#ifndef GENJETBUILDER_H_
#define GENJETBUILDER_H_

#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>


/// Builds collection of generator-level jets
class GenJetBuilder : public CollectionBuilder {
 public:
  GenJetBuilder(TTreeReader &reader, Options const &);
  std::vector<GenJet> const &Get() const;

 private:
  /// Constructs generator-level jets in the current event
  void Build() const;

  /// Returns momentum of jet with given index
  TLorentzVector const &GetMomentum(size_t index) const override;

  /// Returns number of jets
  size_t GetNumMomenta() const override;

  /// Collection of generator-level jets
  mutable std::vector<GenJet> jets_;
  
  /// An object to facilitate caching
  EventCache cache_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_;
};


inline TLorentzVector const &GenJetBuilder::GetMomentum(size_t index) const {
  return jets_.at(index).p4;
}


inline size_t GenJetBuilder::GetNumMomenta() const {
  return jets_.size();
}

#endif  // GENJETBUILDER_H_

