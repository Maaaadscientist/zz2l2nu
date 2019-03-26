#ifndef ELECTRONBUILDER_H_
#define ELECTRONBUILDER_H_

#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <PhysicsObjects.h>
#include <Options.h>


/**
 * \brief Constructs collections of reconstructed electrons
 *
 * For each event two collections of electrons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one. Method \ref GetMomenta is tied
 * to the tight electrons.
 */
class ElectronBuilder : public CollectionBuilder {
 public:
  /// Constructor
  ElectronBuilder(TTreeReader &reader, Options const &);

  /// Returns collection of loose electrons
  std::vector<Electron> const &GetLoose() const;

  /// Returns collection of tight electrons
  std::vector<Electron> const &GetTight() const;

 private:
  /// Constructs electrons for the current event
  void Build() const;

  /// Returns momentum of tight electron with given index
  TLorentzVector const &GetMomentum(size_t index) const override;

  /// Returns the number of tight electrons
  size_t GetNumMomenta() const override;

  /// Minimal pt for loose electrons, GeV
  double minPtLoose_;

  /// Minimal pt for tight electrons, GeV
  double minPtTight_;

  /// Collection of electrons passing loose selection
  mutable std::vector<Electron> looseElectrons_;

  /// Collection of electrons passing tight selection
  mutable std::vector<Electron> tightElectrons_;

  /// An object to facilitate caching
  EventCache cache_;

  TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_, srcEtaSc_;
  TTreeReaderArray<float> srcCharge_;
  TTreeReaderArray<unsigned> srcId_;
};


inline TLorentzVector const &ElectronBuilder::GetMomentum(size_t index) const {
  return tightElectrons_.at(index).p4;
}


inline size_t ElectronBuilder::GetNumMomenta() const {
  return tightElectrons_.size();
}

#endif  // ELECTRONBUILDER_H_

