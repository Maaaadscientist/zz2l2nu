#ifndef ELECTRONBUILDER_H_
#define ELECTRONBUILDER_H_

#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <EventCache.h>
#include <PhysicsObjects.h>
#include <Options.h>


/**
 * \brief Constructs collections of reconstructed electrons
 *
 * For each event two collections of electrons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one.
 */
class ElectronBuilder {
 public:
  /// Constructor
  ElectronBuilder(TTreeReader &reader, Options const &);

  /// Returns collection of loose electrons
  std::vector<Electron> const &GetLooseElectrons() const;

  /// Returns collection of tight electrons
  std::vector<Electron> const &GetTightElectrons() const;

 private:
  /// Constructs electrons for the current event
  void Update() const;

  /// Minimal pt for loose electrons, GeV
  double minPtLoose;

  /// Minimal pt for tight electrons, GeV
  double minPtTight;

  /// Collection of electrons passing loose selection
  mutable std::vector<Electron> looseElectrons;

  /// Collection of electrons passing tight selection
  mutable std::vector<Electron> tightElectrons;

  /// An object to facilitate caching
  EventCache cache_;

  TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_, srcEtaSc_;
  TTreeReaderArray<float> srcCharge_;
  TTreeReaderArray<unsigned> srcId_;
};

#endif  // ELECTRONBUILDER_H_

