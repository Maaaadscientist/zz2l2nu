#ifndef ELECTRONBUILDER_H_
#define ELECTRONBUILDER_H_

#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <PhysicsObjects.h>
#include <Options.h>


/**
 * \brief Constructs collections of reconstructed electrons
 *
 * For each event two collections of electrons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one.
 *
 * In each event operator() must be called on this object before the collections
 * can be accessed.
 */
class ElectronBuilder {
 public:
  /// Constructor
  ElectronBuilder(TTreeReader &reader, Options const &);

  /// Returns collection of loose electrons
  std::vector<Electron> const &GetLooseElectrons() const;

  /// Returns collection of tight electrons
  std::vector<Electron> const &GetTightElectrons() const;

  /// Constructs electrons for the current event
  void operator()();

 private:
  /// Minimal pt for loose electrons, GeV
  double minPtLoose;

  /// Minimal pt for tight electrons, GeV
  double minPtTight;

  /// Collection of electrons passing loose selection
  std::vector<Electron> looseElectrons;

  /// Collection of electrons passing tight selection
  std::vector<Electron> tightElectrons;

  TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_, srcEtaSc_;
  TTreeReaderArray<float> srcCharge_;
  TTreeReaderArray<unsigned> srcId_;
};


inline std::vector<Electron> const &ElectronBuilder::GetLooseElectrons() const {
  return looseElectrons;
}

inline std::vector<Electron> const &ElectronBuilder::GetTightElectrons() const {
  return tightElectrons;
}

#endif  // ELECTRONBUILDER_H_

