#ifndef MUONBUILDER_H_
#define MUONBUILDER_H_

#include <memory>
#include <optional>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <EventCache.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <RoccoR.h>

class TRandom;


/**
 * \brief Constructs collections of reconstructed muons
 *
 * For each event two collections of muons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one. Method \ref GetMomenta is tied
 * to the tight collection.
 *
 * Rochester corrections for muon momenta are applied.
 */
class MuonBuilder : public CollectionBuilder {
 public:
  /**
   * \brief Constructor
   *
   * \param[in] reader   Reader object that provides access to the source tree.
   * \param[in] options  Configuration options for the job.
   * \param[in] randomGenerator  Reference to a common random number generator.
   */
  MuonBuilder(TTreeReader &reader, Options const &options,
              TRandom &randomGenerator);

  /// Returns collection of loose muons
  std::vector<Muon> const &GetLoose() const;

  /// Returns collection of tight muons
  std::vector<Muon> const &GetTight() const;

 private:
  /**
   * \brief Applies Rochester correction to momentum of the muon
   *
   * \param[in,out] muon  Muon to be corrected.
   * \param[in] trackerLayers  Number of tracker layers with measurements for
   *   the given muon.
   * 
   * The <a href="https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon">
   * Rochester correction</a> is applied in place. The momentum before the
   * correction is stored in Muon::uncorrP4.
   */
  void ApplyRochesterCorrection(Muon *muon, int trackerLayers) const;

  /// Constructs muons for the current event
  void Build() const;
  
  /**
   * \brief Finds matching generator-level muon using (eta, phi) metric
   *
   * \param[in] muon  Recontructed muon for which a generator-level match needs
   *   to be found.
   * \param[in] maxDR  Maximal allowed distance in the (eta, phi) metric.
   * \return  An optional that contains the matching GenParicle. If no match is
   *   found in the specified cone, the optional is empty.
   */
  std::optional<GenParticle> FindGenMatch(Muon const &muon, double maxDR) const;

  /// Returns momentum of tight moun with given index
  TLorentzVector const &GetMomentum(size_t index) const override;

  /// Returns the number of tight mouns
  size_t GetNumMomenta() const override;

  /// Minimal pt for loose electrons, GeV
  double minPtLoose_;

  /// Minimal pt for tight electrons, GeV
  double minPtTight_;

  /// Collection of muons passing loose selection
  mutable std::vector<Muon> looseMuons_;

  /// Collection of muons passing tight selection
  mutable std::vector<Muon> tightMuons_;

  /// An object to facilitate caching
  EventCache cache_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  /// Object to compute Rochester correction to muon pt
  std::unique_ptr<RoccoR> rochesterCorrection_;

  /// Reference to common random number generator
  TRandom &randomGenerator_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcE_;
  mutable TTreeReaderArray<float> srcCharge_, srcIsolation_;
  mutable TTreeReaderArray<unsigned> srcId_, srcIdTight_;
  mutable TTreeReaderArray<int> srcTrackerLayers_;
  mutable TTreeReaderArray<int> genLeptonId_;
  mutable TTreeReaderArray<float> genLeptonPt_, genLeptonEta_, genLeptonPhi_;
};


inline TLorentzVector const &MuonBuilder::GetMomentum(size_t index) const {
  return tightMuons_.at(index).p4;
}


inline size_t MuonBuilder::GetNumMomenta() const {
  return tightMuons_.size();
}

#endif  // MUONBUILDER_H_

