#ifndef MUONBUILDER_H_
#define MUONBUILDER_H_

#include <memory>
#include <optional>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <CollectionBuilder.h>
#include <Options.h>
#include <PhysicsObjects.h>
#include <RoccoR.h>

class TRandom;


/**
 * \brief Lazily builds collections of reconstructed muons
 *
 * For each event two collections of muons are constructed: tight and loose.
 * They differ in the minimal pt cut as well as identification requirements. The
 * tight collection is a subset of the loose one.
 *
 * Rochester corrections for muon momenta are applied.
 */
class MuonBuilder : public CollectionBuilder<Muon> {
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

  /// Alias for \ref GetTight
  std::vector<Muon> const &Get() const override;

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
  void Build() const override;
  
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

  /// Minimal pt for loose electrons, GeV
  double minPtLoose_;

  /// Minimal pt for tight electrons, GeV
  double minPtTight_;

  /// Collection of muons passing loose selection
  mutable std::vector<Muon> looseMuons_;

  /// Collection of muons passing tight selection
  mutable std::vector<Muon> tightMuons_;

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


inline std::vector<Muon> const &MuonBuilder::Get() const {
  return GetTight();
}

#endif  // MUONBUILDER_H_

