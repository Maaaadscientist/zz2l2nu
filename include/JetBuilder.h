#ifndef JETBUILDER_H_
#define JETBUILDER_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <CollectionBuilder.h>
#include <Dataset.h>
#include <GenJetBuilder.h>
#include <JetCorrector.h>
#include <Options.h>
#include <PhysicsObjects.h>


/**
 * \brief Lazily builds collection of reconstructed jets
 *
 * Reads parameters from section "jets" of the master configuration. Systematic
 * variations in jets are implemented with the help of JetCorrector, which also
 * applies JER smearing. To follow the standard smearing algorithm, this builder
 * needs to be made aware of generator-level jets via method
 * \ref SetGenJetBuilder.
 *
 * Jet with fully corrected pt > 15 GeV are aggregated for
 * \ref GetSumMomentumShift to be used for the type 1 correction of missing pt.
 */
class JetBuilder : public CollectionBuilder<Jet> {
 public:
  JetBuilder(Dataset &dataset, Options const &options,
             TabulatedRngEngine &rngEngine);

  /// Returns collection of jets
  std::vector<Jet> const &Get() const override;

  /**
   * \brief Specifies an object that provides generator-level jets
   *
   * Generator-level jets are used for JER smearing. However, if no such object
   * is given, stochastic version of the smearing will be performed. When given,
   * the object must have an approriate life time.
   */
  void SetGenJetBuilder(GenJetBuilder const *genJetBuilder);

 private:
  /// Constructs jets in the current event
  void Build() const override;

  /**
   * \brief Finds matching generator-level jet
   *
   * Returns a nullptr if no match is found within the allowed cone.
   */
  GenJet const *FindGenMatch(TLorentzVector const &p4,
                             double ptResolution) const;

  /**
   * \brief Non-owning pointer to an object that produces generator-level jets
   *
   * May be nullptr.
   */
  GenJetBuilder const *genJetBuilder_;

  /// Minimal pt for jets, GeV
  double minPt_;

  /// Maximal |eta| for jets
  double maxAbsEta_;

  /// Collection of jets
  mutable std::vector<Jet> jets_;

  /// Indicates whether running on simulation or data
  bool isSim_;

  /// Index of the bit with jet ID decision
  int jetIdBit_;

  /// Object that computes JEC
  JetCorrector jetCorrector_;

  mutable TTreeReaderArray<float> srcPt_, srcEta_, srcPhi_, srcMass_;
  mutable TTreeReaderArray<float> srcArea_, srcRawFactor_;
  mutable TTreeReaderArray<float> srcBTag_;
  mutable TTreeReaderArray<int> srcId_;
  mutable TTreeReaderValue<float> puRho_;
  mutable std::unique_ptr<TTreeReaderArray<int>> srcHadronFlavour_;

  // Properties of soft jets, which are used in the type 1 correction of ptmiss
  mutable TTreeReaderArray<float> softRawPt_, softEta_, softPhi_, softArea_;
};

#endif  // JETBUILDER_H_

