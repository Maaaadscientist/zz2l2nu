#ifndef BTAGGER_H_
#define BTAGGER_H_

#include <Options.h>
#include <PhysicsObjects.h>

/**
 * \brief Evaluates whether a jet is b-tagged or not
 *
 * Please use this class instead of relying on the raw value of the btagging
 * discriminator.
 */
class BTagger {
 public:
  /// Constructor from configuration options
  BTagger(Options const &options);

  /// Checks whether a given jet is b-tagged(default:loose) if the jet is taggable
  bool operator()(Jet const &jet) const;

  /// Checks whether a given jet is b-tagged(loose) if the jet is taggable
  bool Loose(Jet const &jet) const;

  /// Checks whether a given jet is b-tagged(medium) if the jet is taggable
  bool Medium(Jet const &jet) const;

  /// Checks whether a given jet is b-tagged(tight) if the jet is taggable
  bool Tight(Jet const &jet) const;

  /// Checks whether a given jet is taggable
  bool IsTaggable(Jet const &jet) const;

 private:
  /// b tag loose threshold
  double bTagCutLoose_;

  /// b tag medium threshold
  double bTagCutMedium_;

  /// b tag tight threshold
  double bTagCutTight_;

  /// Minimum jet pt which b-tagging is reliable
  double ptCut_;

  /// Maximum magnitude of jet eta in order to be within the tracker coverage
  double etaCut_;
};

#endif  // BTAGGER_H_
