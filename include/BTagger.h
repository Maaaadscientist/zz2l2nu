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

  /// Checks whether a given jet is b-tagged if the jet is taggable
  bool operator()(Jet const &jet) const;
  
  /// Checks whether a given jet is taggable
  bool IsTaggable(Jet const &jet) const;

 private:
  /// b tag threshold
  double bTagCut_;

  /// Minimum jet pt which b-tagging is reliable
  double ptCut_;

  /// Maximum magnitude of jet eta in order to be within the tracker coverage
  double etaCut_;
};

#endif  // BTAGGER_H_
