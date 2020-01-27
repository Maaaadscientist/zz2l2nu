#ifndef HZZ2L2NU_INCLUDE_EVENTTREES_H_
#define HZZ2L2NU_INCLUDE_EVENTTREES_H_

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include <AnalysisCommon.h>
#include <Dataset.h>
#include <Options.h>


/**
 * \brief Base class for an analysis that produces trees with per-event entries
 *
 * Derived class must create branches in the output tree using method AddBranch.
 * For each event the tree is filled by calling FillTree, which also sets event
 * weights.
 *
 * Normally only the default event weight is saved. If the command line option
 * <tt>--syst=weights</tt> is provided, nominal weight as well as weights for
 * all registered weight-based systematic variations are stored. The latter ones
 * are saved as full as opposed to relative weights.
 */
class EventTrees : public AnalysisCommon {
 public:
  /**
   * \brief Contructor
   *
   * Argument \c treeName specifies the name for the tree that will be produced.
   */
  EventTrees(Options const &options, Dataset &dataset,
             std::string const treeName = "Vars");

  /// Writes the output file
  void PostProcessing();

 protected:
  /// Adds a new branch to the underlying tree
  template<typename... Args>
  void AddBranch(Args... args) {
    tree_->Branch(args...);
  }

  /**
   * \brief Fills the underlying tree
   *
   * Event weights are set automatically.
   */
  void FillTree();

 private:
  /// Indicates whether this is simulation or real data
  bool isSim_;

  /// Indicates whether variations in event weights should be stored
  bool storeWeightSyst_;

  /// Output file
  TFile outputFile_;

  /// Non-owning pointer to the output tree
  TTree *tree_;

  /// Buffer to save the nominal event weight
  Float_t weight_;

  /**
   * \brief Buffers to save alternative event weights
   *
   * These are full weights, i.e. they are not relative with respect to the
   * nominal one.
   */
  std::vector<Float_t> systWeights_;
};

#endif  // HZZ2L2NU_INCLUDE_EVENTTREES_H_

