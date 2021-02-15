#ifndef VBFDISCRIMINANT_H
#define VBFDISCRIMINANT_H

#include <TSpline.h>

#include <Mela.h>
#include <GMECHelperFunctions.h>

#include <JetBuilder.h>
#include <Options.h>


/**
 * \brief Computes the VBF discriminant. 
 *
 * VBF discriminant is computed using
 * <a href="https://github.com/JHUGen/JHUGenMELA/tree/master/MELA">MELA</a> and
 * <a href="https://github.com/MELALabs/MelaAnalytics">MelaAnalyctics</a>
 * packags based on Reco information.
 * Discriminants are cached on a per-event basis.
 */
class VBFDiscriminant {
  public:
    // Constructor
    VBFDiscriminant(Options const &options);

    // Destructor
    ~VBFDiscriminant();

    enum MEP : int {
      kSIGa1 = 0,
      kSIGa2 = 1,
      kSIGa3 = 2,
      kSIGl1 = 3,
      kALT = 4
    };

    enum DjjVBF : int {
      a1 = 0,
      SM = 0,
      a2 = 1,
      a3 = 2,
      L1 = 3
    };

    /**
     * \brief Provides the DjjVBF discriminants w.r.t a1 (SM), a2 and a3.
     *
     * \param[in] p4LL      Four momenta of reconstructed Dilepton system.
     * \param[in] p4Miss    Four momentum of missing pt.
     * \param[in] jets      Reconstructed jets with tight ID.
     */
    std::array<double, 4> const &Get(
        TLorentzVector const &p4LL,
        TLorentzVector const &p4Miss,
        std::vector<Jet> const &jets);

  private:
    /// Resets Mela settings for next event DjjVBF discriminants computation.
    void Reset();

    /// Builds Mela candidate based on the reconstructed Physics Objects.
    void BuildMelaCandidate(TLorentzVector const &p4LL,
        TLorentzVector const &p4Miss, std::vector<Jet> const &jets);

    /// Computes Mela clusters
    void ComputeClusters() const;
    
    /// Facilitates accessing to the Mela object
    Mela *melaHandle_;

    /// Hypotheses registered as requested in the configuration file
    std::vector<std::unique_ptr<MELAHypothesis>> hypotheses_;

    /// A vector of non-owning pointers to hypotheses tagged as Aliased
    std::vector<MELAHypothesis*> aliasedHypos_;

    /// Stores computations registered as requested in the configuration file
    std::vector<std::unique_ptr<MELAComputation>> computers_;

    /// Stores clusters registered as requested in the configuration file
    std::vector<MELACluster*> clusters_;
    
    /// Stores Approximate 4-momentum of ZZ system
    TLorentzVector p4ZZApprox_;

    /// Stores DjjVBF discriminants for a1 (SM), a2 and a3 couplings
    std::array<double, 4> dJJVBF_;

    /// Stores MEs flags
    std::array<std::string, 5> mepFlag_;

    /// Stores g and c constants values used in Djj VBF discriminant formula
    std::unique_ptr<TSpline3> cConstant_;
    std::array<std::unique_ptr<TSpline3>, 4> gConstant_;
};

#endif
