#ifndef PHYSICSOBJECTS_H_
#define PHYSICSOBJECTS_H_

#include <limits>

#include <TLorentzVector.h>


/// Non-polymorphic base class for particle-like physics objects
struct Particle {

  /// Four-momentum, in GeV
  TLorentzVector p4;
};


/**
 * \brief Comparison function to sort collections of particle-like objects in pt
 *
 * \param[in] p1, p2  Particle-like objects.
 * \return True if the pair of given particle-like objects is ordered in pt,
 *   meaning that pt of the first object is greater than pt of the second one.
 */
inline bool PtOrdered(Particle const &p1, Particle const &p2) {
  return (p1.p4.Pt() > p2.p4.Pt());
}


/// Generator-level particle
struct GenParticle : public Particle {

  /// Constructor from PDG ID
  GenParticle(int pdgId) noexcept;

  /// PDG ID code
  int pdgId;
};


inline GenParticle::GenParticle(int pdgId_) noexcept
    : Particle{}, pdgId{pdgId_} {}


/// Generator-level jet
struct GenJet : public Particle {};


/// Reconstructed jet
struct Jet : public Particle {
  
  /// Default constructor
  Jet() noexcept;

  /**
   * \brief Value of CSVv2 b-tagging discriminator
   *
   * NaN if not set.
   */
  double bTag;

  /**
   * \brief Hadron flavour of the jet
   *
   * Possible values are 0, 4, and 5.
   */
  int hadronFlavour;
};


inline Jet::Jet() noexcept
    : Particle{}, bTag{std::numeric_limits<double>::quiet_NaN()},
      hadronFlavour{0} {}


/// Reconstructed missing pt
struct PtMiss : public Particle {

  /// Default constructor
  PtMiss() noexcept;

  /**
   * \brief Value of ptmiss significance
   *
   * NaN if not set.
   */
  double significance;
};


inline PtMiss::PtMiss() noexcept
    : Particle{}, significance{std::numeric_limits<double>::quiet_NaN()} {}


/// Reconstructed photon
struct Photon : public Particle {

  /// Default constructor
  Photon() noexcept;
};

inline Photon::Photon() noexcept
    : Particle{} {}


/// Reconstructed charged lepton
struct Lepton : public Particle {

  /// Lepton flavour
  enum class Flavour {
    Electron,
    Muon
  };

  /// Constructor from flavour
  Lepton(Flavour flavour) noexcept;

  /// Flavour
  Flavour flavour;

  /**
   * \brief Electric charge
   *
   * Allowed values are +-1 and 0, the latter meaning that the charge has not
   * been specified.
   */
  int charge;

  /**
   * \brief Uncorrected momentum
   *
   * If momentum of the lepton is corrected, its old value is expected to be
   * stored in this field. This is useful to access lepton scale factors using
   * the old momentum. If the momentum has not been corrected, set to a null
   * vector.
   */
  TLorentzVector uncorrP4;
};


inline Lepton::Lepton(Flavour flavour_) noexcept
    : Particle{}, flavour{flavour_}, charge{0}, uncorrP4{} {}


/// Reconstructed electron
struct Electron : public Lepton {

  /// Default constructor
  Electron() noexcept;

  /**
   * \brief Pseudorapidity of associated ECAL supercluster
   *
   * NaN if not set.
   */
  double etaSc;
};


inline Electron::Electron() noexcept
    : Lepton{Lepton::Flavour::Electron},
      etaSc{std::numeric_limits<double>::quiet_NaN()} {}


/// Reconstructed muon
struct Muon : public Lepton {

  /// Default constructor
  Muon() noexcept;
};


inline Muon::Muon() noexcept
    : Lepton{Lepton::Flavour::Muon} {}

#endif  // PHYSICSOBJECTS_H_

