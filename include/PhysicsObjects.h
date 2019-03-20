#ifndef PHYSICSOBJECTS_H_
#define PHYSICSOBJECTS_H_

#include <limits>

#include <TLorentzVector.h>


/// Non-polymorphic base class for particle-like physics objects
struct Particle {

  /// Four-momentum, in GeV
  TLorentzVector p4;
};


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
struct Jet : public Particle {};


/// Missing pt
struct PtMiss : public Particle {};


/// Reconstructed photon
struct Photon : public Particle {};


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
};


inline Lepton::Lepton(Flavour flavour_) noexcept
    : Particle{}, flavour{flavour_}, charge{0} {}


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

