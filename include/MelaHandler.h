#ifndef MELAHANDLER_H
#define MELAHANDLER_H

#include <filesystem>
#include <Mela.h>


namespace fs = std::filesystem;


class SmartMela {
 public:
   SmartMela(double sqrts, double smHiggsMass, TVar::VerbosityLevel vl)
       : mela{sqrts, smHiggsMass, vl} {}

   ~SmartMela() {
     for (auto const &path : {"ffwarn.dat", "br.sm1", "br.sm2", "input.DAT",
         "process.DAT", "Pdfdata"})
       fs::remove_all(path);
   }

   Mela* Get() {
     return &mela;
   }

 private:
   Mela mela;
};


// This class will be implemented in Mela Reweighting MR.
/**
 * \brief Computes the weights for MELA reweighing process. 
 *
 * Weights are computed using
 * <a href="https://github.com/JHUGen/JHUGenMELA/tree/master/MELA">MELA</a> and
 * <a href="https://github.com/MELALabs/MelaAnalytics">MelaAnalyctics</a>
 * packags based LHE level information.
 * Computed weights are cached on a per-event basis.
 */
class MelaHandler {
  public:
    MelaHandler();

    /// The Mela object is the core of MELA weights computation
    static SmartMela smartMela_;

};

#endif
