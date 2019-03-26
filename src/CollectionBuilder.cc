#include <CollectionBuilder.h>

#include <cmath>

#include <Utils.h>


bool CollectionBuilder::MomentaWrapper::HasOverlap(
    Momentum const &p4, double maxDR) const {
  double const maxDR2 = std::pow(maxDR, 2);

  for (auto const &curP4 : *this) {
    double const dR2 = utils::DeltaR2(curP4, p4);

    if (dR2 < maxDR2)
      return true;
  }

  return false;
}

