#include <CollectionBuilder.h>

#include <cmath>

#include <Utils.h>


bool CollectionBuilderBase::MomentaWrapper::HasOverlap(
    Momentum const &p4, double maxDR) const {
  double const maxDR2 = std::pow(maxDR, 2);

  for (auto const &curP4 : *this) {
    double const dR2 = utils::DeltaR2(curP4, p4);

    if (dR2 < maxDR2)
      return true;
  }

  return false;
}


void CollectionBuilderBase::EnableCleaning(
    std::initializer_list<CollectionBuilderBase const *> builders) {
  for (auto *b : builders)
    prioritizedBuilders_.emplace_back(b);
}


bool CollectionBuilderBase::IsDuplicate(Momentum const &p4, double maxDR) const {
  for (auto *builder : prioritizedBuilders_) {
    if (builder->GetMomenta().HasOverlap(p4, maxDR))
      return true;
  }
  return false;
}


