#include <WeightCollector.h>


WeightCollector::WeightCollector(
    std::initializer_list<WeightBase const *> computers)
    : computers_{computers} {}


void WeightCollector::Add(WeightBase const *computer) {
  computers_.emplace_back(computer);
}


double WeightCollector::NominalWeight() const {
  double weight = 1.;
  for (auto const c : computers_)
    weight *= c->NominalWeight();
  return weight;
}


int WeightCollector::NumVariations() const {
  int num = 0;
  for (auto const c : computers_)
    num += c->NumVariations();
  return num;
}


double WeightCollector::operator()() const {
  double weight = 1.;
  for (auto const c : computers_)
    weight *= (*c)();
  return weight;
}


double WeightCollector::RelWeight(int variation) const {
  auto const [index, localVariation] = TranslateIndex(variation);
  return computers_[index]->RelWeight(localVariation);
}


std::string_view WeightCollector::VariationName(int variation) const {
  auto const [index, localVariation] = TranslateIndex(variation);
  return computers_[index]->VariationName(localVariation);
}


std::tuple<int, int> WeightCollector::TranslateIndex(int variation) const {
  int compIndex = 0;
  int offset = 0;
  while (variation - offset >= computers_[compIndex]->NumVariations()) {
    offset += computers_[compIndex]->NumVariations();
    ++compIndex;
  }
  return {compIndex, variation - offset};
}

