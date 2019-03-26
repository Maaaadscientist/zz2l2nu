#ifndef COLLECTIONBUILDER_H_
#define COLLECTIONBUILDER_H_

#include <boost/iterator/iterator_facade.hpp>

#include <TLorentzVector.h>


/**
 * \brief Base class for classes that construct collections of physics objects
 *
 * The primary job of this class is to provide a convenient interface to iterate
 * over four-momenta of the physics objects in the collection produce by the
 * derived class. In order to do this, the derived class must implement methods
 * \ref GetMomentum and \ref GetNumMomenta. If the derived class provides
 * several collections, it must choose one.
 */
class CollectionBuilder {
 public:
  using Momentum = TLorentzVector;

  /// A constant iterator for MomentaWrapper
  class MomentumIt : public boost::iterator_facade<
      MomentumIt, Momentum const, boost::random_access_traversal_tag,
      Momentum const &, int> {
   public:
    MomentumIt() = default;
    MomentumIt(CollectionBuilder const &builder, size_t index);

   private:
    friend class boost::iterator_core_access;

    void advance(difference_type n);
    void decrement();
    reference dereference() const;
    difference_type distance_to(MomentumIt const &other) const;
    bool equal(MomentumIt const &other) const;
    void increment();

    CollectionBuilder const &builder_;
    size_t index_;
  };

  /**
   * \brief Facade to access elements of the collection of momenta
   *
   * This auxiliary class is a syntactic sugar to iterate through momenta
   * provided by CollectionBuilder. It mimics a selection of methods of
   * std::vector that provide read-only access. In addition, some convenience
   * methods motivated by physics are provided.
   */
  class MomentaWrapper {
   public:
    MomentaWrapper(CollectionBuilder const &builder);

    Momentum const &at(size_t index) const;
    MomentumIt begin() const;
    MomentumIt end() const;
    Momentum const &operator[](size_t index) const;
    size_t size() const;

    /// Checks for overlaps with the given momentum in the (eta, phi) metric
    bool HasOverlap(Momentum const &p4, double maxDR) const;

   private:
    CollectionBuilder const &builder_;
  };

  /// Provides access to four-momenta of physics objects in the collection
  MomentaWrapper GetMomenta() const;

 private:
  /// Interface to access momentum of the object with the given index
  virtual Momentum const &GetMomentum(size_t i) const = 0;

  /// Interface to access the size of the collection
  virtual size_t GetNumMomenta() const = 0;
};


inline CollectionBuilder::MomentumIt::MomentumIt(
    CollectionBuilder const &builder, size_t index)
    : builder_{builder}, index_{index} {}


inline void CollectionBuilder::MomentumIt::advance(difference_type n) {
  index_ += n;
}


inline void CollectionBuilder::MomentumIt::decrement() {
  --index_;
}


inline CollectionBuilder::MomentumIt::reference
CollectionBuilder::MomentumIt::dereference() const {
  return builder_.GetMomentum(index_);
}


inline CollectionBuilder::MomentumIt::difference_type
CollectionBuilder::MomentumIt::distance_to(MomentumIt const &other) const {
  return int(other.index_) - int(index_);
}


inline bool CollectionBuilder::MomentumIt::equal(
    MomentumIt const &other) const {
  return &builder_ == &other.builder_ and index_ == other.index_;
}


inline void CollectionBuilder::MomentumIt::increment() {
  ++index_;
}


inline CollectionBuilder::MomentaWrapper::MomentaWrapper(
    CollectionBuilder const &builder)
    : builder_{builder} {}


inline CollectionBuilder::Momentum const &
CollectionBuilder::MomentaWrapper::at(size_t index) const {
  return builder_.GetMomentum(index);
}


inline CollectionBuilder::MomentumIt
CollectionBuilder::MomentaWrapper::begin() const {
  return {builder_, 0};
}


inline CollectionBuilder::MomentumIt
CollectionBuilder::MomentaWrapper::end() const {
  return {builder_, builder_.GetNumMomenta()};
}


inline CollectionBuilder::Momentum const &
CollectionBuilder::MomentaWrapper::operator[](size_t index) const {
  return builder_.GetMomentum(index);
}


inline size_t CollectionBuilder::MomentaWrapper::size() const {
  return builder_.GetNumMomenta();
}


inline CollectionBuilder::MomentaWrapper CollectionBuilder::GetMomenta() const {
  return {*this};
}

#endif  // COLLECTIONBUILDER_H_

