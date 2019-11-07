#ifndef COLLECTIONBUILDER_H_
#define COLLECTIONBUILDER_H_

#include <initializer_list>
#include <vector>

#include <boost/iterator/iterator_facade.hpp>
#include <TLorentzVector.h>
#include <TTreeReader.h>

#include <EventCache.h>


/**
 * \brief Base for classes to construct collections of physics objects
 *
 * This class provides an interface to access and iterate over four-momenta of
 * physics objects in the collection produced by a derived class. For this, the
 * derived class must implement methods \ref GetMomentum and \ref GetNumMomenta.
 *
 * This class also provides a mechanism for per-event caching and lazy
 * construction of the collection of physics objects, see methods \ref Update
 * and \ref Build.
 *
 * Methods \ref EnableCleaning and \ref IsDuplicate are helpful to avoid double
 * counting of objects with collections produced by other builders. A derived
 * class must make use of method \ref IsDuplicate to skip objects identified as
 * duplicates.
 *
 * When a derived class changes momenta of physics objects in the collection, it
 * should register the change using method \ref AddMomentumShift. The change in
 * the sum momentum aggregated in an event is provided with method
 * \ref GetSumMomentumShift, which is used by PtMissBuilder to propagate the
 * change into missing pt.
 */
class CollectionBuilderBase {
 public:
  using Momentum = TLorentzVector;

  /// A constant iterator for MomentaWrapper
  class MomentumIt : public boost::iterator_facade<
      MomentumIt, Momentum const, boost::random_access_traversal_tag,
      Momentum const &, int> {
   public:
    MomentumIt() = default;
    MomentumIt(CollectionBuilderBase const &builder, size_t index);

   private:
    friend class boost::iterator_core_access;

    void advance(difference_type n);
    void decrement();
    reference dereference() const;
    difference_type distance_to(MomentumIt const &other) const;
    bool equal(MomentumIt const &other) const;
    void increment();

    CollectionBuilderBase const &builder_;
    size_t index_;
  };

  /**
   * \brief Facade to access elements of the collection of momenta
   *
   * This auxiliary class is a syntactic sugar to iterate through momenta
   * provided by CollectionBuilderBase. It mimics a selection of methods of
   * std::vector that provide read-only access. In addition, some convenience
   * methods motivated by physics are provided.
   */
  class MomentaWrapper {
   public:
    MomentaWrapper(CollectionBuilderBase const &builder);

    Momentum const &at(size_t index) const;
    MomentumIt begin() const;
    MomentumIt end() const;
    Momentum const &operator[](size_t index) const;
    size_t size() const;

    /// Checks for overlaps with the given momentum in the (eta, phi) metric
    bool HasOverlap(Momentum const &p4, double maxDR) const;

   private:
    CollectionBuilderBase const &builder_;
  };

  /**
   * \brief Constructor
   *
   * The reader object provided as the argument is used to implement per-event
   * caching.
   */
  CollectionBuilderBase(TTreeReader &reader);

  /**
   * \brief Requests cleaning with respect to collections produced by given
   * builders
   *
   * It is up to a derived class to actually remove duplicates identified by
   * method \ref IsDuplicate. Objects provided in the argument are not owned by
   * this. If this method is called multiple times, new builders are added to
   * the list of already registered builders.
   */
  void EnableCleaning(
    std::initializer_list<CollectionBuilderBase const *> builders);
  
  /// Provides access to four-momenta of physics objects in the collection
  MomentaWrapper GetMomenta() const;

  /**
   * \brief Returns the accumulated change in momentum in the current event,
   * introduced as a result of calibration and other changes to momenta of
   * individual objects.
   */
  Momentum const &GetSumMomentumShift() const;

  /**
   * \brief Checks if the direction of the given momentum is close to that of a
   * momentum in one of the collections for cleaning.
   *
   * The matching is done in the (eta, phi) metric.
   */
  bool IsDuplicate(Momentum const &p4, double maxDR) const;

 protected:
  /**
   * \brief Adds a shift to the sum four-momentum
   *
   * When a derived class changes momentum of a physics object as a result of an
   * additional calibration or a systematic variation, it should register the
   * change with the help of this method.
   *
   * The accumulated shift is intended to be used to correct missing pt (see
   * PtMissBuilder), and the objects that contribute to it should be chosen
   * accordingly. These are not necessarily the same objects as included in the
   * produced collection.
   */
  void AddMomentumShift(Momentum const &uncorrP4, Momentum const &corrP4) const;

  /**
   * \brief Makes sure that all construction for the current event has been
   * performed
   *
   * When this method is called for the first time for a given event, it calls
   * \ref Build. Otherwise it does nothing. A derived class must call this
   * method before it attempts to return the collection of physics objects in
   * the current event.
   */
  void Update() const;

 private:
  /**
   * \brief A hook to perform construction for the current event
   *
   * This method will be called at maximum once per event. Its implementation in
   * the derived class should construct the collection of physics objects.
   */
  virtual void Build() const = 0;

  /// Interface to access momentum of the object with the given index
  virtual Momentum const &GetMomentum(size_t i) const = 0;

  /// Interface to access the size of the collection
  virtual size_t GetNumMomenta() const = 0;

  /// An object to facilitate caching
  EventCache cache_;

  /**
   * \brief Collection of non-owning pointers to objects that produce
   * collections against which the cleaning should be performed.
   */
  std::vector<CollectionBuilderBase const *> prioritizedBuilders_;

  /**
   * \brief Total change in four-momentum accumulated in the current event with
   * \ref AddP4Shift.
   */
  mutable Momentum sumP4Shift_;
};


inline CollectionBuilderBase::MomentumIt::MomentumIt(
    CollectionBuilderBase const &builder, size_t index)
    : builder_{builder}, index_{index} {}


inline void CollectionBuilderBase::MomentumIt::advance(difference_type n) {
  index_ += n;
}


inline void CollectionBuilderBase::MomentumIt::decrement() {
  --index_;
}


inline CollectionBuilderBase::MomentumIt::reference
CollectionBuilderBase::MomentumIt::dereference() const {
  return builder_.GetMomentum(index_);
}


inline CollectionBuilderBase::MomentumIt::difference_type
CollectionBuilderBase::MomentumIt::distance_to(MomentumIt const &other) const {
  return int(other.index_) - int(index_);
}


inline bool CollectionBuilderBase::MomentumIt::equal(
    MomentumIt const &other) const {
  return &builder_ == &other.builder_ and index_ == other.index_;
}


inline void CollectionBuilderBase::MomentumIt::increment() {
  ++index_;
}


inline CollectionBuilderBase::MomentaWrapper::MomentaWrapper(
    CollectionBuilderBase const &builder)
    : builder_{builder} {}


inline CollectionBuilderBase::Momentum const &
CollectionBuilderBase::MomentaWrapper::at(size_t index) const {
  return builder_.GetMomentum(index);
}


inline CollectionBuilderBase::MomentumIt
CollectionBuilderBase::MomentaWrapper::begin() const {
  return {builder_, 0};
}


inline CollectionBuilderBase::MomentumIt
CollectionBuilderBase::MomentaWrapper::end() const {
  return {builder_, builder_.GetNumMomenta()};
}


inline CollectionBuilderBase::Momentum const &
CollectionBuilderBase::MomentaWrapper::operator[](size_t index) const {
  return builder_.GetMomentum(index);
}


inline size_t CollectionBuilderBase::MomentaWrapper::size() const {
  return builder_.GetNumMomenta();
}


inline CollectionBuilderBase::CollectionBuilderBase(TTreeReader &reader)
    : cache_{reader} {}


inline CollectionBuilderBase::MomentaWrapper
CollectionBuilderBase::GetMomenta() const {
  return {*this};
}


inline CollectionBuilderBase::Momentum const &
CollectionBuilderBase::GetSumMomentumShift() const {
  Update();
  return sumP4Shift_;
}


inline void CollectionBuilderBase::AddMomentumShift(
    Momentum const &uncorrP4, Momentum const &corrP4) const {
  sumP4Shift_ += (corrP4 - uncorrP4);
}


inline void CollectionBuilderBase::Update() const {
  if (cache_.IsUpdated()) {
    sumP4Shift_ = Momentum{};
    Build();
  }
}


/**
 * \brief Intermediate base for classes to construct collections of physics
 * objects
 *
 * \tparam T  Type of physics objects in the collection.
 *
 * Provides implementation for methods to access momenta of physics objects
 * utilizing method \ref Get.
 */
template <typename T>
class CollectionBuilder : public CollectionBuilderBase {
 public:
  /**
   * \brief Constructor
   *
   * Directly forwards its argument to the base class.
   */
  CollectionBuilder(TTreeReader &reader);

  /// Interface to access the collection of physics objects
  virtual std::vector<T> const &Get() const = 0;

 private:
  /// Returns momentum of object with given index
  CollectionBuilderBase::Momentum const &GetMomentum(
    size_t index) const final;

  /// Returns the size of the collection
  size_t GetNumMomenta() const final;
};


template <typename T>
CollectionBuilder<T>::CollectionBuilder(TTreeReader &reader)
    : CollectionBuilderBase{reader} {}


template <typename T>
CollectionBuilderBase::Momentum const &
CollectionBuilder<T>::GetMomentum(size_t index) const {
  return Get().at(index).p4;
}


template <typename T>
size_t CollectionBuilder<T>::GetNumMomenta() const {
  return Get().size();
}

#endif  // COLLECTIONBUILDER_H_

