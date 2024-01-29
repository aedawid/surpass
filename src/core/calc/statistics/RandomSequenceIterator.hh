#ifndef CORE_CALC_STATISTICS_RandomSequenceIterator_HH
#define CORE_CALC_STATISTICS_RandomSequenceIterator_HH

#include <random>
#include <vector>
#include <iterator>

#include <core/calc/statistics/Random.hh>

namespace core {
namespace calc {
namespace statistics {

/** @brief Provides an iterator that iterates through a given vector in random order.
 *
 * This iterator is designed as a faster replacement for std::shuffle()
 */
template<typename T>
class RandomSequenceIterator : public std::iterator<std::forward_iterator_tag, T> {
public:

  /** @brief Returns <code>begin</code> iterator, pointing to the first element of the randomized sequence
   *
   * @param source - data to be iterated in a random order
   * @param sequence_size - how long should the sequence be? It might be actually longer or shorter than <code>source.size()</code>
   */
  static RandomSequenceIterator<T> begin(const std::vector<T> &source, int sequence_size) {
    return RandomSequenceIterator(source, sequence_size);
  }

  /// Provides the <code>end</code> iterator
  static RandomSequenceIterator<T> end() {  RandomSequenceIterator<T> e; return e;}

  /// Provides access to the data
  T operator*() const { return source_->operator[](which_element_); }

  /// Advance iterator to the next element
  RandomSequenceIterator<T> operator++() {
    if (sequence_size_ == 0) return *this;
    --sequence_size_;
    which_element_ = generator() % source_->size();
    return *this;
  }

  /// Postincrement operator
  RandomSequenceIterator<T> operator++(int) {
    RandomSequenceIterator<T> tmp = *this;
    ++*this;
    return tmp;
  }

  /// is-equal operator
  int operator==(const RandomSequenceIterator<T> &other) const {
    if (other.sequence_size_ == 0) return sequence_size_ == 0; //"end" iterator
    return source_ == other.source_;
  }

  /// is-not-equal operator
  int operator!=(const RandomSequenceIterator<T> &other) const { return !(*this == other); }

private:
  const std::vector<T> *source_;
  int sequence_size_;
  index2 which_element_;
  core::calc::statistics::Random &generator = core::calc::statistics::Random::get();

  //Creates "end" iterator
  RandomSequenceIterator() : source_(nullptr), sequence_size_(0) {}

  //Creates random "start" iterator
  RandomSequenceIterator(const std::vector<T> &source, index4 nOutputCount) :
    source_(&source), sequence_size_(nOutputCount + 1), m_distribution(0, source.size() - 1) {
    operator++(); //make new random value
  }

  std::uniform_int_distribution<index2> m_distribution;
};

}
}
}

#endif
