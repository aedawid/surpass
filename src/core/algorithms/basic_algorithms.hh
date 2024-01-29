/** @file basic_algorithms.hh
 * @brief provides some basic algorithms that are not defined by the STL
 *
 */
#ifndef CORE_ALGORITHMS_basic_algorithms_H
#define CORE_ALGORITHMS_basic_algorithms_H

#include <vector>
#include <iostream>
#include <algorithm>

#include <core/index.hh>

namespace core {
namespace algorithms {

/** @brief Binary search routine that returns the index of an element rather than an interator.
 *
 * In general, this method returns an index of the lower bound for a requested element.
 *     - if the given element is smaller than <code>container</code>, returns -1
 *     - if the given element is smaller than <code>container.back()</code>, returns container.size()
 *     - if the given element equals to <code>container[i]</code>, returns i
 *     - if the given element should be inserted between to <code>container[i]</code> and <code>container[i+1]</code>, returns i
 * @tparam E - the type of the elements
 * @tparam T - data container type, e.g. <code>std::vector<E></code>
 * @param container - container to be searched through
 * @param element - what are we looking for
 * @return the index of the lower bound for the searched element
 */
template<typename T, typename E>
static core::index4 binary_search_index(const T& container, const E element) {

  if (container[0] >= element) {
    if (container[0] == element) return 0;
    else return -1;
  }
  if (container.back() <= element) {
    if (container.back() > element) return container.size();
    else return container.size() - 1;
  }
  core::index4 klo = 0;
  core::index4 khi = container.size() - 1;
  while (khi - klo > 1) {
    core::index4 k = (khi + klo) >> 1;
    if (container[k] > element) khi = k;
    else klo = k;
  }

  return klo;
}

/** @brief Removes duplicates from a given range
 *
 * @param begin - iterator pointing to the first element of the range
 * @param end - pass-the-end iterator marking the end of the range
 * @tparam It - iterator type
 */
template<typename It>
static It uniquify(It begin, It const end) {

  std::vector<It> v;
  v.reserve(static_cast<size_t>(std::distance(begin, end)));
  for (It i = begin; i != end; ++i) {
    v.push_back(i);
  }
  std::sort(v.begin(), v.end(), [](const It & i1,const It & i2) {return *i1<*i2;});
  v.erase(std::unique(v.begin(), v.end(), [](const It & i1,const It & i2) {return *i1==*i2;}), v.end());
  std::sort(v.begin(), v.end());
  size_t j = 0;
  for (It i = begin; i != end && j != v.size(); ++i) {
    if (*i == *v[j]) {
      using std::iter_swap;
      iter_swap(i, begin);
      ++j;
      ++begin;
    }
  }

  return begin;
}

/** @brief Removes duplicates from a given range
 *
 * @param begin - iterator pointing to the first element of the range
 * @param end - pass-the-end iterator marking the end of the range
 * @param is_equal - tests whether two elements are equal
 * @param less_than - tests whether one element is less than another
 * @tparam It - iterator type
 * @tparam EqOp - is-equal operator type
 * @tparam LtOp - less-than operator type
 */
template<typename It, typename EqOp, typename LtOp>
static It uniquify(It begin, It const end, EqOp is_equal, LtOp less_than) {

  std::vector<It> v;
  v.reserve(static_cast<size_t>(std::distance(begin, end)));
  for (It i = begin; i != end; ++i) {
    v.push_back(i);
  }
  std::sort(v.begin(), v.end(), less_than);
  v.erase(std::unique(v.begin(), v.end(), is_equal), v.end());
  std::sort(v.begin(), v.end());
  size_t j = 0;
  for (It i = begin; i != end && j != v.size(); ++i) {
    if (is_equal(i, v[j])) {
      using std::iter_swap;
      iter_swap(i, begin);
      ++j;
      ++begin;
    }
  }
  return begin;
}

/** @brief Find intersection of two sorted ranges.
 *
 * Duplicates are allowed. Usage of this method is very simple:
 * \include ex_intersect_sorted.cc
 *
 * @param begin1 - iterator pointing to the first element of the first range
 * @param end1 - pass-the-end iterator marking the end of the first range
 * @param begin2 - iterator pointing to the first element of the second range
 * @param end2 - pass-the-end iterator marking the end of the second range
 * @param sink_container - intersection of the two ranges (i.e. the common elements) will be pushed back to this container
 */
template<typename It, typename Out>
void intersect_sorted(It  begin1, It  end1, It  begin2, It  end2, Out & sink_container) {

  It it1 = begin1, it2 = begin2;

  //while either of the two indices reaches end
  while (it1 != end1 && it2 != end2) {
    //if first array element is lesser, advance that index by one
    if (*it1 < *it2) ++it1;
    // otherwise advance second index
    else {
      if (*it1 > *it2) ++it2;
      //both elements are same, print it, and advance both the pointers
      else {
        sink_container.push_back(*it1);
        ++it1;
        ++it2;
      }
    }
  }
}

}
}

/**
 * \example ex_intersect_sorted.cc
 */
#endif
