#ifndef CORE_ALORITHMS_UnionFind_HH
#define CORE_ALORITHMS_UnionFind_HH

#include <map>
#include <stdexcept>

#include <core/index.hh>
#include <utils/string_utils.hh>

namespace core {
namespace algorithms {

/** \brief The class provides union-find data container for integer data (a disjoint-set data structure).
 *
 * <p>With the container it is possible to efficiently maintain partition of the data into a number of separate, non-overlapping sets.</p>
 *
 * <p>A union-find algorithm performs two useful operations on such a data structure:</p>
 * <dl>
 * <dt><code>find_set()</code></dt><dd> Determine which set a particular element is in. Also useful for determining if two elements are in the same set.</dd>
 * <dt><code>union_set</code></dt><dd> combine or merge two sets into a single set.</dd>
 * </dl>
 * <p>It is also possible to check, how many elements belong to a given group and to retrieve all groups of elements.</p>
 * @tparam T - the type of data to attach to each element (arbitrary)
 * @tparam I - integer type used for indexing elements; use core::index1, core::index2 or core::index4 for this purpose
 */
template<typename T, typename I>
class UnionFind {
public:

  /// Create an empty container
  UnionFind() : set_count_(0) { }

  /** @brief Create an empty container with a given initial capacity.
   * This constructor allocates memory for the internal data structures so they may hold the given number of elements
   */
  UnionFind(I capacity) : set_count_(0) {

    elements.reserve(capacity);
    parents.reserve(capacity);
    ranks.reserve(capacity);
  }

  /** @brief  Adds a single element.
   * Constant time operation
   * @param value - element to be inserted in the set
   */
  I add_element(const T &value) {

    I i = elements.size();
    elements.push_back(value);
    index_4_element[value] = i;
    parents.push_back(i);
    ranks.push_back(1);
    ++set_count_;

    return i;
  }

  /// Returns the number of objects maintained by this container
  size_t count_elements() const { return elements.size(); }

  /** @brief  Returns the current number of disjoint sets in the forest (i.e. the current number of trees).
   */
  I count_sets() const { return set_count_; }

  /** @brief  Returns the size of a disjoint set the given <code>element_index</code> belongs to.
   * @param element_index - index of an element; does not have to be the root of its disjoint set (tree)
   * @return size of a disjoint set of elements represented by a given element index
   */
  I set_size(const I element_index) const { return ranks[find_set(element_index)]; }

  /// Clears the container. After that it contains no elements.
  void clear() {

    elements.clear();
    parents.clear();
    index_4_element.clear();
    ranks.clear();
    set_count_ = 0;
  }

  /// Clears all the connections between elements, but does not actually remove them
  void disconnect() {
    set_count_ = elements.size();
    for(I i =0;i<set_count_;++i) {
      ranks[i] = 1;
      parents[i] = i;
    }
  }

  /** @brief  Returns the index of the root element of the disjoint tree <code>element_index</code> belongs to.
   * @param index - the element whose set should be determined
   * @return index of the root element representing a set of elements
   */
  I find_set(const I index) const {
    I parent = parents[index];
    if (parent != index) { parent = find_set(parent); }
    return parent;
  }

  /** @brief Merges two disjoint sets containing elements index_i and index_j.
   * @param index_i - index of the first element
   * @param index_j - index of the second element
   * @returns true if the two elements or sets were merged; false only if they already belonged to the same cluster
   */
  bool union_set(const I index_i, const I index_j) {

    index4 set_i = find_set(index_i);
    index4 set_j = find_set(index_j);
    if (set_i == set_j) return false;
    if (ranks[set_i] > ranks[set_j]) {
      ranks[set_i] += ranks[set_j];
      parents[set_j] = set_i;
    } else {
      ranks[set_j] += ranks[set_i];
      parents[set_i] = set_j;
    }
    --set_count_;

    return true;
  }

  /** @brief Returns the object associated with a given index.
   * @param index - index referring to the requested element
   * @return requested element
   */
  T &element(I index) { return elements[index]; }

  /** @brief returns index pointing to a given element.
   *  @param element - an element that has been already inserted to this container
   *  @return an index of that element
   */
  I index(const T & element) { return index_4_element[element]; }

private:
  std::map<T, I> index_4_element;
  std::vector<T> elements;
  std::vector<I> parents;
  std::vector<I> ranks;
  I set_count_;
};

}
}

#endif // CORE_ALORITHMS_UnionFind_HH
