#ifndef CORE_DATA_SEQUENCE_SecondaryStructureAnnotation_H
#define CORE_DATA_SEQUENCE_SecondaryStructureAnnotation_H

#include <string>
#include <vector>
#include <memory>

#include <core/index.hh>
#include <core/real.hh>

namespace core {
namespace data {
namespace sequence {

/** @brief a helper structure to pass around all three probabilities in a single data structure
 */
struct HecFractions {

  const core::real p_H; ///< probability of the helical state
  const core::real p_E; ///< probability of the extended state
  const core::real p_C; ///< probability of the coil state

  /// Returns 'H', 'E' or 'C' symbol for which the probability is the highest.
  char ss() const {
    if (p_C >= p_E) {
      if (p_C >= p_H) return 'C';
    } else {
      if (p_E >= p_H) return 'E';
      else return 'H';
    }
    return 'C';
  }

  /** @brief Operator to facilitate array-style access to the probabilities.
   *
   * @param pos - must be either 0 (helix), 1 (extended) or 2 (coil). If <code> pos > 2 </code>, the method returns <code> p_C </code>
   */
  core::real operator[](const core::index2 pos) const { return (pos == 0) ? p_H : ((pos == 1) ? p_E : p_C); }
};

/** @brief Class that holds secondary structure probabilities.
 *
 * For practical application one should use SecondaryStructure class which is derived from this one
 */
class SecondaryStructureAnnotation {
public:

  /** @brief Creates a secondary structure of a coiled chain.
   *
   * This secondary structure will have all secondary structure set to coil (C)
   * @param length - length of the annotate sequence
   */
  SecondaryStructureAnnotation(const core::index4 length) {
    ss_string.resize(length, 'C');
    h.resize(length, 0.0);
    e.resize(length, 0.0);
  }

  /** @brief Creates a secondary structure from a string-type definition.
   *
   * This secondary structure will have all probabilities equal to 1.0 or 0.0, according to the given string
   * @param ss_string - secondary structure as a string
   */
  SecondaryStructureAnnotation(const std::string & ss_string) : ss_string(ss_string) {

    h.resize(ss_string.size(), 0.0);
    e.resize(ss_string.size(), 0.0);
    for (core::index2 i = 0; i < ss_string.size(); i++) {
      if (ss_string[i] == 'H') h[i] = 1.0;
      if (ss_string[i] == 'E') e[i] = 1.0;
    }
  }

  /// Bare virtual destructor to satisfy a compiler
  virtual ~SecondaryStructureAnnotation() {}

  /** @brief Says whether this object bears information about secondary structure for this sequence.
   * @return always true
   */
  virtual bool has_ss() const { return true; }

  /** @brief returns the character denoting the most popular secondary structure type at a given position
   *
   * @param pos - residue index in the sequence
   * @returns the most popular  secondary structure type at that position (E, H or C)
   */
  inline char ss(const core::index2 pos) const { return ss_string[pos]; }

  ///< Returns the secondary structure as a string (containing H, E and C)
  inline const std::string & str() const { return ss_string; }

  ///< Returns the probability of finding a given position in the helical state
  inline core::real fraction_H(const core::index2 pos) const { return h[pos]; }

  ///< Returns the probability of finding a given position in the beta state
  inline core::real fraction_E(const core::index2 pos) const { return e[pos]; }

  ///< Returns the probability of finding a given position in a loop
  inline core::real fraction_C(const core::index2 pos) const { return real(1.0 - h[pos] - e[pos]); }

  ///< Returns the probability to find a certain residue in each of the three secondary structure type
  inline const HecFractions fractions(const core::index2 pos) const {
    return HecFractions({h[pos], e[pos], 1.0 - h[pos] - e[pos]});
  }

  ///< Sets the probability of finding a given position in the alpha, beta and coil state
  inline void fractions(const core::index2 pos, const HecFractions & f) { fractions(pos, f.p_H,f.p_E,f.p_C); }

  ///< Sets the probability of finding a given position in the alpha, beta and coil state
  inline void fractions(const core::index2 pos, const core::real new_fraction_H, const core::real new_fraction_E,
      const core::real new_fraction_C) {
    const core::real sum = new_fraction_H + new_fraction_E + new_fraction_C;
    if (sum != 1.0) {
      h[pos] = new_fraction_H / sum;
      e[pos] = new_fraction_E / sum;
    } else {
      e[pos] = new_fraction_E;
      h[pos] = new_fraction_H;
    }
    update_ss(pos);
  }

private:
  std::vector<core::real> h;
  std::vector<core::real> e;
  std::string ss_string;

  inline void update_ss(const core::index2 pos) {
    if (h[pos] + e[pos] < 0.5) ss_string[pos] = 'C';
    else if (h[pos] > e[pos]) ss_string[pos] = 'H';
    else ss_string[pos] = 'E';
  }
};

}
}
}

#endif
