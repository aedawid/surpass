#ifndef SIMULATIONS_FORCEFIELDS_ByResidueEnergy_HH
#define SIMULATIONS_FORCEFIELDS_ByResidueEnergy_HH

#include <string>
#include <core/index.hh>

#include <simulations/forcefields/CalculateEnergyBase.hh>

namespace simulations {
namespace forcefields {

/** @brief Base class for an energy function that may be decomposed into per-residue interactions
 * Since this class inherits from CalculateEnergyBase, the derived classes must also implement
 * <code>CalculateEnergyBase.calculate()</code> method. It's not possible however to get energy introduced by a single atom.
 */
class ByResidueEnergy : public virtual CalculateEnergyBase {
public:

  /// Calculate energy of a given residue
  virtual double calculate_by_residue(const core::index2 which_residue) = 0;

  /** @brief Calculate energy of a given contiguous chunk of residues.
   *
   * This method assumes that residues [chunk_from, chunk_to] (both inclusive) have been moved in respect to the remaining part
   * of s system but also in respect to each other. This method therefore evaluates interactions within the chunk plus
   * interactions of the chunk with the remaining part of a system.
   * @param chunk_from - index of the first moved residue
   * @param chunk_to  - index of the last moved residue (inclusive)
   * @return energy of the chunk interacting with itself and with the remaining part of a system
   */
  virtual double calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to) = 0;

  /// Virtual destructor
  virtual ~ByResidueEnergy() { }
};

} // ~ simulations
} // ~ ff
#endif
