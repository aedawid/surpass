/** @file ShortRangeEnergyBase.hh
 * @brief provides ShortRangeEnergyBase base class
 */
#ifndef SIMULATIONS_FORCEFIELDS_ShortRangeEnergyBase_HH
#define SIMULATIONS_FORCEFIELDS_ShortRangeEnergyBase_HH

#include <vector>

#include <core/real.hh>
#include <core/index.hh>
#include <utils/Logger.hh>
#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/systems/ResidueChain.hh>

namespace simulations {
namespace forcefields {

using core::real;

/** @brief The base class for energy functions that assess the geometry of a chain fragment along its bonds.
 *
 * Provides a method to evaluate energy by chunk, using the per-residue energy provided by the derived class
 */
template<typename C>
class ShortRangeEnergyBase : public ByResidueEnergy {
public:

  ShortRangeEnergyBase(const systems::ResidueChain <C> &system, const core::index1 property_span) :
    property_span(property_span),last_positions_scored(system.count_residues() - property_span), the_system(system) { }

  virtual ~ShortRangeEnergyBase() { }

  /** @brief Number of residues required to calculate a single local property value.
   *
   * E.g. this is 2 for a protein backbone dihedral angle, 5 for \f$R_{15}\f$ distance etc.
   */
  const unsigned char property_span;
  /** @brief Index of the very last position this score may be applied to.
   * E.g. for \f$R_{15}\f$ distance this is \f$N_{res} - 5\f$ distance
   */
  const core::index2 last_positions_scored;

  /** @brief Evaluates energy of a contiguous fragment of a polymer chain.
   * This method assumes that beads within the chain may move in respect to each other and evaluates all
   * interactions within the fragment.
   * @param first - the first residue of the chunk
   * @param last - the last residue of the chunk (inclusive)
   */
  inline virtual double calculate_by_chunk(const residue_index first, const residue_index last) {

    double en = 0;
    const residue_index first_i = std::max(0, first - property_span);
    const residue_index last_i = std::min(last_positions_scored, last);

    for (residue_index i = first_i; i <= last_i; ++i)
      en += calculate_by_residue(i);

    return en;
  }

  /** @brief Evaluates energy of the whole chain.
 *
 * This method simply calls <code>calculate_by_rigid_chunk(0, last_positions_scored)</code>
 */
  inline virtual double calculate() { return calculate_by_chunk(0, last_positions_scored); }

protected:
  const systems::ResidueChain <C> &the_system;
};


}
}
#endif
