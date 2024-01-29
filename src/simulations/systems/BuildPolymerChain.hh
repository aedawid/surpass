#ifndef SIMULATIONS_SYSTEMS_BuildPolymerChain_HH
#define SIMULATIONS_SYSTEMS_BuildPolymerChain_HH

#include <memory>
#include <core/real.hh>
#include <core/index.hh>

namespace simulations {
namespace systems {

/** @brief Generates a random conformation of a simple polymer in the Cartesian space
 *
 * @tparam C - object used to represent a single atom, bead, particle, etc.
 */
template<class C>
class BuildPolymerChain {
public:
  //< Number of atoms in each chain to be generated
  const core::index4 n_atoms;

  /** @brief Constructs the chain generator for a fixed length of polymers.
   *
   * Generated coordinates will be stored in the given array. The new chain will start where <code>system[0]</code> is placed.
   * The first generated atom is placed at index 1.
   *
   * @param system - a unique pointer to an array where resulting coordinates will be written
   * @param n_atoms - the length of a polymer to be created.
   */
  BuildPolymerChain(std::unique_ptr<C[]>  & system, const core::index4 n_atoms);

  /** @brief Generates a new polymer conformation.
   *
   * This method places randomly each atom <code>bond_length</code> apart from the previous atom and tests excluded
   * volume test of the new atom with the upstream part of a chain. The procedure is repeated until
   * the newly created atom does not violate the excluded volume test. This might be repeated up to <code>n_bead_attempts</code>
   * times. If still failing, a new chain is started over. After <code>n_chain_attempts</code> unsuccessful chain attempts,
   * <code>false</code> is returned
   * @param bond_length - length of each bond between beads, i.e. the distance between subsequent atoms.
   * @param cutoff - excluded volume distance; any two atoms cannot get closer than this value.
   * @param n_bead_attempts - how many times placement of each bead is repeated. After that, this method starts a new chain from scratch.
   * @param n_chain_attempts - how many times a chain generation procedure is repeated. After that, false is returned.
   */
  bool generate(const core::real bond_length, const core::real cutoff,
      const core::index2 n_bead_attempts = 1000, const core::index2 n_chain_attempts = 1000);

private:
  core::real box_width;
  C tmp;
  std::unique_ptr<C[]>  & system;
  std::uniform_real_distribution<core::real> rand_coordinate;
  core::calc::statistics::Random & generator;

  bool try_chain(const core::real bond_length, const core::index2 n_attempts, const core::real cutoff);
  bool is_good_point(const core::index4 n_atoms_so_far, const C & candidate, const core::real min_distance_square);
};

}
}

#endif
