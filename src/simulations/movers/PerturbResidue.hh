#ifndef SIMULATIONS_CARTESIAN_MOVERS_PerturbResidue_HH
#define SIMULATIONS_CARTESIAN_MOVERS_PerturbResidue_HH

#include <random>
#include <memory>

#include <core/real.hh>
#include <core/calc/statistics/Random.hh>

#include <utils/Logger.hh>

#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/movers/Mover.hh>
#include <simulations/sampling/AbstractAcceptanceCriterion.hh>

#include <simulations/systems/ResidueChain.hh>

namespace simulations {
namespace movers {

/** @brief Moves atoms of a randomly selected residue by a random vector in the Cartesian space.
 *
 * This mover will displace all atoms from a randomly selected residue by shifting each atom independently by
 * a random vector. The following example creates a random polymer chain and makes 100 residue perturbation moves:
 *
 * @tparam C - the type used to express coordinates; Vec3 for simple models, Vec3Cubic for modeling in periodic boundary conditions
 * \include ex_PerturbResidue.cc
 */
template<class C>
class PerturbResidue : public simulations::movers::Mover {
public:
  /** @brief Create a mover for a given system.
   *
   * @param system - a ResidueChain<C> object to be altered by this mover
   * @param energy - energy function used to accept / deny a move
   */
  PerturbResidue(systems::ResidueChain<C> & system, forcefields::ByResidueEnergy & energy);

  /** @brief Make MC moves.
   *
   * This method calls <code>move(simulations::generic::sampling::AbstractAcceptanceCriterion &)</move> method
   * multiple times.
   * @param n_moves - how many moves to make
   * @param mc_scheme - an acceptance criterion, used to accept or deny a move based on a system's energy before and after the move
   */
  core::index2 move(const core::index2 n_moves, simulations::sampling::AbstractAcceptanceCriterion & mc_scheme);

  void max_move_range(core::real step);

  /** @brief Make MC moves.
   *
   * This method:
   *     - selects a residue randomly (uniform distribution)
   *     - for each atom of the residue, it creates a random vector and move the atom by the vector
   *     - move is accepted (or cancelled) according to Monte Carlo criterion
   * @param mc_scheme - an acceptance criterion, used to accept or deny a move based on a system's energy before and after the move
   * @return true if the move got accepted; false otherwise
   */
  bool move(simulations::sampling::AbstractAcceptanceCriterion & mc_scheme);

  /** @brief Cancels the most recent move.
   *
   * This method is automatically called by move(simulations::generic::sampling::AbstractAcceptanceCriterion &) when
   * the respective AbstractAcceptanceCriterion said so. User doesn't need to call <code>undo()</code> by himself.
   */
  void undo();

  /// Returns the name of this mover
  virtual const std::string &  name() const { return name_; }

private:
  /// Maximum range of a move for each coordinate.
  core::real max_step_;
  int i_moved = -1;
  core::index4 last_moved_from = 0, last_moved_to = 0;
  systems::ResidueChain<C> & the_system;
  forcefields::ByResidueEnergy & the_energy;
  std::uniform_int_distribution<int> rand_residue_index;
  std::uniform_real_distribution<core::real> rand_coordinate;
  core::calc::statistics::Random & generator = core::calc::statistics::Random::get();
  std::unique_ptr<C[]> backup;
  static const std::string name_;
  utils::Logger logger;
};


}
}

#endif
/**
 * \example ex_PerturbResidue.cc
 */
