#ifndef SIMULATIONS_CARTESIAN_MOVERS_PerturbChainFragment_HH
#define SIMULATIONS_CARTESIAN_MOVERS_PerturbChainFragment_HH

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

/** @brief Moves N beads of a polymer chain
 *
 * This mover will displace N subsequent beads
 *
 * @tparam C - the type used to express coordinates; Vec3 for simple models, Vec3Cubic for modeling in periodic boundary conditions
 */
template<class C>
class PerturbChainFragment : public simulations::movers::Mover {
public:
  /** @brief Create a mover for a given system.
   *
   * @param system - a ResidueChain<C> object to be altered by this mover
   * @param n_moved - the number of atom affected by a single move (constant for all moves)
   * @param energy - energy function used to accept / deny a move
   */
  PerturbChainFragment(systems::ResidueChain<C> &system, core::index2 n_moved, forcefields::ByResidueEnergy &energy);

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
   *     - selects a a fragment to be moved
   *     - for each atom of the residue, it creates a random vector and move the atom by the vector
   *     - move is accepted (or cancelled) according to Monte Carlo criterion
   * @param mc_scheme - an acceptance criterion, used to accept or deny a move based on a system's energy before and after the move
   * @return true if the move got accepted; false otherwise
   */
  virtual bool move(simulations::sampling::AbstractAcceptanceCriterion & mc_scheme);

  /** @brief Cancels the most recent move.
   *
   * This method is automatically called by move(simulations::generic::sampling::AbstractAcceptanceCriterion &) when
   * the respective AbstractAcceptanceCriterion said so. User doesn't need to call <code>undo()</code> by himself.
   */
  virtual void undo();

  /// Returns the name of this mover
  virtual const std::string &  name() const { return name_; }

private:
  /// Maximum range of a move for each coordinate.
  core::real max_step_;
  core::index2 n_moved_;
  core::index4 last_moved_from = 0, last_moved_to = 0;
  systems::ResidueChain<C> & the_system;
  forcefields::ByResidueEnergy & the_energy;
  std::uniform_int_distribution<core::index4> rand_bead_index;
  std::uniform_real_distribution<core::real> rand_coordinate;
  core::calc::statistics::Random & generator = core::calc::statistics::Random::get();
  std::unique_ptr<C[]> backup;
  static const std::string name_;
  utils::Logger logger;
};


}
}

#endif
