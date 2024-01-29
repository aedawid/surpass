/** @file Mover.hh
 * @brief Provides Mover base class
 */
#ifndef SIMULATIONS_GENERIC_MOVERS_Mover_HH
#define SIMULATIONS_GENERIC_MOVERS_Mover_HH

#include <memory> // for shared pointers

#include <core/real.hh>
#include <core/index.hh>
#include <simulations/sampling/AbstractAcceptanceCriterion.hh>

namespace simulations {
namespace movers {

/** Mover represents a Monte Carlo - type update of a system: modifies its coordinates and accepts the change according to a given criterion
 */
class Mover {
public:

  /** @brief Propose a move and accept it according to a given criterion.
   *
   * @param mc_scheme - acceptance criterion object
   * @return true if the move was successful, false otherwise
   */
  virtual bool move(sampling::AbstractAcceptanceCriterion & mc_scheme) = 0;

  /** @brief Attempts a number of MC moves
   *
   * @param n_moves - the number of moves to be attempted
   * @param mc_scheme - acceptance criterion object
   * @return the number of successful moves
   */
  virtual core::index2 n_moves(const core::index2 n_moves, sampling::AbstractAcceptanceCriterion & mc_scheme) {
    core::index2 ret = 0;
    for (core::index2 i = 0; i < n_moves; i++) ret += move(mc_scheme);
    return ret;
  }

  /// Steps back the most recent move
  virtual void undo() = 0;

  /// Returns the name of this mover, so the name may appear in the output when required
  virtual const std::string &  name() const = 0;

  /// Virtual destructor (empty)
  virtual ~Mover() { }

  /** @brief Advances both counters of attempetd and accepted moves by one
   * This method is called by <code>move()</code> to maintain the ratio of successfull moves
   */
  inline void inc_move_counter() {
    n_attempted++;
    n_successful++;
  }

  /** @brief Decrements the counter of successfull moves by one
   * This method is called by <code>undo()</code> to maintain the ratio of successfull moves
   */
  inline void dec_move_counter() {
    n_successful--;
  }

  /// Clears both counters of attempetd and accepted moves by one
  inline void clear_move_counter() {
    n_successful = 0;
    n_attempted = 0;
  }

  /** @brief Returns the success rate for this mover.
   *
   * @return success rate computed since the last reset of the counters
   */
  inline core::real get_success_rate() { return ((core::real) n_successful) / ((core::real) n_attempted); }

  /** @brief Returns the success rate for this mover and resets the counters.
   *
   * @return success rate computed since the last reset of the counters
   */
  inline core::real get_and_clear_success_rate() {

    core::real s = ((core::real) n_successful) / ((core::real) n_attempted);
    clear_move_counter();

    return s;
  }

private:

  int n_attempted = 0;
  int n_successful = 0;
};

/// Type representing a shared pointer to a Mover
typedef std::shared_ptr<Mover> Mover_SP;

}
}

#endif
