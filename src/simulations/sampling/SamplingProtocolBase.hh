#ifndef SIMULATIONS_GENERIC_SAMPLING_SamplingProtocolBase_HH
#define SIMULATIONS_GENERIC_SAMPLING_SamplingProtocolBase_HH

#include <vector>

#include <core/index.hh>

#include <simulations/evaluators/Evaluator.hh>
#include <simulations/observers/ObserverInterface.hh>
#include <simulations/sampling/ReplicaExchangeMC.fwd.hh>

namespace simulations {
namespace sampling {

/** @brief Base class for a sampling protocol handles observers and evaluators.
 *
 */
class SamplingProtocolBase {
public:

  /// Returns the number of outer (big) cycles of sampling
  core::index4 outer_cycles() const { return n_outer_cycles; }

  /// Returns the number of inner (small) cycles of sampling
  core::index4 inner_cycles() const { return n_inner_cycles; }

  /// Returns how many MC sweeps are done within each MC cycle - by default 1
  core::index4 cycle_size() const { return n_cycle_size; }

  /** @brief Set the number of outer cycles
   *  @param n_cycles - the number of outer (big) cycles
   */
  void outer_cycles(const core::index4 n_cycles) { n_outer_cycles = n_cycles; }

  /** @brief Set the number of inner cycles
   *  @param n_cycles - the number of inner (small) cycles
   */
  void inner_cycles(const core::index4 n_cycles) { n_inner_cycles = n_cycles; }

  /** @brief Set the number of MC sweeps per cycle
   *  @param n_sweeps - the number of  MC sweeps to be done per cycle
   */
  void cycle_size(const core::index4 n_sweeps) { n_cycle_size = n_sweeps; }

  /** @brief Set both the number of inner and outer cycles in a single call
   *  @param inner_cycles - the number of inner (small) cycles
   *  @param outer_cycles - the number of outer (big) cycles
   *  @param n_sweeps - the number of  MC sweeps to be done per cycle
   */
  void cycles(const core::index4 inner_cycles, const core::index4 outer_cycles, const core::index4 n_sweeps = 1) {

    n_inner_cycles = inner_cycles;
    n_outer_cycles = outer_cycles;
    n_cycle_size = n_sweeps;
  }

  /** @brief Adds a new Evaluator to be called after every inner cycle.
   * This evaluator will be called  \f$ N_T \times N_O \times N_I\f$ times
   * @param e - shared pointer to an object inheriting Evaluator interface
   */
  void inner_cycle_evaluator(evaluators::Evaluator_SP e) { evaluate_every_inner_cycle.push_back(e); }

  /** @brief Adds a new Evaluator to be called after every outer cycle.
   * This evaluator will be called  \f$ N_T \times N_O \f$ times
   * @param e - shared pointer to an object inheriting Evaluator interface
   */
  void outer_cycle_evaluator(evaluators::Evaluator_SP e) { evaluate_every_outer_cycle.push_back(e); }

  /** @brief Adds a new ObserverInterface instance to be called after every inner cycle.
 * This observer will be called  \f$ N_P \times N_O \times N_I\f$ times
 * @param o - shared pointer to an object inheriting Evaluator interface
 */
  void inner_cycle_observer(observers::ObserverInterface_SP o) { observe_every_inner_cycle.push_back(o); }

  /** @brief Adds a new ObserverInterface instance to be called after every outer cycle.
   * This observer will be called  \f$ N_P \times N_O \f$ times
   * @param o - shared pointer to an object inheriting Evaluator interface
   */
  void outer_cycle_observer(observers::ObserverInterface_SP o) { observe_every_outer_cycle.push_back(o); }

  /// Evaluate all outer cycle evaluators
  void call_outer_cycle_evaluators() { for (const auto &e : evaluate_every_outer_cycle)  e->evaluate(); }

  /// Evaluate all inner cycle evaluators
  void call_inner_cycle_evaluators() { for (const auto &e : evaluate_every_inner_cycle) e->evaluate(); }

  /// Call all outer cycle observers
  void call_outer_cycle_observers() { for (const auto &e : observe_every_outer_cycle) e->observe();  }

  /// Call all inner cycle observers
  void call_inner_cycle_observers() { for (const auto &e : observe_every_inner_cycle) e->observe(); }

protected:
  core::index4 n_outer_cycles;
  core::index4 n_inner_cycles;
  core::index4 n_cycle_size;
  std::vector<evaluators::Evaluator_SP> evaluate_every_inner_cycle;
  std::vector<evaluators::Evaluator_SP> evaluate_every_outer_cycle;
  std::vector<observers::ObserverInterface_SP> observe_every_inner_cycle;
  std::vector<observers::ObserverInterface_SP> observe_every_outer_cycle;
private:
  friend class ReplicaExchangeMC; // This is necessary so REMC can exchange also observers (thus observations are isothermal)
};

/// Define a shared pointer to SamplingProtocolBase as a new type
typedef std::shared_ptr<SamplingProtocolBase> SamplingProtocolBase_SP;

} // ~ sampling
} // ~ simulations

#endif
