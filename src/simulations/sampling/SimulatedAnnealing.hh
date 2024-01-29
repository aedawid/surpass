#ifndef SIMULATIONS_GENERIC_SAMPLING_SimulatedAnnealing_HH
#define SIMULATIONS_GENERIC_SAMPLING_SimulatedAnnealing_HH

#include <core/real.hh>

#include <utils/Logger.hh>

#include <simulations/movers/MoversSet.hh>
#include <simulations/sampling/IsothermalMC.hh>

namespace simulations {
namespace sampling {

/** @brief Simulated annealing protocol.
 *
 * This protocol executes a Monte Carlo sweep defined by the given MoversSet instance
 * \f$ N_I \times N_O \f$ times for each of  \f$ N_T \f$ temperature where  \f$ N_I \f$ and  \f$ N_O \f$
 * is the number of inner and outer cycles, respectively. Thus the total number of Monte Carlo sweeps in this protocol is
 * \f$ N_I \times N_O \times N_T \f$.
 *
 * \include ex_SimulatedAnnealing.cc
 */
class SimulatedAnnealing : public IsothermalMC {
public:

  SimulatedAnnealing(movers::MoversSet_SP ms, const std::vector<core::real> &temperatures)
    : IsothermalMC(ms), logger("SimulatedAnnealing"), temperatures(temperatures) { }

  /// Virtual destructor
  ~SimulatedAnnealing() {}

  /** Brief Runs the sampling protocol.
   * For each temperature, the method makes  \f$ N_O \f$ = <code>outer_cycles()</code> of
   * \f$ N_I \f$ = <code>inner_cycles()</code> of Monte Carlo steps.
   * The size of each MC sweep is defined within the MoversSet instance given to constructor.
   */
  void run();

  /// Returns current simulation temperature
  core::real temperature() const { return temperature_; }

private:
  utils::Logger logger;
  const std::vector<core::real> temperatures;
};

} // ~ sampling
} // ~ simulations

#endif

/**
 * \example ex_SimulatedAnnealing.cc
 */