#ifndef SIMULATIONS_GENERIC_SAMPLING_IsothermalMC_HH
#define SIMULATIONS_GENERIC_SAMPLING_IsothermalMC_HH

#include <core/real.hh>

#include <simulations/movers/MoversSet.hh>
#include <simulations/sampling/SamplingProtocolBase.hh>

namespace simulations {
namespace sampling {

/** @brief Simple isothermal MC protocol.
 *
 * This protocol executes a Monte Carlo sweep defined by the given MoversSet instance
 * \f$ N_I \times N_O \f$ times.
 */
class IsothermalMC : public SamplingProtocolBase {
public:

  /** @brief Creates a new isothermal sampler.
   *
   * @param ms - set of movers used for sampling
   * @param temperature - temperature of the isothermal run
   */
  IsothermalMC(movers::MoversSet_SP & ms, const core::real temperature = 1.0) : movers(ms) { temperature_ = temperature; }

  /// Virtual destructor
  ~IsothermalMC() {}

  /** Brief Runs the sampling protocol.
   * At every call the method makes  \f$ N_O \f$ = <code>outer_cycles()</code> of
   * \f$ N_I \f$ = <code>inner_cycles()</code> of Monte Carlo steps.
   * The size of each MC sweep is defined within the MoversSet instance given to constructor.
   */
  void run();

  /// Sets the new temperature value and calls <code>run()</code>
  void run(const core::real new_temperature) { temperature_ = new_temperature; run(); }

  /// Returns current simulation temperature
  core::real temperature() const { return temperature_; }

  /// Sets the new value of the simulation temperature
  void temperature(const core::real new_temperature) { temperature_ = new_temperature; }

protected:
  movers::MoversSet_SP movers; ///< Movers to be called to sample
  core::real temperature_ = 0; ///< Current temperature
};

typedef std::shared_ptr<IsothermalMC> IsothermalMC_SP;

} // ~ sampling
} // ~ simulations

#endif
