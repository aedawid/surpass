#ifndef SIMULATIONS_GENERIC_SAMPLING_MetropolisAcceptanceCriterion_HH
#define SIMULATIONS_GENERIC_SAMPLING_MetropolisAcceptanceCriterion_HH

#include <random>

#include <core/real.hh>
#include <core/calc/statistics/Random.hh>

#include <simulations/sampling/AbstractAcceptanceCriterion.hh>

namespace simulations {
namespace sampling {

/** @brief Provides a Metropolis acceptance criterion to be used in Monte Carlo sampling.
 */
class MetropolisAcceptanceCriterion : public AbstractAcceptanceCriterion {
public:

  /** @brief Creates an acceptance criterion for a given temperature
   *
   * @param temperature - temperature defines the canonical distribution
   */
  MetropolisAcceptanceCriterion(const core::real temperature) : rando(0.0, 1.0) {
    this->temperature = temperature;
  }


  /// Returns the temperature used by this criterion
  core::real get_temperature() const { return temperature; }

  /// Sets the new value of the temperature used by this criterion
  void set_temperature(const core::real new_temperature) { this->temperature = new_temperature; }

  /** @brief Performs the Monte Carlo test
   *
   * @param old_energy - energy before the considered move
   * @param new_energy - energy after the considered move
   * @return true if the move should be accepted; false otherwise
   */
  inline bool test(const core::real old_energy, const core::real new_energy) {

    core::real delta_E = new_energy - old_energy;
    if (delta_E > 0) { if (rando(generator) > exp(-delta_E / temperature)) return false; }
    return true;
  }

private:
  core::calc::statistics::Random &generator = core::calc::statistics::Random::get();
  std::uniform_real_distribution<float> rando;
  core::real temperature = 0;
};

} // ~ movers
} // ~ simulations

#endif
