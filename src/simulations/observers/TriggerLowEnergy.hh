/** @file TriggerLowEnergy.hh
 *  @brief Provides TriggerLowEnergy object that triggers observations when energy is low enough
 */
#ifndef SIMULATIONS_OBSERVERS_TriggerLowEnergy_HH
#define SIMULATIONS_OBSERVERS_TriggerLowEnergy_HH

#include <memory>
#include <utils/Logger.hh>

#include <simulations/observers/ObserverTrigger.hh>
#include <simulations/forcefields/TotalEnergy.hh>

namespace simulations {
namespace observers {

/** @brief Triggers observation when energy is low enough.
 */
class TriggerLowEnergy : public ObserverTrigger {
public:

  /** @brief Creates a trigger that observes only low-energy events
   * All other will be neglected.
   *
   * @param n - energy function used to guide this trigger. The energy function is evaluated at every
   * <code>observe()</code> call
   */
  TriggerLowEnergy(simulations::forcefields::CalculateEnergyBase & energy, const core::real low_energy_value = 0,
    const core::real fraction = 0.1) : energy_(energy),
    low_energy_value_(low_energy_value), fraction_(fraction), logs("TriggerLowEnergy") {}

  /** @brief Accepts only low-energy observations
   *
   * @return true if the current energy value is low enough
   */
  virtual bool operator()() {

    ++count_observe_calls_;
    core::real en = energy_.calculate();
    core::real cutoff = (low_energy_value_ < 0) ? (1.0 - fraction_) * low_energy_value_ : (1.0 + fraction_) *
                                                                                          low_energy_value_;
    if (en < cutoff) {
      if (en < low_energy_value_) {
        low_energy_value_ = en;
        logs << utils::LogLevel::FINE << "trigger min-energy set to " << en << "\n";
      }
      return true;
    }
    return false;
  }

private:
  simulations::forcefields::CalculateEnergyBase & energy_;
  core::real low_energy_value_;
  core::real fraction_;
  utils::Logger logs;
};

}
}

#endif
