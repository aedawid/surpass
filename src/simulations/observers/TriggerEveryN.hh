/** @file TriggerEveryN.hh
 *  @brief Provides TriggerEveryN object that triggers observations
 */
#ifndef SIMULATIONS_OBSERVERS_TriggerEveryN_HH
#define SIMULATIONS_OBSERVERS_TriggerEveryN_HH

#include <memory>

#include <simulations/observers/ObserverTrigger.hh>

namespace simulations {
namespace observers {

/** @brief Triggers every N-th observation.
 */
class TriggerEveryN : public ObserverTrigger {
public:

  /** @brief Creates a trigger that takes only every n-th observation.
   * All other will be neglected.
   *
   * @param n - how often observations should actually be observed. By default  <code>n = 1</code>
   * which means that everything is observed;  <code>n = 2</code> meas that every second call is effective
   */
  TriggerEveryN(core::index2 n = 1) : observe_every_n_(n) {}

  /** @brief Accepts every N-th observation.
   *
   * This trigger returns true every N calls; every other time returns false. One should not call this operator
   * manually because this will affect the internal counter for triggering observations;
   *
   * @return true if an observation should be actually taken
   */
  virtual bool operator()() { return (((++count_observe_calls_) % observe_every_n_) == 0); }

  /** @brief Sets the new frequency of observations
   *
   * @param n - how often observations should actually be observed. By default  <code>n = 1</code>
   * which means that everything is observed;  <code>n = 2</code> meas that every second call is effective
   */
  void observe_every_n(core::index2 n) { observe_every_n_ = n; }

  /** @brief Check how often observations are being taken.
   *
   * @return how often observations are taken? 1 means everything is observed
   */
  core::index2 observe_every_n() const { return observe_every_n_; }

private:
  core::index2 observe_every_n_;
};

}
}

#endif
