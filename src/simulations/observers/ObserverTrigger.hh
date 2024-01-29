/** @file ObserverTrigger.hh
 *  @brief Provides ObserverTrigger
 */
#ifndef SIMULATIONS_OBSERVERS_ObserverTrigger_HH
#define SIMULATIONS_OBSERVERS_ObserverTrigger_HH

namespace simulations {
namespace observers {

/** @brief Object used to tell an observer that it actually should take an observation.
 *
 * Any observer calls a ObserverTrigger instance; when the trigger returns true, the observer actually takes
 * its observations. Otherwise it pass by silently.
 */
class ObserverTrigger {
public:

  /// This base class implementation always returns true i.e. every observation is actually recorded
  virtual bool operator()() { return true; }

  /// Says how many times a trigger has been called
  core::index4 count_observe_calls() const { return count_observe_calls_; }

protected:
  core::index4 count_observe_calls_;
};

typedef std::shared_ptr<ObserverTrigger> ObserverTrigger_SP;

}
}

#endif
