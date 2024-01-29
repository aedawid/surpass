/** @file ObserverInterface.hh
 *  @brief Provides ObserverInterface and ObserverInterface_SP types
 */
#ifndef SIMULATIONS_OBSERVERS_ObserverInterface_HH
#define SIMULATIONS_OBSERVERS_ObserverInterface_HH

#include <memory>

#include <simulations/systems/CartesianAtomsSimple.hh>
#include <simulations/observers/ObserverTrigger.hh>

namespace simulations {
namespace observers {

/** @brief Observer interface provides <code>bool observe()</code> method.
 *
 * The purpose of an observer is to record a value or values during sampling. The observed value may be written to
 * a stream, stored in a histogram, etc.
 */
class ObserverInterface {
public:

  /// Create an observer which records every observations
  ObserverInterface(ObserverTrigger_SP trigger = std::make_shared<ObserverTrigger>()) : trigger(trigger) {}

  /// This method is called to take observations (to be implemented by a derived class)
  virtual bool observe() = 0;

  /// This method will be called before the program shuts down, e.g. to close open files
  virtual void finalize() = 0;

  /** @brief Replaces the trigger for this observer with a new one.
   *
   * This call will change the way how often observations are taken
   *
   * @param new_trigger - instance of a new trigger object
   */
  void set_trigger(ObserverTrigger_SP new_trigger) { trigger = new_trigger; }

protected:
  ObserverTrigger_SP trigger;
};

/// Declares a shared pointer to ObserverInterface type
typedef std::shared_ptr<ObserverInterface> ObserverInterface_SP;
}
}

#endif
