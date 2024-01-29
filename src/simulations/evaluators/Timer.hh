#ifndef SIMULATIONS_GENERIC_EVALUATORS_Timer_HH
#define SIMULATIONS_GENERIC_EVALUATORS_Timer_HH

#include <memory>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <core/real.hh>

#include <simulations/evaluators/Evaluator.hh>

namespace simulations {
namespace evaluators {

using namespace std::chrono;

/** @brief Evaluates the time elapsed from the creation of this object.
 *
 * At every <code>evaluate()</code> method call the time in seconds is returned.
 */
class Timer : public Evaluator {
public:

  /// Creates a new timer that starts running at the time when this constructor is called
  Timer() {  start = std::chrono::high_resolution_clock::now(); }

  /** @brief Returns the time elapsed since the creation of this timer
   *
   */
  virtual core::real evaluate() {

    std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
    duration<double> time_span = std::chrono::duration_cast<duration<double>>(now - start);
    return time_span.count();
  };

  virtual const std::string & name() const { return name_; }

  virtual core::index1 precision() const { return 2; }

  virtual core::index2 min_width() const { return 8; }

  virtual ~Timer() {}

private:
  std::chrono::high_resolution_clock::time_point start;
  static const std::string name_;
};

} // ~ simulations
} // ~ evaluators

#endif
