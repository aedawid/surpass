#ifndef SIMULATIONS_GENERIC_EVALUATORS_EchoEvaluator_HH
#define SIMULATIONS_GENERIC_EVALUATORS_EchoEvaluator_HH

#include <core/real.hh>

#include <simulations/evaluators/Evaluator.hh>

namespace simulations {
namespace evaluators {

/** @brief EchoEvaluator just returns a value stored under a reference.
 *
 * At construction time, a reference to a value of a generic type T must be provided. Later on, at every <code>evaluate()</code> call
 * the value of this reference is returned.
 *
 * The following example declares two  <code>EchoEvaluator</code> instances. It also declares a
 * <code>CallEvaluator</code> instance
 * \include ex_SimulatedAnnealing.cc
 */
template <typename T>
class EchoEvaluator : public Evaluator {
public:

  /** @brief Creates an evaluator that just returns a value stored in a given reference.
   * @param var - value of this variable will be returned at every <code>evaluate()</code> call
   * @param name - a name assigned to the observed variable
   */
  EchoEvaluator(const T & var,const std::string & name, core::index2 min_width = 8, core::index1 precision = 2) :
    observed(var), name_(name), min_width_(min_width), precision_(precision) { }

  virtual core::real evaluate() { return (core::real) observed; }

  virtual const std::string & name() const { return name_; }

  virtual core::index1 precision() const { return precision_; }

  void precision(const core::index1 precision) { precision_ = precision; }

  virtual core::index2 min_width() const { return min_width_; }

  virtual ~EchoEvaluator() {}

private:
  const T & observed;
  const std::string name_;
  core::index2 min_width_;
  core::index1 precision_;
};

} // ~ evaluators
} // ~ simulations
#endif
