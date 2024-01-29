#ifndef SIMULATIONS_GENERIC_ObserveEvaluators_HH
#define SIMULATIONS_GENERIC_ObserveEvaluators_HH

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>

#include <utils/Logger.hh>

#include <simulations/evaluators/Evaluator.hh>
#include <simulations/observers/ToStreamObserver.hh>

namespace simulations {
namespace observers {

using simulations::evaluators::Evaluator_SP;

/** @brief Creates an observer that writes evaluated values into a given stream as a nice table.
 *
 * At every <code>observe()</code> call this object will call  <code>evaluate()</code> method from each of Evaluator
 * instances gathered in this observer. The evaluated values will be printed as a single row of a table.
 */
class ObserveEvaluators : public virtual ToStreamObserver {
public:

  /** @brief Creates an observer that writes evaluated values into a given stream
   *
   * @param out - output stream where the table will be written
   */
  ObserveEvaluators(std::shared_ptr<std::ostream> out) : logger("ObserveEvaluators"), outstream(out), is_file_(false) {}

  /** @brief Creates an observer that writes evaluated values on the standard output.
   */
  ObserveEvaluators() : ObserveEvaluators(std::shared_ptr<std::ostream>(&std::cout, [](void *) { })) {}

  /** @brief Creates an observer that writes evaluated values into a file
   */
  ObserveEvaluators(const std::string & fname) : ObserveEvaluators(std::make_shared<std::ofstream>(fname)) { is_file_ = true; }

  /// Virtual destructor
  virtual ~ObserveEvaluators() {}

  virtual bool observe();

  virtual void finalize();

  void add_evaluator(evaluators::Evaluator_SP evaluator);

  /// Returns the header of an output table as a string
  std::string header_string() const;

  /** @brief Writes the  table header line to the stream.
   *
   * @param prefix - a string to be printed just ahead of the header line. By default a "#     " string is used
   * so the header line is interpreted as a comment e.g. by gnuplot. Also, <code>observe()</code> method prints
   * observtion counter in the first column and the '#' character works as an id for this column
   */
  void observe_header(const std::string & prefix = "#     ") { (*outstream) << prefix << header_string() << "\n"; }

  std::vector<Evaluator_SP>::iterator begin() { return evaluators.begin(); }

  std::vector<Evaluator_SP>::iterator end() { return evaluators.end(); }

  virtual std::shared_ptr<std::ostream> output_stream() { return outstream; };

  virtual void output_stream(std::shared_ptr<std::ostream> out) { outstream = out; };

  virtual core::index4 count_observe_calls() const { return cnt; }

private:
  std::vector<Evaluator_SP> evaluators;
  utils::Logger logger;
  std::vector<core::index2> sw;
  std::shared_ptr<std::ostream> outstream;
  bool is_file_;
  core::index4 cnt = 0;
};

/// Declare a type of a shared pointer to ObserveEvaluators
typedef std::shared_ptr<ObserveEvaluators> ObserveEvaluators_SP;

} // ~ simulations
} // ~ observers
#endif
