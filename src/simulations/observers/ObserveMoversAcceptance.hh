#ifndef SIMULATIONS_GENERIC_ObserveMoversAcceptance_HH
#define SIMULATIONS_GENERIC_ObserveMoversAcceptance_HH

#include <vector>
#include <memory>
#include <iostream>

#include <utils/Logger.hh>

#include <simulations/evaluators/Evaluator.hh>
#include <simulations/observers/ToStreamObserver.hh>
#include <simulations/forcefields/TotalEnergyByResidue.hh>
#include <simulations/movers/MoversSet.hh>

namespace simulations {
namespace observers {

using simulations::evaluators::Evaluator_SP;

/** @brief Observes acceptance rate for every mover contained in the given MoversSet.
 *
 * Each observation makes a row of acceptance rate values in the output stream.
 */
class ObserveMoversAcceptance : public virtual ToStreamObserver {
public:

  /** @brief Creates an observer that evaluates acceptance rate for each mover and writes into a given stream as a nice table.
   *
   * @param movers_set - MoversSet object
   * @param out - output stream where the table will be written
   */
  ObserveMoversAcceptance(simulations::movers::MoversSet &movers_set, std::shared_ptr<std::ostream> out) :
    logger("ObserveMoversAcceptance"), ms_(movers_set), outstream(out), is_file_(false) {}

  /** @brief Creates an observer that evaluates acceptance rate for each mover and writes into a standard output
   * @param movers_set - MoversSet object
   */
  ObserveMoversAcceptance(simulations::movers::MoversSet &movers_set) :
    ObserveMoversAcceptance(movers_set, std::shared_ptr<std::ostream>(&std::cout, [](void *) {})) {}

  /** @brief Creates an observer that evaluates and writes energy components as a file
   * @param movers_set - MoversSet object
   * @param fname - name of the output file
   */
  ObserveMoversAcceptance(simulations::movers::MoversSet &movers_set, const std::string &fname) :
    logger("ObserveMoversAcceptance"), ms_(movers_set), outstream(std::make_shared<std::ofstream>(fname)),
    is_file_(true) {
  }

  /// Virtual destructor
  virtual ~ObserveMoversAcceptance() {}

  virtual bool observe();

  virtual void finalize();

  /// Returns the header of an output table as a string
  std::string header_string() const { return ms_.header_string(); }

  /** @brief Writes the movers acceptance ratios table header line to the stream.
   *
   * @param prefix - a string to be printed just ahead of the header line. By default a '#' character is used
   * so the header line is interpreted as a comment e.g. by gnuplot
   */
  void observe_header(const std::string &prefix = "#") { (*outstream) << prefix << ms_.header_string() << "\n"; }

  virtual std::shared_ptr<std::ostream> output_stream() { return outstream; };

  virtual void output_stream(std::shared_ptr<std::ostream> out) { outstream = out; };

  virtual core::index4 count_observe_calls() const { return cnt; }

private:
  utils::Logger logger;
  simulations::movers::MoversSet &ms_;
  std::shared_ptr<std::ostream> outstream;
  bool is_file_;
  core::index4 cnt = 0;
};

/// Declare a type of a shared pointer to ObserveMoversAcceptance
typedef std::shared_ptr<ObserveMoversAcceptance> ObserveMoversAcceptance_SP;

} // ~ simulations
} // ~ observers
#endif
