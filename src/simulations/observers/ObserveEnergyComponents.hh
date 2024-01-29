#ifndef SIMULATIONS_OBSERVERS_ObserveEnergyComponents_HH
#define SIMULATIONS_OBSERVERS_ObserveEnergyComponents_HH

#include <vector>
#include <memory>
#include <iostream>

#include <utils/Logger.hh>

#include <simulations/evaluators/Evaluator.hh>
#include <simulations/observers/ToStreamObserver.hh>
#include <simulations/forcefields/TotalEnergyByResidue.hh>

namespace simulations {
namespace observers {

using simulations::evaluators::Evaluator_SP;

/** @brief Observes all energy components of a given TotalEnergyByResidue instance.
 *
 * Each observation makes a row of energy values in the output stream.
 * @tparam E - the type of energy components, e.g. <code>ByResidueEnergy</code> or <code>CalculateEnergyBase</code>
 */
template <typename E>
class ObserveEnergyComponents : public virtual ToStreamObserver {
public:

  /** @brief Creates an observer that evaluates and writes energy components into a given stream as a nice table.
   *
   * @param total_energy - total energy of a system; energy components contained in that object will be evaluated and printed
   * in a single row at each <code>observe()</code> call
   * @param out - output stream where the table will be written
   */
  ObserveEnergyComponents(forcefields::TotalEnergy<E> & total_energy, std::shared_ptr<std::ostream> out);

  /** @brief Creates an observer that evaluates and writes energy components as a nice table on the standard output.
   * @param total_energy - total energy of a system; energy components contained in that object will be evaluated and printed
   * in a single row at each <code>observe()</code> call
   */
  ObserveEnergyComponents(forcefields::TotalEnergy<E> & total_energy) :
    ObserveEnergyComponents(total_energy, std::shared_ptr<std::ostream>(&std::cout, [](void *) {})) {}

  /** @brief Creates an observer that evaluates and writes energy components as a nice table into a file
   * @param total_energy - total energy of a system; energy components contained in that object will be evaluated and printed
   * in a single row at each <code>observe()</code> call
   * @param fname - name of the output file
   */
  ObserveEnergyComponents(forcefields::TotalEnergy<E> & total_energy, const std::string & fname) :
    ObserveEnergyComponents(total_energy,std::make_shared<std::ofstream>(fname)) { is_file_ = true; }

  /// Virtual destructor
  virtual ~ObserveEnergyComponents() {}

  virtual bool observe();

  virtual void finalize();

  /// Returns the header of an energy table as a string
  std::string header_string() const;

  /** @brief Writes the energy table header line to the stream.
   *
   * @param prefix - a string to be printed just ahead of the header line. By default a "#     " string is used
   * so the header line is interpreted as a comment e.g. by gnuplot. Also, <code>observe()</code> method prints
   * observation counter in the first column and the '#' character works as an id for this column
   */
  void observe_header(const std::string & prefix = "#      ");

  virtual std::shared_ptr<std::ostream> output_stream() { return outstream; };

  virtual void output_stream(std::shared_ptr<std::ostream> out) { outstream = out; };

private:
  utils::Logger logger;
  forcefields::TotalEnergy<E> & total_energy_;
  std::shared_ptr<std::ostream> outstream;
  bool is_file_;
  core::index4 cnt = 0;
};

} // ~ simulations
} // ~ observers
#endif
