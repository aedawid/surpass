#ifndef SIMULATIONS_GENERIC_ObserveTopologyMatrix_HH
#define SIMULATIONS_GENERIC_ObserveTopologyMatrix_HH

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <unordered_map>

#include <utils/Logger.hh>

#include <simulations/evaluators/Evaluator.hh>
#include <simulations/observers/ToStreamObserver.hh>
#include <simulations/forcefields/surpass/SurpassHydrogenBond.hh>

namespace simulations {
namespace observers {

/** @brief Creates an observer which for each observation writes topology matrix of a surpass model
 *
 * At every <code>observe()</code> call this object will write topology matrix : all its elements in a single line
 */
template <typename C>
class ObserveTopologyMatrix : public virtual ToStreamObserver {
public:

  /** @brief Creates an observer that writes topology matrix elements into a given stream
   *
   * @param out - output stream where the data will be written
   */
  ObserveTopologyMatrix(std::shared_ptr<simulations::forcefields::surpass::SurpassHydrogenBond<C>> hb_energy,
    std::shared_ptr<std::ostream> out) : system_(hb_energy), logger("ObserveTopologyMatrix"), outstream(out), is_file_(false) {
    n_topologies_ = 0;
  }

  /** @brief Creates an observer that writes topology matrix elements into a given file
 *
 * @param out_fname - name of the output file
 */
  ObserveTopologyMatrix(std::shared_ptr<simulations::forcefields::surpass::SurpassHydrogenBond<C>> hb_energy, std::string out_fname)
    : ObserveTopologyMatrix(hb_energy,std::make_shared<std::ofstream>(out_fname)) { is_file_ = true; }

  /// Virtual destructor
  virtual ~ObserveTopologyMatrix() {}

  virtual bool observe();

  virtual void finalize();

  virtual std::shared_ptr<std::ostream> output_stream() { return outstream; };

  virtual void output_stream(std::shared_ptr<std::ostream> out) { outstream = out; };

  virtual core::index4 count_observe_calls() const { return cnt; }

private:
  std::unordered_map<std::string,core::index4> observed_topologies;
  std::vector<core::index4> topology_counts;
  core::index4 n_topologies_;
  std::shared_ptr<simulations::forcefields::surpass::SurpassHydrogenBond<C>> system_;
  utils::Logger logger;
  std::shared_ptr<std::ostream> outstream;
  bool is_file_;
  core::index4 cnt = 0;
};

template <typename C>
void ObserveTopologyMatrix<C>::finalize() {

  if(is_file_) {
    std::shared_ptr<std::ofstream> of = std::static_pointer_cast<std::ofstream>(outstream);
    of->close();
  } else outstream->flush();
}

template <typename C>
bool ObserveTopologyMatrix<C>::observe() {

  ++cnt;
  if (!trigger->operator()()) return false;

  std::stringstream topo_stream;
  for (core::index1 i = 0; i < system_->beta_topology_matrix().count_rows(); ++i)
    for (core::index1 j = 0; j < system_->beta_topology_matrix().count_columns(); ++j) topo_stream << int(
        system_->beta_topology_matrix()(i, j));

  std::string topo = topo_stream.str();
  if (observed_topologies.find(topo) == observed_topologies.end()) {
    observed_topologies[topo] = n_topologies_;
    topology_counts.push_back(1);
    ++n_topologies_;
  }
  else
    topology_counts[observed_topologies[topo]]++;

  (*outstream) << std::setw(6) << cnt << " " << topo << " " << observed_topologies[topo] << "\n";
  outstream->flush();

  return true;
}

} // ~ simulations
} // ~ observers
#endif
