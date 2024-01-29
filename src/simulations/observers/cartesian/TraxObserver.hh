#ifndef SIMULATIONS_OBSERVERS_TraxObserver_HH
#define SIMULATIONS_OBSERVERS_TraxObserver_HH

#include <iostream>
#include <fstream>

#include <simulations/observers/ObserveEvaluators.hh>
#include <simulations/systems/ResidueChain.hh>
#include <simulations/forcefields/TotalEnergyByResidue.hh>

namespace simulations {
namespace observers {
namespace cartesian {

/** @brief Writes a trajectory for an observed system in TRAX format.
 */
template<typename C>
class TraxObserver : public ObserveEvaluators {
public:

  const simulations::systems::ResidueChain<C> &observed_object;

  /** @brief Creates an observer that writes a conformation to a file in TRAX format
   * @param observed_object - system whose coordinates will be stored in the file
   * @param max_atom_in_residue - maximum number of atoms a residue may have
   * @param out_fname - name of the output file
   */
  TraxObserver(const systems::ResidueChain<C> &observed_object, const core::index1 max_atom_in_residue, const std::string & out_fname);

  /// Append a new frame to the file
  bool observe();

  void finalize() { (std::static_pointer_cast<std::ofstream>(outstream))->close(); }

  /// The trajectory file will hold also energy data
  void observe_energy(std::shared_ptr<forcefields::TotalEnergyByResidue> en) { energy_ = en; }

  /// Returns the header of an output table as a string
  std::string header_string() const;

  void observe_header(const std::string & prefix = "#") { (*outstream) << prefix << header_string() << "\n"; }

private:
  core::index4 n_frames; ///< Counter of frames written so far
  core::index4 cnt = 0; ///< counts how many times the observe() method was called
  core::index1 atoms_per_line;
  std::string out_fname;
  std::shared_ptr<std::ostream> outstream;
  std::shared_ptr<forcefields::TotalEnergyByResidue> energy_;
  std::shared_ptr<std::stringstream> evaluator_stream;
};

}
}
}
#endif
