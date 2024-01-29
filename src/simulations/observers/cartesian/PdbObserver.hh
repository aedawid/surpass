#ifndef SIMULATIONS_OBSERVERS_PdbObserver_HH
#define SIMULATIONS_OBSERVERS_PdbObserver_HH

#include<iostream>
#include <fstream>

#include <core/data/basic/Vec3.hh>
#include <simulations/systems/CartesianAtomsSimple.hh>
#include <simulations/observers/cartesian/AbstractPdbObserver.hh>

namespace simulations {
namespace observers {
namespace cartesian {

/** @brief Writes a PDB trajectory for an observed system.
 *
 */
template<typename C>
class PdbObserver : public AbstractPdbObserver<C> {
public:

  /** @brief Creates an observer that writes a conformation to a file in PDB format
   * @param observed_object - system whose coordinates will be stored in the file
   * @param pdb_format_source - biomolecular structure that corresponds to the system. PDB format data:
   *   - atom names
   *   - residue names
   *   - atom and residue numbering
   *   - chain ID
   * will be extracted from this structure. Coordinates will be taken from <code>observed_object</code>
   * @param out_fname - name of the output file
   */
  PdbObserver(const systems::CartesianAtomsSimple<C> &observed_object,
              const core::data::structural::Structure &pdb_format_source, const std::string & out_fname) :
    AbstractPdbObserver<C>(observed_object, pdb_format_source), out_fname(out_fname) {
    outstream =   std::make_shared<std::ofstream>(out_fname);
  }

  /// Append PDB - formatted data to the file
  bool observe();

  /** @brief Observe another set of coordinates to this stream.
   *
   * The given system should be a conformation of the same system! Otherwise the PDB format will be totally wrong.
   * This method has been introduced to save replica conformations to a single file.
   *
   * @param system - a conformation corresponding to the observed system
   * @return true
   */
  bool observe(const simulations::systems::CartesianAtomsSimple<C> & system);

  void finalize() { (std::static_pointer_cast<std::ofstream>(outstream))->close(); }

  virtual std::shared_ptr<std::ostream> output_stream() { return outstream; };

  virtual void output_stream(std::shared_ptr<std::ostream> out) { outstream = out; };

  virtual core::index4 count_observe_calls() const { return cnt; }

private:
  std::string out_fname;
  core::index4  cnt = 0;
  std::shared_ptr<std::ostream> outstream;
};

template<typename C>
bool PdbObserver<C>::observe() {

  ++cnt;
  if(!ObserverInterface::trigger->operator()()) return false;

  AbstractPdbObserver<C>::observed_object.write_pdb(*outstream, AbstractPdbObserver<C>::format_lines, cnt);
  outstream->flush();

  return true;
}

template<typename C>
bool PdbObserver<C>::observe(const simulations::systems::CartesianAtomsSimple<C> & system) {

  ++cnt;
  system.write_pdb(*outstream, AbstractPdbObserver<C>::format_lines, cnt);
  outstream->flush();

  return true;
}

/** @brief A simple utility method that create a PdbObserver just to write a single conformation into a PDB file
 * @param observed_object - system whose coordinates will be stored in the file
 * @param pdb_format_source - biomolecular structure that corresponds to the system. PDB format data:
 * @param out_fname - name of the output file
 */
template<typename C>
void write_pdb_conformation(const systems::CartesianAtomsSimple<C> &observed_object,
                            const core::data::structural::Structure &pdb_format_source, const std::string out_fname) {

  PdbObserver<C> o(observed_object, pdb_format_source, out_fname);
  o.observe();
}

}
}
}
#endif
