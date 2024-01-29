#ifndef SIMULATIONS_OBSERVERS_AbstractPdbObserver_HH
#define SIMULATIONS_OBSERVERS_AbstractPdbObserver_HH

#include<iostream>
#include <fstream>

#include <core/data/basic/Vec3.hh>
#include <simulations/systems/CartesianAtomsSimple.hh>
#include <simulations/observers/ToStreamObserver.hh>

namespace simulations {
namespace observers {
namespace cartesian {

/** @brief Provides an observer that converts a conformation to a PDB format.
 *
 * This is an abstract class since the sink of the PDB - formatted data is still to be defined.
 * Currently there are two derived classes that make use of this PDB strings: PdbTrajectory and PymolObserver
 *
 * The output PDB format is based on an existing Structure object; the following PDB format data will be borrowed from the template structure:
 *   - atom names
 *   - residue names
 *   - atom and residue numbering
 *   - chain ID
 * @see PdbObserver<C>  PymolObserver<C>
 */
template<typename C>
class AbstractPdbObserver : public ToStreamObserver {
public:

  const simulations::systems::CartesianAtomsSimple<C> &observed_object;

  /** @brief Creates an observer that converts a conformation to a PDB format.
   *
   * @param observed_object - system whose coordinates will be stored in the file
   * @param pdb_format_source - biomolecular structure that corresponds to the system. PDB format data:
   *   - atom names
   *   - residue names
   *   - atom and residue numbering
   *   - chain ID
   * will be extracted from this structure. Coordinates will be taken from <code>observed_object</code>
   * @see PdbTrajectory<C>  PymolObserver<C>
   */
  AbstractPdbObserver(const systems::CartesianAtomsSimple<C> &observed_object,
      const core::data::structural::Structure &pdb_format_source) : observed_object(observed_object) {

    std::string pdb_atom_fmter = "ATOM  %5d %s %s %c%4d    %%8.3f%%8.3f%%8.3f  1.00 99.99\n";
    for (auto atom_it = pdb_format_source.first_const_atom();
         atom_it != pdb_format_source.last_const_atom(); ++atom_it) {
      const auto a = **atom_it;
      const auto r = *(a.owner());
      format_lines.push_back(
        utils::string_format(pdb_atom_fmter, a.id(), a.atom_name().c_str(), r.residue_type().code3.c_str(),
          r.owner()->id(), r.id()));
    }
  }

protected:
  std::vector<std::string> format_lines; ///< formatting string for every atom in the structure
};

}
}
}
#endif
