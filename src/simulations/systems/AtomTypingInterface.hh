#ifndef SIMULATIONS_SYSTEMS_AtomTypingInterface_HH
#define SIMULATIONS_SYSTEMS_AtomTypingInterface_HH

#include <memory>

#include <core/index.hh>
#include <core/data/basic/Vec3.fwd.hh>
#include <core/data/structural/PdbAtom.fwd.hh>
#include <simulations/systems/AtomTypingVariants.hh>

namespace simulations {
namespace systems {

/** @brief AtomTyping converts atoms (identified by their names as they are stored in a PDB file) into internal numbering used by force fields.
 *
 * This interface declares ability to convert from atom name into an internal index and back. Each atom
 * has two names assigned:
 *   - one according to PDB naming conventions, e.g. alpha carbon in proteins in " CA ", and
 *   - name that is specific for a force field.
 */
class AtomTypingInterface  {
public:

  /// empty destructor
  virtual ~AtomTypingInterface() {}

  /** @brief Returns an internal index for an atom identified by its PDB name, its residue name and a requested variant.
   * @param atom_name - name of the requested atom, as written in the PDB file used to start a simulation
   * @param residue_name - name of a residue
   * @param variant - variand of the residue; typically AtomTypingVariants::STANDARD should be used
   *    unless its a special case for this residue, e.g. it is the N-terminal or C-terminal one
   */
  virtual core::index2 atom_type(const std::string & atom_name, const std::string & residue_name,
                                 const AtomTypingVariants variant) const = 0;

  /// Returns an internal index for an atom identified by its PDB name and its residue name
  virtual core::index2 atom_type(const core::data::structural::PdbAtom & a) const = 0;

  /** @brief Returns an internal atom name for a given internal atom type.
   *  @param atom_id - index assigned to a given atom type according to this atom typing
   *  @return name for a given atom as it is assigned by this atom typing
   */
  virtual const std::string & atom_internal_name(const core::index2 atom_id) const = 0;

  /** @brief Returns an internal atom name based on <code>atom_id</code> field that is stored in Vec3 type
   *  @param atom - atom data
   *  @return name for a given atom as it is assigned by this atom typing
   */
  virtual const std::string & atom_internal_name(const core::data::basic::Vec3 & atom) const = 0;

  /// @brief Count how many different atom types are defined by this typing
  virtual core::index2 count_atom_types() const = 0;

  /** @brief Returns atom type index for a given internal atom name.
   * @param atom_internal_name - name of this atom type, as it is defined in this atom typing. This could not be
   * the PDB atom name, because in general more than one atom type may be assigned to a given PDB atom name
   * (e.g. CG1 from VAL may be a different atom than  CG1 in TRP)
   */
  virtual core::index2 atom_type(const std::string & atom_internal_name) const = 0;
};

/// Define a shared pointer type
typedef std::shared_ptr<AtomTypingInterface> AtomTypingInterface_SP;

}
}
#endif
