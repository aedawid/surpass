#ifndef SIMULATIONS_SYSTEMS_SingleAtomType_HH
#define SIMULATIONS_SYSTEMS_SingleAtomType_HH

#include <core/index.hh>

#include <simulations/systems/AtomTypingInterface.hh>

namespace simulations {
namespace systems {

/** @brief There is only one atom type in this AtomTyping.
 *
 * By default every atom is identified by this typing as an alpha carbon (named <code>" CA "</code>),
 * but the name may ba changed. Its code is always 0
 */
class SingleAtomType : public AtomTypingInterface {
public:

  SingleAtomType(const std::string & atom_type_name = " CA ") : ca_name(atom_type_name) {}

  /// @brief Virtual destructor to satisfy a compiler
  virtual ~SingleAtomType() {}

  /// @brief Always returns 0, because there is only one atom type according to this AtomTyping
  virtual core::index2 atom_type(const std::string & atom_internal_name) const { return 0; }

  /// @brief Always returns 0, because there is only one atom type according to this AtomTyping
  virtual core::index2 atom_type(const std::string & atom_name, const std::string & residue_name) const { return 0; }

  /// @brief Always returns 0, no matter what the atom is <code>a</code>
  virtual core::index2 atom_type(const core::data::structural::PdbAtom & a) const { return 0; }

  /// @brief Always returns 0, because there is only one atom type according to this AtomTyping
  virtual core::index2 atom_type(const std::string & atom_name, const std::string & residue_name,const AtomTypingVariants variant) const { return 0; }

  /// @brief Always returns the atom type used to construct this object
  virtual const std::string & atom_internal_name(const core::index2 atom_id) const { return ca_name; }

  /// @brief Always returns the atom type used to construct this object
  virtual const std::string & atom_internal_name(const core::data::basic::Vec3 & atom) const { return ca_name; }

  /// @brief Always returns 0the atom type used to construct this object
  virtual core::index2 atom_id(const std::string & atom_internal_name) const { return 0; }

  /// @brief Just return 1
  virtual core::index2 count_atom_types() const { return 1; }

private:
  const std::string ca_name;
};

}
}
#endif
