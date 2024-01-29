#ifndef SIMULATIONS_SYSTEMS_AtomTypingBase_HH
#define SIMULATIONS_SYSTEMS_AtomTypingBase_HH

#include <unordered_map>

#include <core/index.hh>

#include <core/data/basic/Vec3.hh>
#include <core/data/structural/Structure.fwd.hh>
#include <core/data/structural/PdbAtom.hh>

#include <simulations/systems/AtomTypingInterface.hh>
#include <simulations/systems/AtomTypingVariants.hh>

#include <utils/Logger.hh>

namespace simulations {
namespace systems {

/** @brief Binds PDB-style atom names to internal atom names used by a given force field.
 */
struct AtomNamePair {
  std::string pdb_atom_name; ///< Atom name as seen in a PDB file
  std::string internal_name; ///< Atom name according to a given force field

  AtomNamePair() : pdb_atom_name(""), internal_name("") {}

  /// Constructor makes a new pair
  AtomNamePair(const std::string & name,const std::string & ename) : pdb_atom_name(name), internal_name(ename) {}

  /// Compare a pair to another one
  bool operator<(const AtomNamePair & p) const {

    if (pdb_atom_name.compare(p.pdb_atom_name) < 0) return true;
    else if (pdb_atom_name.compare(p.pdb_atom_name) > 0) return false;
    return (internal_name.compare(p.internal_name) < 0);
  }

};

/** @brief Atom typing based on PDB file
 *
 * Here each unique atom name found in a given PDB structure gets its own ID. There is no particular
 * order of the IDs enforced. Internal (typed) atom names are set to element names
 *
 * \include ex_UniquePdbTyping.cc
 */
class AtomTypingBase : public AtomTypingInterface {
public:

  /** @brief Creates atom typing based on a given biomolecular structure.
   *
   * Resulting atom typing will assing a distinct type to each unique atom name. Atoms of the same name comming from
   * different residue types will share the same atom type, e.g. " CA " atom from GLY and from ALA will be of the same type.
   * Atom type indexes are assigned in the order atoms appear in the PDB file
   * @param structure - a structure
   */
  AtomTypingBase(core::data::structural::Structure & structure) { compile_types(structure); }

  /** @brief Create new atom typing from a given atom name pairs.
   *
   * This constructor supports only AtomTypingVariants::STANDARD atom types
   * @param atom_types - a vector of AtomNamePair objects that define the mapping
   */
  AtomTypingBase(const std::vector<AtomNamePair> & atom_types);

  /** @brief Create new atom typing from a given atom names.
   *
   * This constructor supports only AtomTypingVariants::STANDARD atom types.
   * @param atom_types - a vector of atom names objects that define the mapping. The same name will be used both
   * as internal name and for PDB output. The given names therefore must comply to the PDB format
   * i.e. must always comprise four characters
   */
  AtomTypingBase(const std::vector<std::string> & atom_names);

  /// Virtual destructor to satisfy a compiler
  virtual ~AtomTypingBase() {}

  virtual core::index2 atom_type(const core::data::structural::PdbAtom &a) const {
    return atom_type(a.atom_name(),a.owner()->residue_type().code3);
  }

  /** @brief Returns an internal index for an atom identified by its PDB name, its residue name and a requested variant.
   * This default implementation ignores the residue type variant.
   * @param atom_name - name of the atom as given in the PDB file
   * @param residue_name - three letter code of the residue
   * @param variant - defines residue variant, e.g. C-terminal. This class doesn't recognise variants -
   * all atom types are treated as AtomTypingVariants::STANDARD. This however may be changed in a derived class
   */
  virtual core::index2 atom_type(const std::string & atom_name, const std::string & residue_name,
      const AtomTypingVariants variant = AtomTypingVariants::STANDARD) const {

#ifdef DEBUG
    if(internal_name_to_index.find(atom_name)==internal_name_to_index.end()) {
      l << utils::LogLevel::SEVERE << "Can't find atom type for reside::atom "<<residue_name<<":"<<atom_name
        <<"\nKnown atom types are:";
      int i=0;
      for (const auto &n : index_to_internal_name) l << (((++i) % 5 == 1) ? "\n" : "") << n;
    }
#endif
    return internal_name_to_index.at(atom_name);
  }

  /** @brief Returns integer atom type index for a given force field internal atom name
   *
   * @param atom_internal_name - internal name of an atom type, e.g. "CT" in AMBER force field
   * @return - integer index of the requested atom type
   */
  virtual core::index2 atom_type(const std::string & atom_internal_name) const { return internal_name_to_index.at(atom_internal_name); }

  /** @brief Returns a string that identifies an atom type in a given force field
   *
   * @param atom_id - integer index of the requested atom type
   * @return - string assigned to an atom type according to a force field
   */
  virtual const std::string & atom_internal_name(const core::index2 atom_id) const { return index_to_internal_name[atom_id]; }

  /** @brief Returns a string that identifies an atom type in a given force field
   *
   * @param atom - atom object
   * @return - string assigned to an atom type according to a force field
   */
  virtual const std::string & atom_internal_name(const core::data::basic::Vec3 & atom) const { return index_to_internal_name[atom.atom_type]; }

  /// Returns the number of atom types known to this typing
  virtual core::index2 count_atom_types() const { return index_to_internal_name.size(); }

  void show_known_atom_types(std::ostream & output) const;

protected:
  /// Mapping from  atom type name to atom type index
  std::unordered_map<std::string, core::index2> internal_name_to_index;
  /// Mapping from  atom type index to atom type name
  std::vector<std::string> index_to_internal_name;

  /// Default constructor is protected so only derived classes may initialize an empty typing
  AtomTypingBase()  {}

private:
  static utils::Logger l;
  void compile_types(const core::data::structural::Structure & structure);
};

}
}
#endif

/**
 * \example ex_UniquePdbTyping.cc
 */
