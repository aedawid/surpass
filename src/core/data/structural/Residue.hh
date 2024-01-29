/** @file Residue.hh
 *  @brief Defines Residue data structure and less-than operators for it.
 */
#ifndef CORE_DATA_STRUCTURAL_Residue_H
#define CORE_DATA_STRUCTURAL_Residue_H

#include <string>
#include <memory>
#include <algorithm>

#include <core/real.hh>
#include <core/index.hh>

#include <utils/exceptions/AtomNotFound.hh>
#include <utils/string_utils.hh> // for string_format()

#include <core/data/basic/Vec3.hh>

#include <core/data/structural/Residue.fwd.hh>
#include <core/data/structural/Chain.fwd.hh>
#include <core/data/structural/PdbAtom.fwd.hh>
#include <core/data/structural/structure_selectors.fwd.hh>

#include <core/chemical/Monomer.hh>

namespace core {
namespace data {
namespace structural {

/** \brief Single residue : an amino acid, nucleic base, ligand, ion, etc.
 *
 * Residue is a container for atoms.
 *
 * The following example reads a PDB file and tests if each of the amino acid residues has all the five backbone heavy atoms:
 * \include ex_Residue.cc
 */
class Residue : public std::enable_shared_from_this<Residue>, public std::vector<std::shared_ptr<PdbAtom>> {

friend std::ostream& operator<<(std::ostream &out, const Residue & r);
friend utils::Logger &operator <<(utils::Logger &logger, const Residue & r);

public:
  // ---------- C-tors ----------
  /// Constructor creates  a new residue of a given type
  Residue(const int id, const core::chemical::Monomer &residue_type) : id_(id), residue_type_(residue_type) {
	  reserve(residue_type.n_atoms);
	  insertion_code_ = ' ';
	  ss_type_ = 'C';
  }

  /// Constructor creates  a new residue of a given type defined by its 3-letter code
  Residue(const int id, const std::string &residue_type) : id_(id), residue_type_(core::chemical::Monomer::get(residue_type)) {
	  reserve(residue_type_.n_atoms);
	  insertion_code_ = ' ';
	  ss_type_ = 'C';
  }

  /// Constructor creates  a new residue of a given type defined by its one-letter code
  Residue(const int id, const char &residue_type) : id_(id), residue_type_(core::chemical::Monomer::get(residue_type)) {
	  reserve(residue_type_.n_atoms);
	  insertion_code_ = ' ';
	  ss_type_ = 'C';
  }

  /** @brief Creates a new Residue as a deep copy of this object.
   *
   * This method makes also a deep copy of atoms that belong to this residue providing that
   * they satisfy a given selector.
   *
   * @param which_atoms - clone also atoms which satisfy the given selector
   */
  Residue_SP clone(const AtomSelector &which_atoms) const;

  // ---------- Getters ----------
  /// Returns and integer identifying this residue
  inline int id() const { return id_; }

  /// Returns an insertion code character
  inline char icode() const { return insertion_code_; }

  /// Returns a character denoting secondary structure assigned to this residue.
  inline char ss() const { return ss_type_; }

  /// Returns a PDB-convention string identifying this residue (id + icode)
  inline std::string &residue_id() { return residue_id_; }

  /// Returns the chemical residue type, i.e. the Monomer object assigned to this residue
  inline const core::chemical::Monomer &residue_type() const { return residue_type_; }

  // ---------- Misc operations ----------
  /// Returns the number of atoms that belong to this residue
  inline core::index2 count_atoms() const { return size(); }

  /// Returns the number of heavy (i.e. non-hydrogen) atoms that belong to this residue
  core::index2 count_heavy_atoms() const;

  // ---------- Setters ----------
  /// Sets the new integer index identifying this residue
  inline void id(const int new_id) {
    id_ = new_id;
    residue_id_ = utils::string_format("%d%c", id_, insertion_code_);
  }

  /// Sets the new insertion code character
  inline void icode(const char new_icode) {
    insertion_code_ = new_icode;
    residue_id_ = utils::string_format("%d%c", id_, insertion_code_);
  }

	/// Sets secondary structure type for this residue
  inline void ss(const char new_ss) { ss_type_ = new_ss; }

  // ---------- Atom tree operations ----------
  /// Returns a const-pointer to the chain owning this residue
  const std::shared_ptr<Chain> owner() const {
    std::shared_ptr<Chain> c = owner_.lock();
    if (c) return c;
    else return nullptr;
  }

  /// Returns a pointer to the chain owning this residue
  std::shared_ptr<Chain> owner() {
    std::shared_ptr<Chain> c = owner_.lock();
    if (c) return c;
    else return nullptr;
  }

  /** @brief Returns a pointer to the residue that directly precedes this residue in the chain.
   *
   * If this residue is N-terminal, <code>nullptr</code> is returned
   */
  Residue_SP previous() {
    std::shared_ptr<Residue> c = previous_.lock();
    if (c) return c;
    else return nullptr;
  }

  /** @brief Returns a const-pointer to the residue that directly precedes this residue in the chain.
   *
   * If this residue is N-terminal, <code>nullptr</code> is returned
   */
  const Residue_SP previous() const {
    std::shared_ptr<Residue> c = previous_.lock();
    if (c) return c;
    else return nullptr;
  }

  /** @brief Defines a residue that precedes this residue in a chain.
   * @param next_residue - points to a residue that precedes this residue in a biomolecular chain
   */
  void previous(Residue_SP previous_residue) { previous_ = previous_residue; }

  /** @brief Returns a pointer to the residue that directly follows this residue in the chain.
   *
   * If this residue is C-terminal, <code>nullptr</code> is returned
   */
  Residue_SP next() {
    std::shared_ptr<Residue> c = next_.lock();
    if (c) return c;
    else return nullptr;
  }

  /** @brief Defines a residue that follows this residue in a chain.
   * @param next_residue - points to a residue that follows this residue in a biomolecular chain
   */
  void next(Residue_SP next_residue) { next_ = next_residue; }

  /** @brief Returns a const-pointer to the residue that directly follows this residue in the chain.
   *
   * If this residue is C-terminal, <code>nullptr</code> is returned
   */
  const Residue_SP next() const {
    Residue_SP c = next_.lock();
    if (c) return c;
    else return nullptr;
  }

  /// Sets the new owner (i.e. a chain) that owns this residue
  void owner(std::shared_ptr<Chain> new_owner) { owner_ = new_owner; }

  /// Adds an atom to this residue
  void push_back(PdbAtom_SP a);

  friend Chain;

  /** @brief Computes the closest distance between any atom of this residue and any atom from another residue
   *
   * @param another_residue - the second residue for the distance calculation
   * @return a distance value
   *
   * The example program reads a PDB file and looks for a specific ligand. Then it finds
   * all residues from a structure whose <code>min_distance()</code> is shorter than a given cutoff.
   *
   * \include ./ex_LigandNeighbors.cc
   */
  const core::real min_distance(const Residue & another_residue) const;

  /** @brief Computes the closest distance between any atom of this residue and any atom from another residue
   *
   * @param another_residue - the second residue for the distance calculation
   * @return a distance value
   */
  const core::real min_distance(const Residue_SP another_residue) const { return min_distance(*another_residue); }

  /** @brief Finds atom by name; throws an exception if not found
   *
   * @param atom_name - the name of an atom to look for
   * @return a pointer to an atom
   */
  const PdbAtom_SP find_atom_safe(const std::string & atom_name) const {
    if (find_atom(atom_name) != nullptr) return find_atom(atom_name);
    throw utils::exceptions::AtomNotFound(atom_name, residue_id_);
  }

  /** @brief Finds atom by name.
   *
   * @param atom_name - the name of an atom to look for
   * @return a pointer to an atom; <code>nullptr</code> when not found
   */
  PdbAtom_SP find_atom(const std::string & atom_name);

  /** @brief Finds atom by name.
   *
   * @param atom_name - the name of an atom to look for
   * @return a pointer to an atom; <code>nullptr</code> when not found
   */
  PdbAtom_SP find_atom(const char* atom_name)  { return find_atom(std::string(atom_name)); }

  /** @brief Finds atom by name.
   *
   * @param atom_name - the name of an atom to look for
   * @return a pointer to an atom; <code>nullptr</code> when not found
   */
  const PdbAtom_SP find_atom(const std::string & atom_name) const;

  /** @brief Finds atom by name.
   *
   * @param atom_name - the name of an atom to look for
   * @return a pointer to an atom; <code>nullptr</code> when not found
   */
  const PdbAtom_SP find_atom(const char* atom_name) const  { return find_atom(std::string(atom_name)); }

  /// Sort atoms of this residue by their ID
  void sort();

private:
  int id_;    ///< Integer identifying this residue
  const core::chemical::Monomer &residue_type_; ///< The chemical type of this residue
  char insertion_code_;    ///< PDB-style insertion code
  char ss_type_; ///< secondary structure type
  std::string residue_id_;
  std::weak_ptr<Chain> owner_;
  std::weak_ptr<Residue> previous_;
  std::weak_ptr<Residue> next_;
};

/** @brief less-then operator provides PDB-like ordering.
 *
 * A given residue \f$r_i\f$ should go before residue \f$r_j\f$ if:
 *     - the ID of \f$r_i\f$ is less than the ID of \f$r_j\f$
 *     - IDs of the two residues are equal and the insertion code of \f$r_i\f$ is lexically before the insertion code of \f$r_j\f$
 */
bool operator<(const Residue & ri, const Residue & rj);

/** @brief less-then operator provides PDB-like ordering.
 *
 * A given residue \f$r_i\f$ should go before residue \f$r_j\f$ if:
 *     - the ID of \f$r_i\f$ is less than the ID of \f$r_j\f$
 *     - IDs of the two residues are equal and the insertion code of \f$r_i\f$  lexicographically precedes the insertion code of \f$r_j\f$
 */
bool operator<(const Residue_SP ri, const Residue_SP rj);

}
}
}

#endif
/**
 * \example ex_Residue.cc
 * \example ./ex_LigandNeighbors.cc
 */
