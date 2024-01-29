#ifndef CORE_DATA_STRUCTURAL_PdbAtom_H
#define CORE_DATA_STRUCTURAL_PdbAtom_H

#include <string>
#include <memory>
#include <functional>

#include <core/real.hh>
#include <core/index.hh>

#include <core/data/structural/PdbAtom.fwd.hh>
#include <core/data/structural/Residue.fwd.hh>
#include <core/chemical/AtomicElement.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/io/Pdb.hh>

#include <utils/string_utils.hh>

namespace core {
namespace data {
namespace structural {

using core::data::basic::Vec3;

/** @brief Method that attempts to convert an ungapped atom name to a PDB-style atom name.
 * @param atom_name - ungapped atom name, e.g. <code>"N"</code> or <code>"CB"</code>
 * @return correctly padded atom name, e.g. <code>" N  "</code> or <code>" CB "</code>
 */
std::string pdb_atom_name(const std::string &atom_name);

/** \brief Represents an atom coming from a PDB data.
 *
 * PdbAtom has its coordinates, chemical atom type (i.e. the element) as well as PDB-specific information: occupancy and temperature factor.
 *
 * To write a PdbAtom object in the PDB format, use <code>to_pdb_line()</code> method.
 */
class PdbAtom: public Vec3, public std::enable_shared_from_this<PdbAtom> {
public:

  /** \brief Creates a new atom.
   *
   * This constructor assigns (0,0,0) to coordinates.
   *
   * Example creates a GLY amino acid (non-hydrogen atoms only).
   * @code
   * PdbAtom n(1," N  ",7);
   * PdbAtom ca(2," CA "); // Carbon is the default element, we there is no need to provide atomic number
   * PdbAtom c(3," C  ");
   * PdbAtom o(4," O  ",8);
   * PdbAtom oxt(5," OXT", 8);
   * @endcode
   * @param id - integer index of an atom; according to the specification of the PDB format counting starts from 1
   * @param atom_name - four-character string (PDB-style) encoding the name of the atom, e.g. " CA ", " N  ", "1HB1"
   * @param element_index - atom number of the element, e.g. 6 for carbon (Carbon is the default element)
   */
  PdbAtom(const core::index4 id = 0, const std::string &atom_name = " CA ",
          const core::index2 element_index = core::chemical::AtomicElement::CARBON.z);

  /** \brief Creates a new atom from the given data.
   *
   * @code
   * PdbAtom n(2," N  ",0.0,0.0,0.0,1.0,21.78,7);
   * @endcode
   * @param id - integer index of an atom; according to the specification of the PDB format counting starts from 1
   * @param atom_name - four-character string (PDB-style) encoding the name of the atom, e.g. " CA ", " N  ", "1HB1"
   * @param cx - X coordinate of the atom
   * @param cy - Y coordinate of the atom
   * @param cz - Z coordinate of the atom
   * @param occupancy - occupancy value
   * @param b_factor - temperature factor
   * @param element_index - atom number of the element, e.g. 6 for carbon (Carbon is the default element)
   */
  PdbAtom(const core::index4 id, const std::string &atom_name,
		  const core::real cx, const core::real cy, const core::real cz, const core::real occupancy = 1.0,
		  const core::real b_factor = 99.99, const core::index2 element_index = core::chemical::AtomicElement::CARBON.z);

  /** \brief Makes a deep copy of this atom.
   *
   * The returned copy however does not belong to any residue!
   * @return a shared pointer to a newly made deep copy of this object
   */
  PdbAtom_SP clone() const;

  // ---------- Getters ----------

  /// Returns integer index of this atom
  core::index4 id() const { return id_; }

  /// Returns alternate locator for this atom
  char alt_locator() const { return alt_locator_; }

  /// Returns occupancy of this atom
  core::real occupancy() const { return occupancy_; }

  /// Returns temperature factor of this atom
  core::real b_factor() const { return b_factor_; }

  /// Returns the name of this atom
  const std::string &atom_name() const { return atom_name_; }

  /// Returns the index of the element this atom represents
  unsigned char element_index() const { return element_index_; }

  /// Returns true if this is a hetero-atom (<code>HETATM</code> line rather than <code>ATOM</code>)
  bool is_heteroatom() const { return is_heteroatom_; }

  // ---------- Setters ----------
  /// Sets integer index of this atom
  void id(const core::index4 id) { id_ = id; }

  /// Sets new atom name for this atom
  void atom_name(const std::string atom_name) { atom_name_ = atom_name; }

  /// Sets occupancy for this atom
  void occupancy(const core::real occ) { occupancy_ = occ; }

  /// Sets alternate locator for this atom
  void alt_locator(const char alt_loc) { alt_locator_ = alt_loc; }

  /// Sets temperature factor for this atom
  void b_factor(const core::real bf) { b_factor_ = bf; }

  /// Sets the new content for this atom instance
  inline void set(const core::index4 new_id, std::string &new_atom_name, const core::real new_x, const core::real new_y, const core::real new_z) {
    id_ = new_id;
    x = new_x;
    y = new_y;
    z = new_z;
    atom_name_ = new_atom_name;
  }

  /// Sets the new coordinates for this atom instance
  inline void set(const Vec3 & new_v) {
    x = new_v.x;
    y = new_v.y;
    z = new_v.z;
  }

  /// Sets a hetero-atom flag: ifset to true, this atom will be printed as <code>HETATM</code> line,  <code>ATOM</code> otherwise
  void is_heteroatom(const bool flag) { is_heteroatom_ = flag; }

	// ---------- Atom tree operations ----------
	/// Returns a const-pointer to the residue owning this atom
  const std::shared_ptr<Residue> owner() const {
    std::shared_ptr<Residue> r = owner_.lock();
    if (r) return r;
    else return nullptr;
  }

  /// Returns a pointer to the residue owning this atom
  std::shared_ptr<Residue> owner() {
    std::shared_ptr<Residue> r = owner_.lock();
    if (r) return r;
    else return nullptr;
  }

  /// Sets the new owner (i.e. a residue) that owns this atom
  void owner(std::shared_ptr<Residue> new_owner) { owner_ = new_owner; }

	// ---------- Other stuff ----------
  /** \brief Creates a PDB-formatted string from this atom.
   *
   * @returns a PDB string that encodes this atom, e.g:
   * <pre>ATOM    161  CB  VAL A  21     -14.040  14.093  -3.468  1.00 25.78           C
</pre>
   */
  std::string to_pdb_line() const;

  friend Residue;

  /// Less-than operator returns <code>true </code> when ID of this atom is lower than a.id()
  bool operator<(const PdbAtom & a) const { return id_ < a.id_; }

  /// Two atoms are equal if their IDs are equal
  bool operator==(const PdbAtom & a) const { return id_ == a.id_; }

private:
  core::index4 id_;
  std::string atom_name_;
  unsigned char element_index_;
  char alt_locator_;
  core::real occupancy_;
  core::real b_factor_;
  bool is_heteroatom_;
  std::weak_ptr<Residue> owner_;
};

/** @brief Two atoms are equal if their IDs are equal
 *  @param ai - a pointer to the first of the two atoms being compared
 *  @param aj - a pointer to the second of the two atoms being compared
 */
static inline bool operator==(const PdbAtom_SP & ai,const PdbAtom_SP & aj)  { return ai->id() == aj->id(); }

/** @brief Atom <code>ai</code> is lower than <code>aj</code> if and only if its ID number is lower
 *  @param ai - a pointer to the first of the two atoms being compared
 *  @param aj - a pointer to the second of the two atoms being compared
 */
//bool operator<(const PdbAtom_SP & ai,const PdbAtom_SP & aj);
static inline bool operator<(const PdbAtom_SP & ai,const PdbAtom_SP & aj)  { return ai->id() < aj->id(); }

}
}
}

namespace std {

/** @brief Hash of an atom is simply equal to this atom ID number
 */
template<>
struct hash<core::data::structural::PdbAtom> {

  /// Returns ID of the given atom
  std::size_t operator()(const core::data::structural::PdbAtom &a) const { return (a.id()); }
};

/** @brief Hash of an atom is simply equal to this atom ID number
 */
template<>
struct hash<core::data::structural::PdbAtom_SP> {

  /// Returns ID of the given atom
  std::size_t operator()(const core::data::structural::PdbAtom_SP &a) const { return (a->id()); }
};

}

#endif
