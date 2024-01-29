/** \file structure_selectors.hh
 * @brief Provides functors that selects a fragment of a biomolecular structure.
 *
 * Selector object tests whether a given element of a molecular structure (an atom, a residue or a chain)
 * satisfies particular condition. The following example shows how to use atom selectors:
 *
 * \include ex_AtomSelector.cc
 *
 * and here is another example that uses SelectResidueRange selector:
 *
 * \include ex_SelectResidueRange.cc
 *
 * The functors may be used in combination with STL algorithms and containers to
 * remove a part of a structure or to iterate over a selection. A few examples are given below:
 * @code
 * // Calculate the number of atoms in the first ten residues from a Structure_SP object
 * Structure_SP str;  // init the structure somehow
 * core::index4 n = std::count_if(str->first_atom(),str->last_atom(),SelectResidueRange("1-10"));
 *
 * // Erase residues from -1 to 10A (inclusive) from chain A
 * Chain_SP chain_a = str->get_chain('A');
 * chain_a->erase(std::remove_if(chain_a->begin(), chain_a->end(), SelectResidueRange("-1-10B")), chain_a->end());
 *
 * // Erase amino acid residues that lack alpha-carbon
 * chain_a->erase(std::remove_if(chain_a->begin(), chain_a->end(), ResidueHasCA()), chain_a->end());
 * @endcode
 *
 * The following example shows how to use selectors in cloning a part of a structure:
 * \include ex_StructureSelector.cc
 */

#ifndef CORE_DATA_STRUCTURAL_structure_selectors_H
#define CORE_DATA_STRUCTURAL_structure_selectors_H

#include <string>
#include <limits>
#include <iostream>
#include <algorithm> // for std::replace
#include <cctype>

#include <core/index.hh>

#include <core/data/structural/structure_selectors.fwd.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/PdbAtom.hh>

#include <utils/string_utils.hh>

namespace core {
namespace data {
namespace structural {

/** @brief A base class for any selector working on atoms.
 *
 * AtomSelector is a virtual class, do not instantiate. Use other selector instead. The example below demonstrates
 * a few of them:
 *
 * \include ex_AtomSelector.cc
 */
class AtomSelector {
public:
  AtomSelector() : selection_string_("") {}
  virtual bool operator()(const PdbAtom & c) const { return true; }
  bool operator()(const PdbAtom_SP c) const { return operator()(*c); }
  virtual ~AtomSelector() { }
  virtual const std::string & selector_string() const { return selection_string_; }
  virtual void set(const std::string & new_selection) { }

private:
  const std::string selection_string_;
};

/** @brief Always returns true.
 */
class SelectAllAtoms : public AtomSelector {
public:
  /** @brief Always returns true.
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const { return true; }

  /** @brief Returns the selection string which is always <code>"*"</code>
   */
  virtual const std::string & selector_string() const { return selection_string_; }

private:
  static const std::string selection_string_;
};

/** @brief Returns true if a given atom is an \f$\alpha\f$-carbon.
 *
 * Example that shows how to use AtomSelector objects, in particular an IsCA instance:
 * \include ex_AtomSelector.cc
 */
class IsCA : public AtomSelector {
public:
  /** @brief Returns true if a given atom is an \f$\alpha\f$-carbon.
   *
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const { return a.atom_name() == selection_string_; }

  /** @brief Returns the selection string which is always <code>" CA "</code>
   */
  virtual const std::string & selector_string() const { return selection_string_; }

private:
  static const std::string selection_string_;
};

/** @brief Returns true if a given atom is an \f$\beta\f$-carbon.
 *
 * Example that shows how to use AtomSelector objects, in particular an IsCB instance:
 * \include ex_AtomSelector.cc
 */
class IsCB : public AtomSelector {
public:
  /** @brief Returns true if a given atom is a \f$\beta\f$-carbon.
   *
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const { return a.atom_name() == selection_string_; }

  /** @brief Returns the selection string which is always <code>" CB "</code>
   */
  virtual const std::string & selector_string() const { return selection_string_; }

private:
  static const std::string selection_string_;
};

/** @brief Returns true for an atom which appears to be hydrogen.
 * This selector tests whether the element type or atom name suggests a given atom is hydrogen
 */
class IsHydrogen : public AtomSelector {
public:

  IsHydrogen()  {}

  /** @brief Returns true if a given atom is a hydrogen
   *
   * @param a - reference to an atom
   * @return true is any of the following is true:
   *   - the given atom's element index equals to 1 (i.e. hydrogen)
   *   - the given atom's name contains 'H' at the second position (according PDB convention this should be hydrogen)
   */
  virtual bool operator()(const PdbAtom & a) const {
    return (a.atom_name()[1] == 'H') || (a.element_index() == 1);
  }
};

/** @brief Returns true always when <code>IsHydrogen</code> would return false
 */
class NotHydrogen : public AtomSelector {
public:

  NotHydrogen() : is_h() {}

  /** @brief Returns true if a given atom is not a hydrogen
   *
   * @param a - reference to an atom
   * @return true is any of the following is true:
   *   - the given atom's element index is not 1 (i.e. not a hydrogen)
   *   - the given atom's name does not contain 'H' at the second position (according PDB convention this should be hydrogen)
   */
  virtual bool operator()(const PdbAtom & a) const { return ! is_h(a); }
private :
  IsHydrogen is_h;
};

/** @brief Returns true if a given atom is an alternate location
 *
 */
class IsAlternateLocation : public AtomSelector {
public:
  /** @brief Returns true if a given atom is an alternate location
   *
   * @param a - reference to an atom
   * @return true if the locator char for this atom is neither ' ' nor 'A' character
   */
  virtual bool operator()(const PdbAtom & a) const {
    return (a.alt_locator() !=  ' ') && (a.alt_locator() !=  'A');
  }
};

/** @brief Returns true if a given atom belongs to protein backbone.
 *
 * Example that shows how to use AtomSelector objects, in particular an IsBB instance:
 * \include ex_AtomSelector.cc
 *
 */
class IsBB : public AtomSelector {
public:

  IsBB() : selection_string_("_N__+_CA_+_O__+_C__+_H__+_HA_+_HA1+_HA2+_HA3+_OXT") {}

  /** @brief Returns true if a given atom belongs to a protein backbone
   *
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const {
    return (a.atom_name() == " CA ") || (a.atom_name() == " N  ") || (a.atom_name() == " C  ")
        || (a.atom_name() == " O  ") || (a.atom_name() == " OXT") || (a.atom_name() == " H  ")
        || (a.atom_name() == " HA ") || (a.atom_name() == " HA3") || (a.atom_name() == " HA1")
        || (a.atom_name() == " HA2");
  }

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selection_string_; }

private:
  const std::string selection_string_;
};

/** @brief Returns true if a given atom belongs to protein backbone or is a \f$\beta\f$-carbon.
 *
 * Example that shows how to use AtomSelector objects, in particular an IsBBCB instance:
 * \include ex_AtomSelector.cc
 */
class IsBBCB : public AtomSelector {
public:

  IsBBCB() : selection_string_("_N__+_CA_+_O__+_C__+_H__+_HA_+_HA1+_HA2+_HA3+_OXT+_CB_") {}

  /** @brief Returns true if a given atom belongs to a protein backbone or a \f$\beta\f$-carbon.
   *
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const {
    return (a.atom_name() == " CA ") || (a.atom_name() == " N  ") || (a.atom_name() == " C  ") ||
           (a.atom_name() == " O  ") || (a.atom_name() == " CB ") || (a.atom_name() == " OXT") ||
           (a.atom_name() == " H  ") || (a.atom_name() == " HA ") || (a.atom_name() == " HA3") ||
           (a.atom_name() == " HA1") || (a.atom_name() == " HA2");
  }

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selection_string_; }

private:
  const std::string selection_string_;
};

/** @brief Matches an atom by its name, e.g. <code>"NZ1 "</code>.
 *
 * Example that shows how to use AtomSelector objects, in particular an IsNamedAtom instance:
 * \include ex_AtomSelector.cc
 */
class IsNamedAtom : public AtomSelector {
public:

  IsNamedAtom()  {}

  IsNamedAtom(const std::string & atom_name) :
      atom_name_(atom_name), padded_atom_name_(atom_name) {
    std::replace(atom_name_.begin(), atom_name_.end(), '_', ' ');
  }

  /** @brief Returns true if a given atom is named as declared
   *
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const;

  /** @brief Returns the selection string which is equal to the name of selected atom
   */
  virtual const std::string & selector_string() const { return padded_atom_name_; }

  /** @brief Sets the new atom name to be selected
   *
   * @param atom_name - atom name
   */
  void set(const std::string & atom_name);

private:
  std::string atom_name_;
  std::string padded_atom_name_;
};


/** @brief Matches an atom by its name, e.g. <code>"NZ1 "</code>.
 *
 * Example that shows how to use AtomSelector objects, in particular an IsNamedAtom instance:
 * \include ex_AtomSelector.cc
 */
class IsElement : public AtomSelector {
public:

  IsElement(const std::string & element_name) :
    element_index_(core::chemical::AtomicElement::by_symbol(element_name).z) {}

  /** @brief Returns true if a given atom is named as declared
   *
   * @param a - shared pointer to an atom
   */
  virtual bool operator()(const PdbAtom & a) const;

  /** @brief Returns the selection string which is equal to the name of selected atom
   */
  virtual const std::string & selector_string() const {
    return core::chemical::AtomicElement::periodic_table[element_index_].name;
  }

  /** @brief Sets the new atom name to be selected
   *
   * @param atom_name - atom name
   */
  void set(const std::string & atom_name);

private:
  core::index2 element_index_;
};


/** @name Selectors for chains and residues
 */
///@{


/** @brief A base class for all residue selectors.
 *
 * By default a ResidueSelector instance selects every residue and atom
 */
class ResidueSelector : public AtomSelector {
public:

  /** @brief Tests whether a given residue matches the selection.
   *
   * @param r - points to the tested residue
   * @return always true
   */
  virtual bool operator()(const Residue & r) const { return true; }

  /// Calls bool operator()(const Residue & r)
  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /** @brief Tests whether a given atom belongs to a selected residue
   *
   * @param a - points to a tested atom
   * @return true if and only if this residue selector returns true for the residue that owns the given atom
   */
  virtual bool operator()(const PdbAtom & a) const;

  /** @brief Redefines selection of this selector
   *
   * @param new_selection - string which defines the new selection (details depend on derived class)
   */
  virtual void set(const std::string & new_selection) { }

  virtual ~ResidueSelector() {}
};


/** @brief Returns true when a given residue has an alpha carbon atom.
 */
class ResidueHasCA : public ResidueSelector {
public:

  /** @brief Returns true when this residue has an alpha carbon atom.
   *
   * @param r - points to the tested residue
   * @return true if the given residue has C\f$\alpha\f$ atom
   */
  virtual bool operator()(const Residue & r) const;

  /** @brief Tests whether a given atom belongs to a residue with C\f$\alpha\f$ atom.
   *
   * @param a - points to a tested atom (does not have to be C\f$\alpha\f$ for a successful selection)
   * @return true if the given residue has C\f$\alpha\f$ atom
   */
  virtual bool operator()(const PdbAtom & a) const { return operator ()(*(a.owner())); }

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~ResidueHasCA() {}
};

/** @brief Creates a respective residue selector based on a string.
 *
 * @param selection - string definition of a selector. Valid options are:
 *      - <strong>*</strong> - selects everything (an instance of <code>ResidueSelector</code> is returned)
 *      - <strong>aa</strong> or <strong>AA</strong> - select amino acid residues  (an instance of <code>IsAA</code> is returned)
 *      - <strong>nt</strong> or <strong>NT</strong> - select nucleotide residues  (an instance of <code>IsNT</code> is returned)
 *      - <strong>1-120</strong>,  any range of integer residue IDs  (an instance of <code>SelectResidueRange</code> is returned)
 * @return a newly created residue selector
 */
ResidueSelector_SP residue_selector_from_string(const std::string & selection);

/** @brief Returns true when the given residue's CA atom has neighbors in the proper distance.
 *
 * This selector tests the distance between the CA atom of this residue and the CA of the previous residue. Similarly,
 * distance between the next CA is tested. Both distances must be shorter than 4.2 Angstroms (this is the default cutoff,
 * which may be changed by the constructor. For terminal residues (e.g. the N terminal amino acid) only one distance is tested
 */
class ProperlyConnectedCA : public ResidueSelector {
public:

  ProperlyConnectedCA(core::real cutoff = 4.2) : cutoff_(cutoff) {}

  /** @brief Returns true when this residue has an alpha carbon atom.
   *
   * @param r - points to the tested residue
   * @return true if the given residue has C\f$\alpha\f$ atom
   */
  virtual bool operator()(const Residue & r) const;

  /** @brief Tests whether a given atom belongs to a residue with C\f$\alpha\f$ atom.
   *
   * @param a - points to a tested atom (does not have to be C\f$\alpha\f$ for a successful selection)
   * @return true if the given residue has C\f$\alpha\f$ atom
   */
  virtual bool operator()(const PdbAtom & a) const { return operator ()(*(a.owner())); }

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~ProperlyConnectedCA() {}

private:
  core::real cutoff_;
  PdbAtom_SP get_ca(const Residue & r) const;
};

/** @brief Returns true when the given residue is a terminal one.
 *
 * This selector tests whether <code>next()</code> or <code>previous()</code> methods return a <code>nullptr</code>
 * Returns true also when the following (or preceding) residue is of a different type than this one.
 */
class ResidueIsTerminal : public ResidueSelector {
public:

  /** @brief Returns true when both <code>IsFirstResidue</code> and <code>IsLastResidue</code> returned true
   *
   * @param r - points to the tested residue
   * @return true if the given residue is a terminal one.
   */
  virtual bool operator()(const Residue &r) const;

  virtual bool operator()(const PdbAtom & a) const { return operator ()(*(a.owner())); }

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~ResidueIsTerminal() {}
};

/** @brief Returns true when the given residue is the very first one in its chain.
 *
 * This selector tests whether code>previous()</code> methods return a <code>nullptr</code>
 * Returns true also when the following residue is of a different type than this one.
 */
class IsFirstResidue : public ResidueSelector {
public:

  /** @brief Returns true when this residue is the very first residue in its chain.
   *
   * This method returs also true if the preceding residue is of a different type than this one.
   * The latter condition eliminates e.g ligands stored in an amino acid chain
   *
   * @param r - points to the tested residue
   * @return true if the given residue is a terminal one.
   */
  virtual bool operator()(const Residue &r) const;

  virtual bool operator()(const PdbAtom & a) const { return operator ()(*(a.owner())); }

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~IsFirstResidue() {}
};

/** @brief Returns true when the given residue is the very last one in its chain.
 *
 * This selector tests whether code>previous()</code> methods return a <code>nullptr</code>
 * Returns true also when the following residue is of a different type than this one.
 */
class IsLastResidue : public ResidueSelector {
public:

  /** @brief Returns true when this residue is the very first residue in its chain.
   *
   * This method returs also true if the following residue is of a different type than this one.
   * The latter condition eliminates e.g ligands stored in an amino acid chain
   *
   * @param r - points to the tested residue
   * @return true if the given residue is a terminal one.
   */
  virtual bool operator()(const Residue &r) const;

  virtual bool operator()(const PdbAtom & a) const { return operator ()(*(a.owner())); }

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~IsLastResidue() {}
};

/** @brief Returns true when a given residue has all heavy backbone atoms.
 *
 * This selector selects residues if and only if a residue contains <code>N</code>, <code>CA</code>, <code>C</code> and <code>O</code> atoms.
 */
class ResidueHasBB : public ResidueSelector {
public:

  /** @brief Returns true when this residue has all backbone heavy atoms
   *
   * @param r - points to the tested residue
   * @return true if the given residue has <code>N</code>, <code>CA</code>, <code>C</code> and <code>O</code> atoms.
   */
  virtual bool operator()(const Residue & r) const;

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~ResidueHasBB() {}
};

/** @brief Returns true when a given residue has all heavy backbone atoms and beta carbon.
 *
 * This selector selects residues if and only if a residue contains
 * <code>N</code>, <code>CA</code>, <code>C</code>, <code>O</code>  and <code>CB</code>atoms.
 */
class ResidueHasBBCB : public ResidueSelector {
public:

  /** @brief Returns true when this residue has all backbone heavy atoms
   *
   * @param r - points to the tested residue
   * @return true if the given residue has <code>N</code>, <code>CA</code>, <code>C</code>, <code>O</code>  and <code>CB</code>atoms atoms.
   */
  virtual bool operator()(const Residue & r) const;

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~ResidueHasBBCB() {}
};

/** @brief Returns true when a given residue has all heavy atoms (both backbone and side chain).
 *
 * In order to satisfy this selector, a residue must have the correct number of heavy atoms,
 * as defined in <code>monomers.txt</code> file decreased by one (due to the leaving oxygen).
 * For N terminal amino acid residues as well as for 5' terminal bases the atom count must match exactly.
 */
class ResidueHasAllHeavyAtoms : public ResidueSelector {
public:

  /** @brief Returns true when this residue has all heavy atoms (both backbone and side chain).
   *
   * @param r - points to the tested residue
   * @return true if the given residue is complete
   */
  virtual bool operator()(const Residue & r) const;

  virtual bool operator()(const Residue_SP r) const { return operator()(*r); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~ResidueHasAllHeavyAtoms() {}

  core::index1 expected_heavy_atoms(const Residue & r) const;
};

/** @brief Returns true when this residue is an amino acid.
 */
class IsAA : public ResidueSelector {
public:

  /** @brief Returns true when this residue is an amino acid.
   *
   * @param r - points to the tested residue
   * @return true if the given residue  is an amino acid.
   */
  virtual bool operator()(const Residue & r) const;

  /** @brief Tests whether a given atom belongs to an amino acid residue
   *
   * @param a - points to a tested atom
   * @return true if the owning residue is an amino acid.
   */
  virtual bool operator()(const PdbAtom & a) const { return operator ()(*a.owner()); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  /// Does nothing here
  virtual ~IsAA() {}
};

/** @brief Returns true when this residue is an <strong>aromatic</strong> amino acid.
 */
class IsAromaticAA : public IsAA {
public:

  /** @brief Returns true when this residue is an amino acid.
   *
   * @param r - points to the tested residue
   * @return true if the given residue  is an amino acid.
   */
  virtual bool operator()(const Residue & r) const;

  virtual ~IsAromaticAA() {}
};

/** @brief Returns true when this residue is an <strong>aromatic</strong> amino acid.
 */
class SelectResidueByName : public ResidueSelector {
public:

  /** @brief Creates a selector that selects residues of <code>selected_code3</code>
   *
   * @param selected_code3 - three-letter string denoting a selected residue type, e.g. <code>ALA</code> or <code>HEM</code>
   */
  SelectResidueByName(const std::string & selected_code3);

  /** @brief Creates a selector that selects residues of given codes
   *
   * @param selected_code3 - three-letter strings denoting a selected residue types,
   *   e.g. <code>{"TRP", "TYR", "PHE", "HIS"}</code> selects all aromatic amino acids
   */
  SelectResidueByName(const std::vector<std::string> & selected_code3);

  /** @brief Returns true when this residue's name matches
   *
   * @param r - points to the tested residue
   * @return true if the given residue  is an amino acid.
   */
  virtual bool operator()(const Residue & r) const;

  virtual ~SelectResidueByName() {}
private:
  std::vector<std::string> matching_code3;
};


/** @brief Returns true when this residue is a nucleic acid.
 */
class IsNT : public ResidueSelector {
public:

  /** @brief Returns true when this residue is a nucleic acid.
   *
   * @param r - points to the tested residue
   * @return true if the given residue  is a nucleic acid.
   */
  virtual bool operator()(const Residue & r) const;

  /** @brief Tests whether a given atom belongs to a nucleic acid.
   *
   * @param a - points to a tested atom
   * @return true if the owning residue is a nucleic acid.
   */
  virtual bool operator()(const PdbAtom & a) const { return operator ()(*a.owner()); }

  /// Does nothing here
  virtual void set(const std::string & new_selection) { }

  virtual ~IsNT() {}
};

/** @brief Selects a range of residues within a chain.
 *
 * The functor returns true if an argument residue belongs to a given range, defined by
 * <code>residue_id</code> values.
 * The selection string provides two residue id values separated with a dash, e.g.
 * <code>1-20</code> or <code>-5-110</code> (note the first negative index value)
 * Alternatively, '*' selects all the residues.
 * Here is an example that uses SelectResidueRange selector:
 *
 * \include ex_SelectResidueRange.cc
 */
class SelectResidueRange : public ResidueSelector {
public:
  friend std::ostream & operator<<(std::ostream &out, const SelectResidueRange &selector);

  /** @brief Bare constructor creates a selector that selects everything
   */
  SelectResidueRange() {}

  /** @brief Constructor reads the selected range from a given string.
   *
   * The selection string provides two residue id values separated with a dash, e.g.
   * <code>1-20</code> or <code>-5-110</code> (note the first negative index value)
   * The '*' character selects all the residues.
   *
   * @param selector - the selection string
   */
  SelectResidueRange(const std::string & selector) { set(selector); }

  /** @brief Sets the new selection for this object
   *
   * The method changes the selection selected by this object.
   *
   * @param residue_from - the index of the first selected residue
   * @param residue_to - the index of the last selected residue (inclusive)
   * @param icode_from - icode  of the first selected residue
   * @param icode_to - icode of the last selected residue
   */
  void set(const int residue_from, const int residue_to, const char icode_from = ' ', const char icode_to = ' ') {

    first_residue = residue_from;
    last_residue = residue_to;
    first_residue_icode = icode_from;
    last_residue_icode = icode_to;
    update_selector_string();
  }

  /** @brief Sets the new selection range  definition
   *
   * @param selector - the selection string provides two residue id values separated with a dash, e.g.
   * <code>1-20</code> or <code>-5-110</code> (note the first negative index value)
   */
  void set(const std::string & selector);

  /** @brief Tests whether a given residue matches the selection.
   *
   * @param r - points to the tested residue
   * @return true when a given residue belongs to the selected range
   */
  virtual bool operator()(const Residue & r) const;

 /** @brief Tests whether a given residue matches the selection.
  *
  * @param r - points to the tested residue
  * @return true when a given residue belongs to the selected range
  */
  virtual bool operator()(const Residue_SP r) const { return operator()(*r); };

  /** @brief Tests whether a given atom matches the selection.
   *
   * @param r - points to the tested atom
   * @return true when a given atom belongs to the selected residue range
   */
  inline bool operator()(const PdbAtom & a) const { return operator()(*a.owner()); }

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selector; }

private:
  int first_residue = 0;
  int last_residue = -1;
  char first_residue_icode = ' ';
  char last_residue_icode = ' ';
  std::string selector = "*";

  void update_selector_string() {

    std::stringstream ss;
    ss << first_residue;
    if (first_residue_icode != ' ') ss << first_residue_icode;
    ss << '-';
    ss << last_residue;
    if (last_residue_icode != ' ') ss << last_residue_icode;
    selector = ss.str();
  }
};

/** @brief Selects a particular chain.
 *
 * The functor returns true if an argument chain has particular <code>chain_id</code>.
 * Alternatively, '*' selects all the chains in a structure object
 */
class ChainSelector : public ResidueSelector {
public:
  friend std::ostream & operator<<(std::ostream &out, const ChainSelector &selector);

  /** @brief Creates the selector for a particular <code>chain_id</code> value
   * @param chain_id - a chain identifier (single character). Character '*' selects all chains.
   */
  ChainSelector(const char chain_id = '*') : chain_id_(chain_id){ selector[0] = chain_id; }

  /** @brief Sets the new value for <code>chain_id</code>  field
   * @param c - a  <code>chain_id</code> value to be used by this functor
   */
  void set(const char chain_id) {
    chain_id_ = chain_id;
    selector[0] = chain_id;
  }

  /** @brief Sets the new value for <code>chain_id</code>  field
   * @param c - a  <code>chain_id</code> value to be used by this functor (only the first character from the given string is used)
   */
  virtual void set(const std::string & chain_id);

  /** @brief Returns true if the id of the given chain equals to the <code>chain_id</code>  value declared in this object.
   * @param c - a chain to be filtered
   */
  virtual bool operator()(const Chain & c) const;

  virtual bool operator()(const Chain_SP c) const { return operator()(*c); }

  virtual inline bool operator()(const Residue & r) const { return operator ()(*r.owner());}

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selector; }

private:
  char chain_id_;
  std::string selector = "*";
};

/** @brief Selector that selects every atom from every residue and chain.
 */
class SelectEverything : public ChainSelector {
public:

  /// Selects every atom - always returns true
  virtual inline bool operator()(const PdbAtom & a) const { return true; }

  /// Selects every chain - always returns true
  virtual inline bool operator()(const Chain & c) const { return true; }

  /// Selects every residue - always returns true
  virtual inline bool operator()(const Residue & r) const { return true;}

private:
  std::string selector = "*";
};

/** @brief Selector that combines one selector for chain(s) and one for residues.
 * Valid selection strings are for example:
 *      - <code>B:1-20</code> - selects residues from 1 to 20 (both inclusive) of chain B
 *      - <code>A</code> - selects the whole chain A
 */
class SelectChainResidues : public ChainSelector {
public:
  friend std::ostream & operator<<(std::ostream &out, const SelectChainResidues &selector);

  /** @brief Selector that combines one selector for chain(s) and one for residues
   *
   * @param chain_selector - object used to select chain(s)
   * @param residue_selector - object used to select residues
   */
  SelectChainResidues(ChainSelector_SP chain_selector, ResidueSelector_SP residue_selector) :
      chain_selector(chain_selector), residue_selector(residue_selector) {
    selector_ = chain_selector->selector_string()+":"+residue_selector->selector_string();
  }

  /** @brief Selector that combines one selector for chain(s) and one for residues
   *
   * @param selector - the selection string provides two selectors separated with a colon, e.g.
   * <code>B:1-20</code>. To select the whole chain, say simply <code>C:</code>
   */
  SelectChainResidues(const std::string & selector) { set(selector); }

  virtual void set(const std::string & selector);

  virtual inline bool operator()(const Chain & c) const { return (*chain_selector)(c); }

  virtual inline bool operator()(const Residue & r) const { return (*residue_selector)(r); }

  inline bool operator()(const Residue_SP r) const { return (*residue_selector)(*r); }

  /// Selects an atom if both the chain and the residue it belongs to is selected by this selector
  virtual inline bool operator()(const PdbAtom &a) const {
    return (*chain_selector)(a.owner()->owner()) && (*residue_selector)(a.owner());
  }

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selector_; }

private:
  ChainSelector_SP chain_selector;
  ResidueSelector_SP residue_selector;
  std::string selector_ = "*:*";
};

/** @brief Selector that combines one selector for chain(s) and one for residues
 */
class SelectChainResidueAtom : public ChainSelector {
public:

  friend std::ostream & operator<<(std::ostream &out, const SelectChainResidueAtom &selector);

  /** @brief Selector that combines one selector for chain(s) and one for residues
   *
   * @param chain_selector - object used to select chain(s)
   * @param residue_selector - object used to select residues
   * @param atom_selector - object used to select atoms
   */
  SelectChainResidueAtom(ChainSelector_SP chain_selector, ResidueSelector_SP residue_selector, AtomSelector_SP atom_selector) :
      chain_selector(chain_selector), residue_selector(residue_selector), atom_selector(atom_selector) {
    selector = chain_selector->selector_string() + ":" + residue_selector->selector_string() + ":"
        + atom_selector->selector_string();
  }

  /** @brief Selector that combines one selector for chain(s) and one for residues
   *
   * @param selector - the selection string provides two selectors separated with a colon, e.g.
   * <code>B:1-20</code>. To select the whole chain, say simply <code>C:</code>. To select all chains (or all residues), say '*'
   * For instance the selector <code>*:*</code> selects the whole structure
   */
  SelectChainResidueAtom(const std::string & selector) { set(selector); }

  virtual void set(const std::string & selector);

  virtual inline bool operator()(const Chain & c) const { return (*chain_selector)(c); }

  /** @brief Check if a given residue is selected by this object
   *
   * This call tests <strong>both</strong> the residue and the chain part of this selector
   * @tparam S - the type of selectors being aggregated, e.g. SelectChainResidues, SelectChainResidueAtom or SelectChain
   */
  virtual bool operator()(const Residue & r) const {
    return (*residue_selector)(r) && (*chain_selector)(*r.owner());
  }

  virtual bool operator()(const PdbAtom & a) const {

    return (*atom_selector)(a) && (*residue_selector)(*a.owner()) && (*chain_selector)(*(a.owner()->owner()));
  }

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selector; }

private:
  ChainSelector_SP chain_selector;
  ResidueSelector_SP residue_selector;
  AtomSelector_SP atom_selector;
  std::string selector;
};

/** @brief Selector that aggregates several selectors.
 *
 * The selectors being aggregated must be of the same type.
 * @tparam S - the type of selectors being aggregated, e.g. SelectChainResidues, SelectChainResidueAtom or SelectChain
 */
template<typename S>
class CompositeSelector : public AtomSelector {
public:
  template<typename T>
  friend std::ostream & operator<<(std::ostream &out, const CompositeSelector<T> & selector);

  /** @brief Bare constructor creates an empty selector that does not select anything.
   *
   * Selector components may be added by <code> add(S & selector)</code> method.
   */
  CompositeSelector(const char separator = ',') : separator(separator) { }

  /** @brief Creates a selector that may select several fragments
   *
   * @param selector - the selection string, e.g.
   * <code>A:20:40,B:1-20</code>.
   */
  CompositeSelector(std::string & multiselector,const char separator = ',') : separator(separator) {

    utils::trim(multiselector);
    std::replace(multiselector.begin(), multiselector.end(), '\n', ',');
    std::vector<std::string> tokens = utils::split(multiselector, ',');
    for (std::string & s : tokens)
      add(s);
  }

  size_t count_selectors() {
    return selectors.size();
  }

  void add(std::shared_ptr<S> & selector) {
    selectors.push_back(selector);
    if(selector_.size()>0) selector_ += separator;
    selector_ += selectors.back()->selector_string();
  }

  void add(const std::string & selector_string) {
    selectors.push_back(std::make_shared<S>(selector_string));
    if(selector_.size()>0) selector_ += separator;
    selector_ += selectors.back()->selector_string();
  }

  void add(const char* selector_string) {
    std::string s(selector_string);
    selectors.push_back(std::make_shared<S>(s));
    if(selector_.size()>0) selector_ += separator;
    selector_ += selectors.back()->selector_string();
  }

  inline bool operator()(const Chain & c) const {
    bool out = false;
    for (const  std::shared_ptr<S> s : selectors)
      out = out || (*s)(c);
    return out;
  }

  inline bool operator()(const Residue & r) const {
    bool out = false;
    for (const std::shared_ptr<S> s : selectors) out = out || (*s)(r);
    return out;
  }

  virtual bool operator()(const PdbAtom & a) const {
    bool out = false;
    for (const AtomSelector_SP s : selectors) out = out || s->operator()(a);
    return out;
  }

  /** @brief Returns the selection string.
   */
  virtual const std::string & selector_string() const { return selector_; }

private:
  std::vector<std::shared_ptr<S>> selectors;
  std::string selector_;
  const char separator;
};
/// @}

std::ostream & operator<<(std::ostream &out, const AtomSelector &selector);

std::ostream & operator<<(std::ostream &out, const ChainSelector &selector);

std::ostream & operator<<(std::ostream &out, const SelectResidueRange &selector);

std::ostream & operator<<(std::ostream &out, const SelectChainResidues &selector);

std::ostream & operator<<(std::ostream &out, const SelectChainResidueAtom &selector);

template<typename S>
std::ostream & operator<<(std::ostream &out, const CompositeSelector<S> &selector) {

  if (selector.selector_.size() == 0) {
    out << "<CompositeSelector is empty>";
    return out;
  }
  auto it = selector.selectors.cbegin();
  out << (**it);
  ++it;
  while (it != selector.selectors.cend()) {
    out << ',' << (**it);
    ++it;
  }
  return out;
}

/** @brief Returns true if and only if all the contained selectors return true for a given atom.
 */
class LogicalANDSelector  : public AtomSelector {
public:

  /** @brief Creates a selector combining several AtomSelector objects using logical AND operation
   *
   * Underlying selector may be added by add_selector(AtomSelector_SP) method
   */
  LogicalANDSelector(const std::vector<AtomSelector_SP> & selectors) : selectors(selectors) {}

  /** @brief Creates an empty selector which selects nothing.
   *
   * Underlying selector may be added by add_selector(AtomSelector_SP) method
   */
  LogicalANDSelector() {}

  /** @brief Adds an additional selector to this selection
   *
   * @param selector - a pointer to an AtomSelector object
   */
  void add_selector(AtomSelector_SP selector) { selectors.push_back(selector); }

  /** @brief Returns true if and only if all the contained selectors return true for a given atom.
   *
   * @param c - atom to be tested
   * @return true if all selectors selects atom <code>c</code>
   */
  virtual bool operator()(const PdbAtom & c) const;

private:
  std::vector<AtomSelector_SP> selectors;
};

}
}
}

#endif
/**
 * \example ex_AtomSelector.cc
 * \example ex_SelectResidueRange.cc
 * \example ex_StructureSelector.cc
 */
