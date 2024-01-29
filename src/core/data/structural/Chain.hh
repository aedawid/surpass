/** @file Chain.hh Provides Chain class and some related utility functions
 */
#ifndef CORE_DATA_STRUCTURAL_Chain_H
#define CORE_DATA_STRUCTURAL_Chain_H

#include <string>
#include <memory>
#include <iterator>

#include <core/real.hh>
#include <core/index.hh>

#include <core/data/basic/Vec3.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/structural/Chain.fwd.hh>
#include <core/data/structural/Residue.fwd.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/PdbAtom.fwd.hh>
#include <core/data/structural/Structure.fwd.hh>
#include <core/data/structural/structure_selectors.fwd.hh>
#include <core/data/sequence/SecondaryStructure.hh>

#include <core/chemical/Monomer.hh>

namespace core {
namespace data {
namespace structural {


/** \brief Single biopolymer chain, either a protein or nucleic.
 *
 * This Chain class is derived from <code>std::vector<std::shared_ptr<Residue> > </code>, which allows one
 * use algorithms defined in STL. For example, one can easily remove from a chain these residues which do
 * not have the <code>CA</code> atom:
 *
 * <code>
 * chain.erase(std::remove_if(chain.begin(),chain.end(),[](Residue_SP r){ return r->find_atom(" CA ")==nullptr;}),chain.end());
 * </code>
 */
class Chain :  public std::enable_shared_from_this<Chain>,public std::vector<std::shared_ptr<Residue> > {

  friend std::ostream& operator<<(std::ostream &out, const Chain &r);
public:

  class atom_const_iterator {
public:

    typedef std::vector<std::shared_ptr<Residue> >::const_iterator outer_iterator;
    typedef std::vector<std::shared_ptr<PdbAtom> >::const_iterator inner_iterator;

    typedef std::forward_iterator_tag iterator_category;
    typedef inner_iterator::value_type value_type;
    typedef inner_iterator::difference_type difference_type;
    typedef inner_iterator::pointer pointer;
    typedef inner_iterator::reference reference;

    atom_const_iterator() { }
    atom_const_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
    atom_const_iterator(outer_iterator it, outer_iterator end) : outer_it_(it), outer_end_(end) {

      if (outer_it_ == outer_end_) return; inner_it_ = (*outer_it_)->begin();
      advance_past_empty_inner_containers();
    }

    reference operator*()  const { return *inner_it_;  }
    pointer operator->() const { return &*inner_it_; }

    atom_const_iterator& operator++ () {
      ++inner_it_;
      if (inner_it_ == (*outer_it_)->end()) advance_past_empty_inner_containers();

      return *this;
    }

    atom_const_iterator operator++ (int) {
      atom_const_iterator it(*this);
      ++*this;

      return it;
    }

    friend bool operator==(const atom_const_iterator &a, const atom_const_iterator &b) {
      if (a.outer_it_ != b.outer_it_) return false;

      if (a.outer_it_ != a.outer_end_ && b.outer_it_ != b.outer_end_ && a.inner_it_ != b.inner_it_)
        return false;

      return true;
    }

    friend bool operator!=(const atom_const_iterator &a,const atom_const_iterator &b) { return !(a == b); }

private:

    void advance_past_empty_inner_containers() {
      while (outer_it_ != outer_end_ && inner_it_ == (*outer_it_)->end()) {
        ++outer_it_;
        if (outer_it_ != outer_end_) inner_it_ = (*outer_it_)->begin();
      }
    }

    outer_iterator outer_it_;
    outer_iterator outer_end_;
    inner_iterator inner_it_;
  };  // ~ atom_const_iterator

  class atom_iterator {
public:

    typedef std::vector<std::shared_ptr<Residue> >::iterator outer_iterator;
    typedef std::vector<std::shared_ptr<PdbAtom> >::iterator inner_iterator;

    typedef std::forward_iterator_tag iterator_category;
    typedef inner_iterator::value_type value_type;
    typedef inner_iterator::difference_type difference_type;
    typedef inner_iterator::pointer pointer;
    typedef inner_iterator::reference reference;

    atom_iterator() { }
    atom_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
    atom_iterator(outer_iterator it, outer_iterator end) : outer_it_(it), outer_end_(end) {

      if (outer_it_ == outer_end_) return; inner_it_ = (*outer_it_)->begin();
      advance_past_empty_inner_containers();
    }

    reference operator*()  const { return *inner_it_;  }
    pointer operator->() const { return &*inner_it_; }

    atom_iterator& operator+= (const unsigned int rhs) {

      for (unsigned int i = 0; i < rhs; ++i) {
        ++inner_it_;
        if (inner_it_ == (*outer_it_)->end()) advance_past_empty_inner_containers();
      }
      return *this;
    }

    atom_iterator& operator++ () {
      ++inner_it_;
      if (inner_it_ == (*outer_it_)->end()) advance_past_empty_inner_containers();

      return *this;
    }

    atom_iterator operator++ (int) {
      atom_iterator it(*this);
      ++*this;

      return it;
    }

    friend bool operator==(const atom_iterator &a, const atom_iterator &b) {
      if (a.outer_it_ != b.outer_it_) return false;

      if (a.outer_it_ != a.outer_end_ && b.outer_it_ != b.outer_end_ && a.inner_it_ != b.inner_it_)
    	  return false;

      return true;
    }

    friend bool operator!=(const atom_iterator &a,const atom_iterator &b) { return !(a == b); }

private:

    void advance_past_empty_inner_containers() {
      while (outer_it_ != outer_end_ && inner_it_ == (*outer_it_)->end()) {
        ++outer_it_;
        if (outer_it_ != outer_end_) inner_it_ = (*outer_it_)->begin();
      }
    }

    outer_iterator outer_it_;
    outer_iterator outer_end_;
    inner_iterator inner_it_;
  };  // ~ atom_iterator

  // ---------- C-tors ----------
  /// Constructor creates  a new residue of a given type
  Chain(const char id) : id_(id){}

  /** @brief Creates a new chains as a deep copy of this object.
   *
   * This method makes also a deep copy of residues that belong to this structure
   * (with their atoms, accordingly) providing that they satisfy a given selector.
   *
   * @param which_residues - clone also residues and atoms which satisfy the given selector
   *
   * The following example shows how to clone a fragment of a structure or chain:
   * \include ex_StructureSelector.cc
   */
  Chain_SP clone(const ResidueSelector &which_residues) const;

  // ---------- Getters ----------
  /// Returns a single character identifying this residue, e.g. 'A'
  inline char id() const { return id_; }

  // ---------- Setters ----------
  /// Sets the new character identifying this residue, e.g. 'A'
  inline void id(const char new_id)  { id_ = new_id; }

  // ---------- Atom tree operations ----------
  /// Returns a const-pointer to the structure owning this residue
  const std::shared_ptr<Structure> owner() const {
    std::shared_ptr<Structure> r = owner_.lock();
    if (r) return r;
    else return nullptr;
  }

  /// Sets the new owner (i.e. a structure) that owns this chain
  void owner(std::shared_ptr<Structure> new_owner) { owner_ = new_owner; }

  /// begin() iterator for atoms
  inline Chain::atom_const_iterator  first_const_atom() const { return atom_const_iterator(begin(),end()); };

  /// end() iterator for atoms
  inline Chain::atom_const_iterator  last_const_atom() const { return atom_const_iterator(end()); };

  /// begin() iterator for atoms
  inline Chain::atom_iterator  first_atom() { return atom_iterator(begin(),end()); };

  /// end() iterator for atoms
  inline Chain::atom_iterator  last_atom() { return atom_iterator(end()); };

  /// end() iterator for atoms
  void push_back(Residue_SP r);

  /// Returns residue for a given residue_id number
  Residue_SP get_residue(const core::index2 residue_index);

  // ---------- Misc tree operations ----------
  /// Returns the number of atoms that belong to this chain
  inline core::index4 count_atoms() const {
    core::index4 sum=0;
    for(std::vector<std::shared_ptr<Residue> >::const_iterator it=begin(); it!=end(); ++it) sum +=(*it)->count_atoms();

    return sum;
  }

  /** @brief sorts atoms, residues and chains of this chain according to the PDB ordering.
   *
   * Atoms are sorted by their <code>id</code>. Residues are sorted by the means of
   * <code>bool operator<(const Residue & ri, const Residue & rj)</code> operator
   */
  void sort();

  /// Returns the number of residues that belong to this chain
  inline core::index2 count_residues() const { return size(); }

  /** @brief Creates a Sequence object for this Chain.
   * The index of the first residue of the returned sequence will be set according to the index
   * of the very first residue in this chain. Sequence's comment will hold code of the structure this chain belongs to
   * @return a shared pointer to the newly created Sequence object
   */
  core::data::sequence::SecondaryStructure_SP create_sequence() const;

  /** @brief Counts how many amino acid residues belong to this chain
   * @return the number of amino acid residues in this chain. In the case of a nucleic chain (e.g. a piece of DNA)
   * this method returns 0
   */
  core::index2 count_aa_residues() const;

  /** @brief Counts how many nucleic acid residues belong to this chain
   * @return the number of nucleic acid residues in this chain. In the case of a protein chain
   * this method returns 0
   */
  core::index2 count_na_residues() const;

  friend Structure;

  /** @brief Creates a C-\f$\alpha\f$ - only chain for a given amino acid sequence.
   * @param sequence - requested protein sequence
   * @return a shared pointer to the newly created chain object
   */
  static Chain_SP create_ca_chain(const std::string & sequence, const char chain_id = 'A');

private:
  char id_;             ///< Id for this chain (a single character - according to the PDB convention
  std::weak_ptr<Structure> owner_;
};

/** \brief Copy coordinates of all the atoms of the given chain into a given Coordinates object (which in fact is  std::vector<Vec3>).
 *
 * This method does not check if the destination vector has sufficient size.
 * @param chain - the source of atoms
 * @param coordinates - destination vector
 */
size_t chain_to_coordinates(const Chain_SP structure, core::data::basic::Coordinates & coordinates);

/** \brief Copy coordinates of the matching atoms of the given chain into a given Coordinates object (which in fact is  std::vector<Vec3>).
 *
 * This method does not check if the destination vector has sufficient size.
 * @param chain - the source of atoms
 * @param coordinates - destination vector
 * @param op - atom selector (boolean operator)
 */
size_t chain_to_coordinates(const Chain_SP chain, core::data::basic::Coordinates & coordinates, AtomSelector_SP op);

/// ostream operator just prints the chain code
std::ostream& operator<<(std::ostream &out, const Chain &c  );

}
}
}

#endif
/**
 * \example ex_StructureSelector.cc
 */
