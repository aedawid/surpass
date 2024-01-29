#ifndef CORE_DATA_STRUCTURAL_Structure_H
#define CORE_DATA_STRUCTURAL_Structure_H

#include <stdexcept>

#include <core/index.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/io/PdbField.hh>
#include <core/data/structural/Chain.hh>
#include <core/data/structural/PdbAtom.hh>
#include <core/data/structural/Structure.fwd.hh>
#include <core/data/structural/structure_selectors.fwd.hh>

namespace core {
namespace data {
namespace structural {

/** \brief Represents a biomacromolecular structure.
 *
 * To read a Structure from a PDB file, use Pdb class.
 * Structures may be conveniently loaded at the command line level from a PDB file by using
 * <code>utils::options::structures_from_cmdline()</code> declared in input_utils.hh
 *
 * To write a structure to a PDB-formatted file, use <code>core::data::io::write_pdb()</code> declared in Pdb.hh file
 *
 * The following example reads a PDB file and creates a Structure object
 * \include ex_Pdb.cc
 *
 * The following example iterates over chains, residues and atom; it prints all atoms of a residue in a single line, one line per residue
 * \include ex_Structure.cc
 */
class Structure :  public std::enable_shared_from_this<Structure>, public std::vector<std::shared_ptr<Chain> >{
public:
  class atom_const_iterator {
  public:

      typedef std::vector<std::shared_ptr<Chain> >::const_iterator outer_iterator;
      typedef Chain::atom_iterator inner_iterator;

      typedef std::forward_iterator_tag iterator_category;
      typedef inner_iterator::value_type value_type;
      typedef inner_iterator::difference_type difference_type;
      typedef inner_iterator::pointer pointer;
      typedef inner_iterator::reference reference;

      atom_const_iterator() { }
      atom_const_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
      atom_const_iterator(outer_iterator it, outer_iterator end) : outer_it_(it), outer_end_(end) {

        if (outer_it_ == outer_end_) return; inner_it_ = (*outer_it_)->first_atom();
        advance_past_empty_inner_containers();
      }

      reference operator*()  const { return *inner_it_;  }
      pointer operator->() const { return &*inner_it_; }

      atom_const_iterator& operator++ () {
        ++inner_it_;
        if (inner_it_ == (*outer_it_)->last_atom()) advance_past_empty_inner_containers();

        return *this;
      }

      atom_const_iterator operator++ (int) {
        atom_const_iterator it(*this);
        ++*this;

        return it;
      }

      atom_const_iterator& operator+=(const unsigned int rhs){

        for(unsigned int i=0;i<rhs;++i)  ++*this;
        return *this;
      }

      bool operator==(const atom_const_iterator &b) const {
        if (outer_it_ != b.outer_it_) return false;
        if (outer_it_ != outer_end_ && b.outer_it_ != b.outer_end_ && inner_it_ != b.inner_it_)
          return false;

        return true;
      }

      bool operator!=(const atom_const_iterator &b) const { return !operator==(b); }

  private:

      void advance_past_empty_inner_containers() {
        while (outer_it_ != outer_end_ && inner_it_ == (*outer_it_)->last_atom()) {
          ++outer_it_;
          if (outer_it_ != outer_end_) inner_it_ = (*outer_it_)->first_atom();
        }
      }

      outer_iterator outer_it_;
      outer_iterator outer_end_;
      inner_iterator inner_it_;
    }; // ~ atom_iterator

  class atom_iterator {
	public:

	    typedef std::vector<std::shared_ptr<Chain> >::iterator outer_iterator;
	    typedef Chain::atom_iterator inner_iterator;

	    typedef std::forward_iterator_tag iterator_category;
	    typedef inner_iterator::value_type value_type;
	    typedef inner_iterator::difference_type difference_type;
	    typedef inner_iterator::pointer pointer;
	    typedef inner_iterator::reference reference;

	    atom_iterator() { }
	    atom_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
	    atom_iterator(outer_iterator it, outer_iterator end) : outer_it_(it), outer_end_(end) {

	      if (outer_it_ == outer_end_) return; inner_it_ = (*outer_it_)->first_atom();
	      advance_past_empty_inner_containers();
	    }

	    reference operator*()  const { return *inner_it_;  }
	    pointer operator->() const { return &*inner_it_; }

    atom_iterator& operator+= (const unsigned int rhs) {

        for(unsigned int i=0;i<rhs;++i)  ++*this;
        return *this;
      }


    atom_iterator& operator++ () {
	      ++inner_it_;
	      if (inner_it_ == (*outer_it_)->last_atom()) advance_past_empty_inner_containers();

	      return *this;
	    }

	    atom_iterator operator++ (int) {
	      atom_iterator it(*this);
	      ++*this;

	      return it;
	    }

        bool operator==(const atom_iterator &b) const {
          if (outer_it_ != b.outer_it_) return false;
          if (outer_it_ != outer_end_ && b.outer_it_ != b.outer_end_ && inner_it_ != b.inner_it_)
            return false;

          return true;
        }

        bool operator!=(const atom_iterator &b) const { return !operator==(b); }

  private:

	    void advance_past_empty_inner_containers() {
	      while (outer_it_ != outer_end_ && inner_it_ == (*outer_it_)->last_atom()) {
	        ++outer_it_;
	        if (outer_it_ != outer_end_) inner_it_ = (*outer_it_)->first_atom();
	      }
	    }

	    outer_iterator outer_it_;
	    outer_iterator outer_end_;
	    inner_iterator inner_it_;
	  }; // ~ atom_iterator

  class residue_const_iterator {
  public:

      typedef std::vector<std::shared_ptr<Chain> >::const_iterator outer_iterator;
      typedef std::vector<std::shared_ptr<Residue> >::const_iterator inner_iterator;

      typedef std::forward_iterator_tag iterator_category;
      typedef inner_iterator::value_type value_type;
      typedef inner_iterator::difference_type difference_type;
      typedef inner_iterator::pointer pointer;
      typedef inner_iterator::reference reference;

      residue_const_iterator() { }
      residue_const_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
      residue_const_iterator(outer_iterator it, outer_iterator end) : outer_it_(it), outer_end_(end) {

        if (outer_it_ == outer_end_) return; inner_it_ = (*outer_it_)->begin();
        advance_past_empty_inner_containers();
      }

      reference operator*()  const { return *inner_it_;  }
      pointer operator->() const { return &*inner_it_; }

      residue_const_iterator& operator++ () {
        ++inner_it_;
        if (inner_it_ == (*outer_it_)->cend()) advance_past_empty_inner_containers();

        return *this;
      }

      residue_const_iterator operator++ (int) {
        residue_const_iterator it(*this);
        ++*this;

        return it;
      }

    residue_const_iterator& operator+=(const unsigned int rhs){

      for(unsigned int  i=0;i<rhs;++i)  ++*this;
      return *this;
    }

      friend bool operator==(const residue_const_iterator &a, const residue_const_iterator &b) {
        if (a.outer_it_ != b.outer_it_) return false;

        if (a.outer_it_ != a.outer_end_ && b.outer_it_ != b.outer_end_ && a.inner_it_ != b.inner_it_)
          return false;

        return true;
      }

      friend bool operator!=(const residue_const_iterator &a,const residue_const_iterator &b) { return !(a == b); }

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
    }; // ~ residue_const_iterator

  class residue_iterator {
  public:

      typedef std::vector<std::shared_ptr<Chain> >::iterator outer_iterator;
      typedef std::vector<std::shared_ptr<Residue> >::iterator inner_iterator;

      typedef std::forward_iterator_tag iterator_category;
      typedef inner_iterator::value_type value_type;
      typedef inner_iterator::difference_type difference_type;
      typedef inner_iterator::pointer pointer;
      typedef inner_iterator::reference reference;

      residue_iterator() { }
      residue_iterator(outer_iterator it) : outer_it_(it), outer_end_(it) { }
      residue_iterator(outer_iterator it, outer_iterator end) : outer_it_(it), outer_end_(end) {

        if (outer_it_ == outer_end_) return; inner_it_ = (*outer_it_)->begin();
        advance_past_empty_inner_containers();
      }

      reference operator*()  const { return *inner_it_;  }
      pointer operator->() const { return &*inner_it_; }

    friend residue_iterator operator+(residue_iterator lhs, const unsigned int rhs) {
      lhs += rhs;
      return lhs;
    }

    residue_iterator& operator+=(const unsigned int rhs){

      for(unsigned int  i=0;i<rhs;++i)  ++*this;
      return *this;
    }

    residue_iterator& operator++ () {
        ++inner_it_;
        if (inner_it_ == (*outer_it_)->last_atom()) advance_past_empty_inner_containers();

        return *this;
      }

      residue_iterator operator++ (int) {
        residue_iterator it(*this);
        ++*this;

        return it;
      }

      friend bool operator==(const residue_iterator &a, const residue_iterator &b) {
        if (a.outer_it_ != b.outer_it_) return false;

        if (a.outer_it_ != a.outer_end_ && b.outer_it_ != b.outer_end_ && a.inner_it_ != b.inner_it_)
          return false;

        return true;
      }

      friend bool operator!=(const residue_iterator &a,const residue_iterator &b) { return !(a == b); }

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
    }; // ~ residue_iterator

  // public data
  std::multimap<std::string, std::shared_ptr<core::data::io::PdbField>> pdb_header;

  // ---------- C-tor ----------
  /** \brief Creates a structure with no atoms.
   *
   * @param code - a structure ID, preferably 4-character string so it looks like a PDB code
   */
  Structure(const std::string &code) : code_(code), logger("Structure") { }

  /** @brief Creates a new Structure as a deep copy of this object.
   *
   * This method makes also a deep copy of chains  that belong to this structure (with their residues and atoms, accordingly)
   * providing that they satisfy a given selector.
   *
   * @param which_chains - clone also chains, residues and atoms which satisfy the given selector
   */
  Structure_SP clone(const ChainSelector &which_chains) const;

  // ---------- Getters ----------
  inline const std::string & code() const { return code_; }

  /** @brief Returns true if this Structure has a chain identified by a given character.
   * @param code - chainId
   * @return true if the chain was found in this structure; false otherwise
   */
  bool has_chain(const char code) const;

  /** @brief Returns the chain for the given index.
   * Throws an exception, if the chain index is too high
   * @param index - index of the requested chain
   * @return pointer to the relevant chain
   */
  std::shared_ptr<Chain> get_chain(const core::index2 index) { return at(index); }

  /** @brief Returns the const-pointer to a chain for the given index.
   * Throws an exception, if the chainId is invalid
   * @param code - chainId
   * @return pointer to the relevant chain
   */
  const std::shared_ptr<Chain> get_chain(const char code) const;

  /** @brief Returns the chain for the given code.
   * Throws an exception, if the chain code is invalid i.e. there is no such chain in this Structure object
   * @param code - chainId
   * @return pointer to the relevant chain
   */
  std::shared_ptr<Chain> get_chain(const char code);

  /** @brief Returns a pointer to a residue from this Structure
   * Throws an exception, if the chain code is invalid i.e. there is no such chain in this Structure object
   * @param chain_code - chainId the residue belongs to
   * @param residue_id - resId of the residue
   * @param icode - insertion code for the requested residue
   * @return pointer to the relevant residue
   */
  std::shared_ptr<Residue> get_residue(const char chain_code, const core::index2 residue_id, const char icode);

  /** @brief Returns a const-pointer to a residue from this Structure
   * Throws an exception, if the chain code is invalid i.e. there is no such chain in this Structure object
   * @param chain_code - chainId the residue belongs to
   * @param residue_id - resId of the residue
   * @param icode - insertion code for the requested residue
   * @return const-pointer to the relevant residue
   */
  const std::shared_ptr<Residue> get_residue(const char chain_code, const core::index2 residue_id, const char icode) const;

  /** @brief Returns a const-pointer to a residue requested by its index.
   *
   * @param residue_index - index of the residue of interest; It starts from zero and goes continously through all chains
   * @return a const-pointer to the requested residue
   */
  const std::shared_ptr<Residue> get_residue(const core::index2 residue_index) const;

  /** @brief Returns a pointer to a residue requested by its index.
   *
   * @param residue_index - index of the residue of interest; It starts from zero and goes continously through all chains
   * @return a pointer to the requested residue
   */
  std::shared_ptr<Residue> get_residue(const core::index2 residue_index);

  /// Returns the polymer sequence as stored in SEQRES field in the source PDB file
  std::vector<core::chemical::Monomer> & original_sequence(const char chain_code);

  // ---------- Setters ----------
  inline void code(const std::string & new_code) { code_ = new_code; }

	// ---------- Chain tree operations ----------
  /// begin() iterator for residues
	inline Structure::residue_iterator  first_residue() { return residue_iterator(begin(),end()); }

  /// end() iterator for residues
	inline Structure::residue_iterator  last_residue() { return residue_iterator(end()); }

  /// begin() const iterator for residues
  inline Structure::residue_const_iterator  first_const_residue() const { return residue_const_iterator(begin(),end()); }

  /// end() const iterator for residues
  inline Structure::residue_const_iterator  last_const_residue() const { return residue_const_iterator(end()); }

  /// begin() iterator for atoms
  inline Structure::atom_iterator  first_atom() { return atom_iterator(begin(),end()); }

  /// end() const iterator for atoms
  inline Structure::atom_iterator  last_atom() { return atom_iterator(end()); }

  /// begin() const iterator for atoms
  inline Structure::atom_const_iterator  first_const_atom() const { return atom_const_iterator(begin(),end()); }

  /// end() iterator for atoms
  inline Structure::atom_const_iterator  last_const_atom() const { return atom_const_iterator(end()); }

  void push_back(std::shared_ptr<Chain> c);

  // ---------- Misc operations ----------
	/// Returns the number of atoms that belong to this structure
	inline core::index4 count_atoms() const {
		core::index4 sum=0;
		for(std::vector<std::shared_ptr<Chain> >::const_iterator it=begin();it!=end();++it)
			sum +=(*it)->count_atoms();
		return sum;
	}

	/// Returns the number of residues that belong to this structure
	inline core::index2 count_residues() const {
		core::index2 sum=0;
		for(std::vector<std::shared_ptr<Chain> >::const_iterator it=begin();it!=end();++it)
			sum +=(*it)->count_residues();
		return sum;
	}

	/// Returns the number of chains that belong to this structure
	inline core::index2 count_chains() const { return size(); }

  /** @brief Sort chains, residues and atoms of this structure.
   *
   * Chains are sorted lexicographically. Atoms are sorted by their <code>id</code>. Residues are sorted by the means of
   * <code>bool operator<(const Residue & ri, const Residue & rj)</code> operator.
   */
  void sort();

private:
  std::string code_;
  utils::Logger logger;
};

/** \brief Copy coordinates of all atoms of a structure into a std::unique_ptr<C[]> container
 *
 * This method does not check if the destination vector has sufficient size.
 * @param coordinates - source structure
 * @param structure - destination vector of coordinates
 */
template <typename C>
size_t structure_to_coordinates(const Structure_SP structure, std::unique_ptr<C[]>  & coordinates);

/** \brief Copy coordinates of all atoms of a Coordinates object (which in fact is  std::vector<Vec3>) into a given Structure.
 *
 * This method does not check if the destination vector has sufficient size.
 * @param coordinates - source vector of coordinates
 * @param structure - destination structure
 */
size_t coordinates_to_structure(const core::data::basic::Coordinates & coordinates, Structure & structure);

/** \brief Copy coordinates of all atoms of a given container into a given Structure.
 *
 * This method does not check if the destination vector has sufficient size.
 * @param coordinates - source vector of coordinates
 * @param structure - destination structure
 */
template <typename C>
size_t coordinates_to_structure(const std::unique_ptr<C[]>  & coordinates,const Structure_SP structure);

/** \brief Copy coordinates of all atoms of this structure into a given Coordinates object (which in fact is  std::vector<Vec3>).
 *
 * This method does not check if the destination vector has sufficient size.
 * @param structure - the source of atoms
 * @param coordinates - destination vector
 * @return the number of atoms copied
 */
size_t structure_to_coordinates(const Structure_SP structure, core::data::basic::Coordinates & coordinates);

/** \brief Copy coordinates of the matching atoms of this structure into a given Coordinates object (which in fact is  std::vector<Vec3>).
 *
 * This method does not check if the destination vector has sufficient size.
 * @param structure - the source of atoms
 * @param coordinates - destination vector
 * @param op - atom selector (boolean operator) used to pick the atoms of interest. E.g. use core::data::structural::IsCA to copy coordinates of alpha-carbons only.
 * @return the number of atoms copied
 */
size_t structure_to_coordinates(const Structure_SP structure, core::data::basic::Coordinates & coordinates,
    const AtomSelector_SP op);

/** \brief Copy coordinates of the matching atoms of this structure into a given Coordinates object (which in fact is  std::vector<Vec3>).
 *
 * This method does not check if the destination vector has sufficient size.
 * @param structure - the source of atoms
 * @param coordinates - destination vector
 * @param op - atom selector (boolean operator) used to pick the atoms of interest. E.g. use core::data::structural::IsCA to copy coordinates of alpha-carbons only.
 * @return the number of atoms copied
 */
size_t structure_to_coordinates(const Structure_SP structure, core::data::basic::Coordinates & coordinates,
    const AtomSelector & op);
}
}
}

#endif

/**
 * \example ex_Structure.cc
 */
