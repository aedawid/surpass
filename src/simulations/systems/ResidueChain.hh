#ifndef SIMULATIONS_SYSTEMS_ResidueChain_HH
#define SIMULATIONS_SYSTEMS_ResidueChain_HH

#include <string>
#include <fstream>
#include <memory>
#include <random>
#include <stdexcept>
#include <type_traits>

#include <core/real.hh>
#include <simulations/atom_indexing.hh>
#include <simulations/systems/AtomTypingInterface.hh>
#include <simulations/systems/AtomTypingVariants.hh>
#include <core/data/structural/Chain.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/PdbAtom.hh>

#include <utils/string_utils.hh>
#include <utils/Logger.hh>

#include <simulations/systems/CartesianAtomsSimple.hh>

namespace simulations {
namespace systems {


template<typename T>
class ResidueChain;

/** @brief A range of atoms (both side inclusive) that belong to a residue or a chain.
 *
 * Because this structure provides iterators, it is possible to use this range in STL algorithms or to iterate over it, as in the following example:
 */
template<class C>
struct AtomRange {
  atom_index first_atom; ///< First atom of this range
  atom_index last_atom; ///< Last atom of this range

  /** @brief  Constructor creates a new range.
   *
   * @param ref_system - system where the atoms belong to
   * @param first_atom - the first atom in this range
   * @param last_atom - the last atom in this range
   */
  AtomRange(ResidueChain<C> & ref_system, const atom_index first_atom, const atom_index last_atom) :
    first_atom(first_atom), last_atom(last_atom), my_atoms(ref_system) {}

  /// Returns begin iterator for this range of atoms
  inline C* begin() const { return my_atoms.coordinates.get() + first_atom; }

  /// Returns pass-the-end iterator for this range of atoms
  inline C* end() const { return my_atoms.coordinates.get() + last_atom + 1; }

  /// Counts atoms that belong to this range
  inline atom_index size() const { return last_atom + 1 - first_atom; }

private:
  ResidueChain<C> & my_atoms;
};

/** @brief Atomic system in the Cartesian space composed of atoms, residues and chains.
 *
 * \include ex_ResidueChain.cc
 * @tparam C - the type used to express coordinates; Vec3 for simple models, Vec3Cubic for modeling in periodic boundary conditions
 */
template<class C>
class ResidueChain : public CartesianAtomsSimple<C> {
public:

  /** @brief Creates a new object based on a biomolecular structure
   *
   * @param s - a Structure object, typically read from a PDB file
   */
  ResidueChain(const AtomTypingInterface_SP typing, const core::data::structural::Structure & s);

  /** @brief Creates a system of of <code>n_atoms</code> residues, each of them comprising one atom, placed in a single chain
   *
   * @param typing - object used to translate between PDB atom names and force field atom types
   * @param n_atoms - the total number of atoms in this system
   * @param atom_name - string representing a PDB-style atom name (the same for all atoms)
   */
  ResidueChain(const AtomTypingInterface_SP typing, const core::index4 n_atoms, const std::string atom_name = " CA ");

  /** @brief Creates a system of of multiple chains
   *
   * @param typing - object used to translate between PDB atom names and force field atom types
   * @param n_atoms_in_chains - the number of residues (one atom each) in each chain
   * @param atom_name - string representing a PDB-style atom name (the same for all atoms)
   */
  ResidueChain(const AtomTypingInterface_SP typing, const std::vector<core::index2> & n_atoms_in_chains,
               const std::string atom_name = " CA ");

  /// Necessary virtual destructor
  virtual ~ResidueChain() {}

  /** @brief Returns the number of atoms in a given residue.
   *
   * The total number of atoms in this system may be checked by calling <code>CartesianAtomsSimple<C>::count_atoms();</code>
   * @param residue_index - which residue?
   */
  inline core::index2 count_atoms(const core::index2 residue_index) const { return atoms_for_residue_[residue_index].last_atom - atoms_for_residue_[residue_index].first_atom + 1; }

  /** @brief look for an atom by its internal name.
   *
   * Attempts to find an atom called <code>atom_name</code> in a residue <code>i_res</code>. Throws an exception
   * if nothing was found.
   *
   * @param i_res - residue index
   * @param internal_atom_name - name of an atom to look for
   * @return a reference to the atom of interest
   * @throw std::range_error if the requested atom cannot be found
   */
  C & find_atom(const residue_index i_res, const std::string & internal_atom_name) {

    // --- we get the first atom from the residue just to read the type of that residue
    core::index2 residue_type = CartesianAtomsSimple<C>::coordinates[atoms_for_residue_[i_res].first_atom].residue_type;
    // --- here we obtain the type ot the atom, according to the force field (i.e. the relevant AtomTyping)
    core::index2 type = CartesianAtomsSimple<C>::atom_typing->atom_type(
      internal_atom_name,core::chemical::Monomer::get(residue_type).code3,AtomTypingVariants::STANDARD);

    for (auto &a : atoms_for_residue_[i_res]) if (a.atom_type == type) return a;

    throw new std::range_error(utils::string_format("Can't find atom %s  in residue %d",internal_atom_name.c_str(),i_res));
  }

  /** @brief look for an atom by its internal name.
   *
   * Attempts to find an atom called <code>internal_atom_name</code> in a residue <code>i_res</code>.
   *
   * @param i_res - residue index
   * @param internal_atom_name - name of an atom to look for. The name should follow internal atom names used
   * by the force field (i.e. AtomTyping instance) rather than just PDB codes
   * @return index of a requested atoms or <code>std::numeric_limits<core::index4>::max()</code> when not found
   */
  core::index4 find_atom_index(const core::index2 i_res, const std::string & internal_atom_name) const {

    // --- we get the first atom from the residue just to read the type of that residue
    core::index2 residue_type = CartesianAtomsSimple<C>::coordinates[atoms_for_residue_[i_res].first_atom].residue_type;
    // --- here we obtain the type ot the atom, according to the force field (i.e. the relevant AtomTyping)
    core::index2 type = CartesianAtomsSimple<C>::atom_typing->atom_type(
      internal_atom_name, core::chemical::Monomer::get(residue_type).code3,AtomTypingVariants::STANDARD);
    // --- and finally we iterate over all atoms from the residue and look for the one of the requested type.
    for (atom_index i = atoms_for_residue_[i_res].first_atom; i <= atoms_for_residue_[i_res].last_atom; ++i) {
      if (CartesianAtomsSimple<C>::coordinates[i].atom_type == type) return i;
    }

    return std::numeric_limits<core::index4>::max();
  }

  /** @brief look for an atom by the order number of a residue this atom belongs to and its atom_type.
   *
   * Attempts to find an atom called <code>atom_type</code> in a residue <code>i_res</code>.
   *
   * @param i_res - residue index
   * @param atom_type - atom typing of an atom to look for
   * @return index of a requested atoms or <code>std::numeric_limits<core::index4>::max()</code> when not found
   */
  core::index4 find_atom_index(const core::index2 i_res, const core::index2 & type) const {

    for (atom_index i = atoms_for_residue_[i_res].first_atom; i <= atoms_for_residue_[i_res].last_atom; ++i) {
      if (CartesianAtomsSimple<C>::coordinates[i].atom_type == type) return i;
    }

    return std::numeric_limits<core::index4>::max();
  }

  /** @brief Returns the PDB name of a given atom.
   *
   * @param atom_index - which atom?
   * @return PDB-style name of the atom; always four characters, e.g. " CA "
   */
  const std::string & pdb_atom_name(const core::index4 atom_index) const { return pdb_atom_names_[atom_index]; }

  /** @brief Exposes the vector of all PDB-style atom names.
   *
   * @return a const-reference to the vector of atom names
   */
  const std::vector<std::string> & pdb_atom_names() const { return pdb_atom_names_; }

  /** @brief Returns the type of the residue the given atom belongs to.
   *
   * This method simply returns Vec3::residue_type value for the relevant atom.
   * @param atom_index index of the atom of interest
   */
  inline core::index1 residue_type(const core::index4 atom_index) const { return CartesianAtomsSimple<C>::coordinates[atom_index].residue_type; }

  /// Counts residues of this system
  inline residue_index count_residues() const { return atoms_for_residue_.size(); }

  /// Counts chains of this system
  inline residue_index count_chains() const { return atoms_for_chain_.size(); }

  /// Returns an atom range for a given residue
  inline const AtomRange<C> & atoms_for_residue(const simulations::residue_index residue_index) const {
    return atoms_for_residue_[residue_index];
  }

  /// Returns the index of the residue a given atom belongs to
  inline core::index2 residue_for_atom(const core::index4 atom_index) const {

    return CartesianAtomsSimple<C>::coordinates[atom_index].residue_index;
  }

  /// Returns an atom range for a given chain
  inline const AtomRange<C> & atoms_for_chain(const simulations::residue_index chain_index) const {
    return atoms_for_chain_[chain_index];
  }

  /// Returns a chain index a given atom index
  inline core::index2 chain_for_atom(const core::index4 which_atom) const {

    for (core::index2 i = 0; i < count_chains(); ++i) {
      if ((which_atom >= atoms_for_chain_[i].first_atom) && (which_atom <= atoms_for_chain_[i].last_atom)) return i;
    }
    throw std::out_of_range(utils::string_format(
      "Can't find chain index for atom %d. The system has only %d atom in %d residues, %d chains\n",
      which_atom, CartesianAtomsSimple<C>::n_atoms, count_residues(), count_chains()));
  }

  /// Returns integer ID of this system
  inline core::index2 system_id() const { return system_id_; }

  /// Sets a new integer ID for this system
  inline void system_id(core::index2 id) { system_id_ = id; }

  virtual void write_pdb(std::ostream & where, const std::vector<std::string> & atom_format_lines,
                          core::index4 model_id = 0) const {

    if (model_id > 0) where << utils::string_format("MODEL %6d\n", model_id);
    core::index2 nc = count_chains();
    if (nc > 1) { // --- iterate over chains - periodic boundary conditions
      C o(0.0, 0.0, 0.0); // --- temporary vector to keep the origin of the current chain
      int i_line = -1;
      for (core::index2 ic = 0; ic < nc; ++ic) {
        CartesianAtomsSimple<C>::coordinates[atoms_for_chain(ic).last_atom].wrap(o);
        o -= CartesianAtomsSimple<C>::coordinates[atoms_for_chain(ic).last_atom];
        for (core::index4 ia = atoms_for_chain(ic).first_atom; ia <= atoms_for_chain(ic).last_atom; ++ia) {
          C & at = CartesianAtomsSimple<C>::coordinates[ia];
          where << utils::string_format(atom_format_lines[++i_line], at.x + o.x, at.y + o.y, at.z + o.z);
        }
      }
    } else {
      core::data::basic::Vec3 cm;
      for (core::index4 i = 0; i < CartesianAtomsSimple<C>::count_atoms(); i++)
        cm += CartesianAtomsSimple<C>::coordinates[i];
      cm /= double(CartesianAtomsSimple<C>::count_atoms());
      core::data::io::TVect tv(1,cm.x,cm.y,cm.z);
      where << tv.to_pdb_line() << "\n";
      int i_line = -1;
      for (core::index4 i = 0; i < CartesianAtomsSimple<C>::count_atoms(); i++) {
        C & ic = CartesianAtomsSimple<C>::coordinates[i];
        where << utils::string_format(atom_format_lines[++i_line], ic.x-cm.x, ic.y-cm.y, ic.z-cm.z);
      }
    }
    if (model_id > 0) where << "ENDMDL\n";
  }

  /** @brief A helper static method that initializes system's coordinates as a single homopolymer chain.
   *
   * This method simply calls BuildPolymerChain<C>::generate() method.
   * In brief, this method does not create any new atoms. It takes an initialized system and sets its coordinates so they form a chain.
   * Atoms of the resulting chain are placed <code>bond_length</code> Angstroms apart. The atoms satisfy excluded
   * volume i.e. they are not closer to each other than <code>cutoff</code>. This method makes <code>n_bead_attempts</code>
   * attempts to place each bead (atom) so it does not violate excluded volume. After that a chain is restarted from
   * beginning. After unsuccessful <code>n_chain_attempts</code> chain attempts the method returns false.
   */
  static bool homopolymer(CartesianAtomsSimple<C> & system, const core::real bond_length, const core::real cutoff,
      const core::index2 n_chain_attempts = 1000, const core::index2 n_bead_attempts = 1000) {

    BuildPolymerChain<C> builder(system.coordinates,system.n_atoms);
    return builder.generate(3.8,4.0);
  }

protected:
  std::vector<AtomRange<C>> atoms_for_residue_;
  std::vector<AtomRange<C>> atoms_for_chain_;
  std::vector<std::string> pdb_atom_names_;
private:
  core::index2 system_id_;
  utils::Logger logger;
};

template<class C>
ResidueChain<C>::ResidueChain(const AtomTypingInterface_SP typing, const core::index4 n_atoms,
    const std::string atom_name) : CartesianAtomsSimple<C>(typing, n_atoms), system_id_(0), logger("ResidueChain")  {

  for (residue_index i = 0; i < CartesianAtomsSimple<C>::n_atoms; i++) {
    atoms_for_residue_.emplace_back(*this, i, i);
    pdb_atom_names_.push_back(atom_name);
  }
  atoms_for_chain_.emplace_back(*this, 0, CartesianAtomsSimple<C>::n_atoms - 1);
}

template<class C>
ResidueChain<C>::ResidueChain(const AtomTypingInterface_SP typing, const core::data::structural::Structure & s) :
    ResidueChain<C>(typing, s.count_atoms()) {

  core::index4 first_atom_in_next_residue = 0;
  core::index4 first_atom_in_next_chain = 0;
  core::index4 i_atom = 0;
  core::index2 i_resid = 0;
  atoms_for_chain_.clear();
  atoms_for_residue_.clear();
  pdb_atom_names_.clear();

  for(const auto & cp : s) {
    logger << utils::LogLevel::FINE << "processing chain " << (*cp).id() << "\n";
    for(const auto & rp : *cp) {
      for(const auto & ap : *rp) {
        C & a = CartesianAtomsSimple<C>::coordinates[i_atom];
        a.set(*ap);
        a.chain_id = cp->id();
        pdb_atom_names_.push_back(ap->atom_name());
        try {
          a.atom_type = typing->atom_type(*ap);
        } catch (std::out_of_range &e) {
          logger << utils::LogLevel::SEVERE << "Unknown atom >" << ap->atom_name() << "< found in residue " << rp->residue_type().code3 << "\n";
        }
        a.residue_index = i_resid;
        a.residue_type = rp->residue_type().parent_id;
        ++i_atom;
      }
        ++i_resid;
      atoms_for_residue_.emplace_back(*this, first_atom_in_next_residue,i_atom-1);

      // --- The IF below is just to print detailed message what atoms go to the current residue and what's their typing
      if(logger.is_logable(utils::LogLevel::FINER)) {
        C & a = CartesianAtomsSimple<C>::coordinates[atoms_for_residue_.back().last_atom];
        logger << utils::LogLevel::FINER << "Atoms for residue " << size_t(a.residue_index) << " "
            << core::chemical::Monomer::get(core::index2(a.residue_type)).code3 << " : " <<
            size_t(atoms_for_residue_.back().first_atom) << " - " << size_t(atoms_for_residue_.back().last_atom) << "\n";
        for (atom_index ai = atoms_for_residue_.back().first_atom; ai <= atoms_for_residue_.back().last_atom; ++ai) {
          C & aa = CartesianAtomsSimple<C>::coordinates[ai];
          logger << utils::LogLevel::FINER << "\t" << size_t(ai)  << " is typed "<< aa.atom_type
            <<" "<<CartesianAtomsSimple<C>::atom_typing->atom_internal_name(aa.atom_type)<< "\n";
        }
      }
      first_atom_in_next_residue = i_atom;
    }
    atoms_for_chain_.emplace_back(*this, first_atom_in_next_chain,i_atom-1);
    first_atom_in_next_chain = i_atom;
  }

  for (core::index2 iresid = 0; iresid < count_residues(); ++iresid) {
    for (core::index4 ia = atoms_for_residue_[iresid].first_atom; ia <= atoms_for_residue_[iresid].last_atom; ++ia)
      if (CartesianAtomsSimple<C>::operator[](ia).residue_index != iresid)
        throw std::invalid_argument(utils::string_format("Inconsistent residue indexing for atom %d\n",ia));
  }

  for (core::index2 ichain = 0; ichain < count_chains(); ++ichain) {
    for (core::index4 ia = atoms_for_chain_[ichain].first_atom; ia <= atoms_for_chain_[ichain].last_atom; ++ia)
      if (chain_for_atom(ia) != ichain)
        throw std::invalid_argument(utils::string_format("Inconsistent chain indexing for atom %d\n\tChain cached: %d chain assumed: %d\n"
          ,ia, chain_for_atom(ia), ichain));
  }
}

template<class C>
ResidueChain<C>::ResidueChain(const AtomTypingInterface_SP typing, const std::vector<core::index2> & n_atoms_in_chains,
    const std::string atom_name) : ResidueChain<C>(typing, std::accumulate(n_atoms_in_chains.begin(),
    n_atoms_in_chains.end(), 0),atom_name) {

  atom_index i_atom = 0;
  char chain_char = 'A';
  atoms_for_chain_.clear();
  atoms_for_residue_.clear();
  for (core::index2 n : n_atoms_in_chains) { // --- loop over chains
    residue_index i_res = 0;
    atoms_for_chain_.emplace_back(*this, i_atom, i_atom + n - 1);
    for (atom_index i = i_atom; i < i_atom + n; ++i) { // --- loop over atoms in n-th chain
      CartesianAtomsSimple<C>::coordinates[i].chain_id = chain_char;
      CartesianAtomsSimple<C>::coordinates[i].residue_index = i_res;
      ++i_res;
    }
    i_atom += n;
    ++chain_char;
  }
  for (core::index2 iresid = 0; iresid < count_residues(); ++iresid) {
    for (core::index4 ia = atoms_for_residue_[iresid].first_atom; ia <= atoms_for_residue_[iresid].last_atom; ++ia)
      if (CartesianAtomsSimple<C>::operator[](ia).residue_index != iresid)
        throw std::invalid_argument(utils::string_format("Inconsistent residue indexing for atom %d\n",ia));
  }

  for (core::index2 ichain = 0; ichain < count_chains(); ++ichain) {
    for (core::index4 ia = atoms_for_chain_[ichain].first_atom; ia <= atoms_for_chain_[ichain].last_atom; ++ia)
      if (chain_for_atom(ia) != ichain)
        throw std::invalid_argument(utils::string_format("Inconsistent chain indexing for atom %d\n",ia));
  }
}

} // ~ simulations
} // ~ cartesian

#endif
/**
 * \example ex_ResidueChain.cc
 */
