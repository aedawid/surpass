#ifndef SIMULATIONS_SYSTEMS_CartesianAtomsSimple_HH
#define SIMULATIONS_SYSTEMS_CartesianAtomsSimple_HH

#include <string>
#include <fstream>
#include <memory>
#include <random>

#include <core/real.hh>
#include <core/calc/statistics/Random.hh>
#include <core/calc/structural/angles.hh>
#include <core/data/structural/PdbAtom.hh>

#include <utils/string_utils.hh>
#include <utils/Logger.hh>

#include <simulations/atom_indexing.hh>
#include <simulations/systems/AtomTypingInterface.hh>
#include <simulations/systems/BuildPolymerChain.hh>

namespace simulations {
namespace systems {

/** @brief A simple system of atoms defined in the Cartesian space.
 *
 * This object does not carry information on residues nor chains the atoms belong to. If you need to
 * keep atoms in higher-level object, use derived classes: ResidueChain or ResidueMultiChain
 *
 * @tparam C - object used to represent a single atom, bead, particle, etc. Currently a library
 * provides two classes to be used here: Vec3 and Vec3Cubic.
 */
template<class C>
class CartesianAtomsSimple {
public:
  /// AtomTyping object converts PDB-style atom names into internal numbering used by a force field
  const AtomTypingInterface_SP atom_typing;
  const atom_index n_atoms; ///< The total number of atoms in this system
  std::unique_ptr<C[]> coordinates; ///< Atomic coordinates

  /** @brief Create a new system on <code>n_atoms</code> atoms.
   * @param typing - object used to convert internal atom type indexes to PDB atom names and vice versa
   * @param n_atoms - the number of atoms in this system
   */
  CartesianAtomsSimple(const AtomTypingInterface_SP typing, const size_t n_atoms) :
      atom_typing(typing), n_atoms(n_atoms), coordinates(new C[n_atoms]) {

    for (residue_index i = 0; i < CartesianAtomsSimple<C>::n_atoms; i++) {
      coordinates[i].residue_index = i + 1;
      coordinates[i].atom_type = 0;
      coordinates[i].chain_id = 'A';
    }
  }

  /** @brief Copying constructor.
   *
   * @param source - a system to be duplicated
   */
  CartesianAtomsSimple(const CartesianAtomsSimple& source) :
    atom_typing(source.atom_typing), n_atoms(source.n_atoms), coordinates(new C[n_atoms]) {

    for (core::index4 i = 0; i < n_atoms; ++i) {
      coordinates[i].set(source.coordinates[i]);
      coordinates[i].chain_id = source.coordinates[i].chain_id;
      coordinates[i].atom_type = source.coordinates[i].atom_type;
      coordinates[i].residue_type = source.coordinates[i].residue_type;
      coordinates[i].residue_index = source.coordinates[i].residue_index;
      coordinates[i].register_ = source.coordinates[i].register_;
    }
  }

  /// Necessary virtual destructor
  virtual ~CartesianAtomsSimple() {}

  /** @brief Index-operator returns coordinates of a requested atom
   * @param index - index of an atom; starts from 0 of course
   */
  inline C& operator[](const core::index4 & index) { return coordinates[index]; }

  /** @brief Index-operator returns coordinates of a requested atom
   * @param index - index of an atom; starts from 0 of course
   */
  inline C& operator[](const core::index4 & index) const { return coordinates[index]; }

  /** @brief Counts atoms of this system
   * @return the number of atoms in this system
   */
  inline atom_index count_atoms() const { return n_atoms; }

  /** @brief Returns the box size
   * @return side length of a simulation box the system lives in
   */
  inline core::real box_side() const { return coordinates[0].box_side; }

  /** @brief Writes the state of this system in PDB format.
   * @param where - destination stream
   * @param model_id - if greater than 0, this method will print MODEL / ENDMDL lines in the PDB output. By default
   *    model_id = 0, so the lines are not printed
   */
  virtual void write_pdb(std::ostream & where, const std::vector<std::string> & atom_format_lines,
                         const core::index4 model_id = 0) const {

#ifdef DEBUG
  if (atom_format_lines.size() != n_atoms) {
    logger << utils::LogLevel::SEVERE << "Can't print a system in PDB format baceuse the number of formatting lines differs from the number of atoms to be printed!";
    throw std::length_error("Incorrect number of lines for PDB format");
  }
#endif
    if (model_id != 0) where << utils::string_format("MODEL    %7d\n", model_id);
    for (size_t i = 0; i < n_atoms; i++) {
      C & ic = coordinates[i];
      where << utils::string_format(atom_format_lines[i], ic.wrap_x(), ic.wrap_y(), ic.wrap_z());
    }
    if (model_id != 0) where << "ENDMDL\n";
  }


  /** @brief Creates a system based on atoms from a given chain.
   *
   * The resulting system will have as many atoms as the source chain. The information about residues
   * the atoms belong to is lost.
   * @param c - all atoms of this Chain object will be copied to the newly created system.
   */
  CartesianAtomsSimple(core::data::structural::Chain & c) : n_atoms(c.count_atoms()), coordinates(new C[n_atoms]) {

    core::index4 i = 0;
    for (core::data::structural::Chain::atom_iterator ai = c.first_atom(); ai != c.last_atom(); ai++) {
      coordinates[i].set(*(*ai));
      i++;
    }
  }

  /** @brief Writes the state of this system in PDB format.
   * @param fname - name of the output file
   * @param model_id - model ID
   */
  void write_pdb(const std::string & fname, const core::index4 model_id = 0) const {
    std::ofstream out(fname);
    write_pdb(out, model_id);
    out.close();
  }


  /** @brief Distributes atoms of this system uniformly inside a simulation box.
   *
   * The system must be already created, so the number of its atoms is also defined.
   * The box size is read from an instance of C class (used to hold coordinates). This method makes
   * sense of course only in the case of simulations in periodic boundary conditions
   */
  static void atoms_in_box(CartesianAtomsSimple & system);

  /** @brief Distributes atoms of this system so they form a single chain.
   *
   * This method keeps the constant distance between subsequent beads and ensures that
   * any two atoms are not too close to each other. If a newly placed atom is closer to any other atom
   * that the given <code>cutoff</code>, then the new atom is placed again. In the case of  <code>n_bead_attempts</code>
   * such failures, the whole chain is discarded and build again. If that procedure fails to build a chain
   * <code>n_chain_attempts</code>, <code>false</code> is returned. Otherwise it returns true.
   *
   * @param system - the input system
   * @param bond_length - the distance between subsequent atoms, which is the lenth of a (virtual) bond
   * between them
   * @param cutoff - hard-core repulsion distance; any two atom (except the bonded neighbors) cannot
   * be closer to each other than this distance.
   * @param n_chain_attempts - the number of times the procedure is repeated (in the case of excluded volume violation)
   * @param n_bead_attempts - how many times a single bead is placed
   */
  static bool homopolymer(CartesianAtomsSimple<C> & system, const core::real bond_length, const core::real cutoff,
      const core::index2 n_chain_attempts = 1000, const core::index2 n_bead_attempts = 1000) {

    BuildPolymerChain<C> builder(system.coordinates,system.n_atoms);
    return builder.generate(3.8,4.0);
  }

private:
  std::string chain_name = "A";
  static utils::Logger logger;

};

template<class C>
utils::Logger CartesianAtomsSimple<C>::logger = utils::Logger("CartesianAtomsSimple");

/// Computes the center of mass of the whole system
template<class C>
static inline C cm(const CartesianAtomsSimple<C> & system) {

  C p;
  for (atom_index i = 0; i < system.count_atoms(); i++)
    p += (system.coordinates[i]);
  p /= ((core::real) system.count_atoms());

  return p;
}

template<class C>
void CartesianAtomsSimple<C>::atoms_in_box(CartesianAtomsSimple<C> & system) {

  core::real box_width = system.coordinates[0].get_box_len(); // get the size of the simulation box
  int n_each_side = ceil(pow(system.n_atoms, 1.0 / 3.0)); // count the number of atoms along each axis; some of these sites remain empty
  int ix = 0, iy = 0, iz = 0;
  core::real w = box_width / ((core::real) n_each_side);
  for (size_t i = 0; i < system.n_atoms; i++) {
    if (ix >= n_each_side) {
      ix = 0;
      iy++;
      if (iy >= n_each_side) {
        iy = 0;
        iz++;
        if (iz >= n_each_side) iz = 0;
      }
    }
    system.coordinates[i].x = w * (ix + 0.5);
    system.coordinates[i].y = w * (iy + 0.5);
    system.coordinates[i].z = w * (iz + 0.5);
    ix++;
  }
}

/** @brief Copies given coordinates into this chain.
 */
template<typename It, typename C>
void set_conformation(It begin, It end, CartesianAtomsSimple<C> & destination) {

  atom_index i = 0;
  for(It it=begin;it!=end;++it) {
    destination[i].set(**it);
    ++i;
  }
}

} // ~ simulations
} // ~ cartesian

#endif
