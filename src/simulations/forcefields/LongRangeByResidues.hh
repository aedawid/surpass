#ifndef SIMULATIONS_FORCEFIELDS_LongRangeByResidues_HH
#define SIMULATIONS_FORCEFIELDS_LongRangeByResidues_HH

#include <core/index.hh>
#include <core/data/basic/Array2D.hh>

#include <simulations/atom_indexing.hh>
#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/systems/ResidueChain.hh>

namespace simulations {
namespace forcefields {

/** @brief Base class for a long-range energy, e.g. Van der Waals, electrostatic etc.
 *
 * This base class provides methods that loop over a range of residues and calculate their energy by calling an energy
 * kernel method. The kernel must be provided by a derived class.
 *
 * @see calculate_by_residue(const core::index2) for details how energy is calculated.
 */
template<class C>
class LongRangeByResidues : public ByResidueEnergy {
public:

  /** @brief Constructs an instance that will evaluate energy of a given system.
   *
   * @param system - energy of this system will be calculated by methods of this class
   */
  LongRangeByResidues(const systems::ResidueChain<C> & system) :
      the_system(system), n_residues(the_system.count_residues()), n_atoms(the_system.n_atoms) {}

  virtual const std::string & name() const { return name_; }

  /** @brief Calculates energy of all atoms from <code>which_residue</code> interacting with all other residues in the system.
   *
   * Interactions will be calculated between residues <code>which_residue</code> and <code>another_residue</code> if and only if
   * <code>abs(which_residue-another_residue) >= offset</code>.
   * @param which_residue - the residue of interest
   * @return energy for a given residue
   */
  virtual inline double calculate_by_residue(const core::index2 which_residue) {

    double energy = 0.0;
    // --- Energy with residues preceding the residue of interest; if offset is equal to 0 then self energy will be evaluated here as well
    for (core::index2 ri = 0; ri <= which_residue - offset_; ++ri)
      if (!energy_kernel(which_residue, ri, energy)) return std::numeric_limits<double>::max();
    // --- Energy with residues following the residue of interest; if offset is equal to 0 then additional self-energy evaluation must avoided
    for (core::index2 ri = which_residue + offset_ + correct_for_zero_offset; ri < the_system.count_residues(); ++ri)
      if (!energy_kernel(which_residue, ri, energy)) return std::numeric_limits<double>::max();

    return energy;
  }

  /** @brief Calculates energy of all atoms from residues in the range from <code>chunk_from</code>
   * to  <code>chunk_to</code> (both inclusive) interacting with all other residues in the system. The residues
   * within the specified range interact with themselves.
   *
   * Interactions will be calculated between residues <code>which_residue</code> and <code>another_residue</code> if and only if
   * <code>abs(which_residue-another_residue) >= offset</code>.
   * @param chunk_from - the first residue of the range
   * @param chunk_to - the last residue of the range
   * @return energy for a given residue
   */
  virtual inline double calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to) {
//    std::cerr << "By chunk "<<chunk_from<<" "<<chunk_to<<"\n";

    double energy = 0.0;
    for (residue_index chunk_r = chunk_from; chunk_r <= chunk_to; ++chunk_r) {
      // --- Chunk interacting with upstream (e.g. N-terminal) residues
      for (residue_index ir = 0;
           ir <= std::min(int(chunk_r - correct_for_zero_offset - offset_), int(chunk_from-1)); ++ir) {
        if (!energy_kernel(chunk_r, ir, energy)) return std::numeric_limits<double>::max();
      }

      // Chunk interacting with downstream (e.g. C-terminal)  residues
      for (core::index2 ir = std::max((core::index2)(chunk_to+1), (core::index2) (chunk_r + correct_for_zero_offset + offset_)); ir < n_residues; ++ir) {
        if (!energy_kernel(chunk_r, ir, energy))return std::numeric_limits<double>::max();
      }
    }
    // Chunk interacting with itself
    for (core::index2 ir = chunk_from + offset_ ; ir <= chunk_to; ++ir) {
      for (core::index2 jr = chunk_from; jr <= ir-offset_; ++jr) {
        if (!energy_kernel(jr, ir, energy)) return std::numeric_limits<double>::max();
      }
    }

    return energy;
  }

  /** @brief Calculates energy of all atoms from <code>which_residue</code> interacting with all other residues in the system.
   *
   * Interactions will be calculated between residues <code>which_residue</code> and <code>another_residue</code> if and only if
   * <code>abs(which_residue-another_residue) >= offset</code>.
   *
   * This method also acumulates by-atom energy values for every atom pair considered during this call in the given map
   *
   * @param which_residue - the residue of interest
   * @param energy_map - 2D array used to accumulate by-atom energies
   * @return energy for a given residue
   */
  virtual inline double calculate_by_residue(const core::index2 which_residue, core::data::basic::Array2D<float> & energy_map) {

    double energy = 0.0;
    // --- Energy with residues preceding the residue of interest; if offset is equal to 0 then self energy will be evaluated here as well
    for (core::index2 ri = 0; ri <= which_residue - offset_; ++ri)
      if(!energy_kernel(which_residue, ri, energy, energy_map)) return std::numeric_limits<double>::max();
    // --- Energy with residues following the residue of interest; if offset is equal to 0 then additional self-energy evaluation must avoided
    for (core::index2 ri = which_residue + offset_ + correct_for_zero_offset; ri < the_system.count_residues(); ++ri) {
      if(!energy_kernel(which_residue, ri, energy, energy_map)) return std::numeric_limits<double>::max();
    }

    return energy;
  }

  /** @brief Calculates energy of all atoms from residues in the range from <code>chunk_from</code>
   * to  <code>chunk_to</code> (both inclusive) interacting with all other residues in the system. The residues
   * within the specified range interact with themselves.
   *
   * Interactions will be calculated between residues <code>which_residue</code> and <code>another_residue</code> if and only if
   * <code>abs(which_residue-another_residue) >= offset</code>.
   *
   * This method also acumulates by-atom energy values for every atom pair considered during this call in the given map
   *
   * @param chunk_from - the first residue of the range
   * @param chunk_to - the last residue of the range
   * @param energy_map - 2D array used to accumulate by-atom energies
   * @return energy for a given residue
   */
  virtual inline double calculate_by_chunk(const core::index2 chunk_from,
    const core::index2 chunk_to, core::data::basic::Array2D<float> & energy_map) {

    double energy = 0.0;
    for (residue_index chunk_r = chunk_from; chunk_r <= chunk_to; ++chunk_r) {
      // --- Chunk interacting with upstream (e.g. N-terminal) residues
      for (residue_index ir = 0;
           ir <= std::min(int(chunk_from - correct_for_zero_offset - offset_), int(chunk_from)); ++ir) {
        if(!energy_kernel(chunk_r, ir, energy, energy_map)) return std::numeric_limits<double>::max();
      }

      // Chunk interacting with downstream (e.g. C-terminal)  residues
      for (core::index2 ir = std::max((core::index2)(chunk_to+1), (core::index2) (chunk_r + correct_for_zero_offset + offset_)); ir < n_residues; ++ir) {
        if(!energy_kernel(chunk_r,ir,energy, energy_map)) return std::numeric_limits<double>::max();
      }
    }
    // Chunk interacting with itself
    for (core::index2 ir = chunk_from + offset_ ; ir <= chunk_to; ++ir) {
      for (core::index2 jr = chunk_from; jr <= ir-offset_; ++jr) {
        if(!energy_kernel(jr, ir, energy, energy_map)) return std::numeric_limits<double>::max();
      }
    }

    return energy;
  }

  /** @brief Calls energy kernel for a given pair of residues.
   *
   * @param the_moved_residue - the residue moved by a mover
   * @param the_other_residue - another residue to calculate the pairwise energy
   * @param energy - energy value will be accumulated there
   * @return true if the energy is finite, i.e. there is no 'hard' reason to reject the curernt MC move
   */
  virtual bool energy_kernel(const core::index2 the_moved_residue, const core::index2 the_other_residue, double & energy) = 0;

  /** @brief Calls energy kernel and stores by-atom-pair energy values in a matrix.
   *
   * Default implementation just calls <code>bool energy_kernel(const core::index2 the_moved_residue, const core::index2 the_other_residue, double & energy)</code>
   * and does not write anything to the <code>energy_map</code>. Override this method to change the default behavior
   * @param the_moved_residue - the residue moved by a mover
   * @param the_other_residue - another residue to calculate the pairwise energy
   * @param energy - energy value will be accumulated there
   * @param energy_map -  by-atom-pair energy values will be accumulated in this matrix.
   * @return true if the energy is finite, i.e. there is no 'hard' reason to reject the curernt MC move
   */
  virtual bool energy_kernel(const core::index2 the_moved_residue, const core::index2 the_other_residue,
                             double & energy, core::data::basic::Array2D<float> & energy_map) {
    return energy_kernel(the_moved_residue, the_other_residue, energy);
  };

  virtual inline double calculate() {

    double energy = 0.0;
    for (residue_index k = offset_; k < n_residues; ++k) {
      for (residue_index i = 0; i <= k - offset_; i++) {
        if(!energy_kernel(k, i, energy)) return std::numeric_limits<double>::max();
      }
    }
    return energy;
  }

  virtual inline double calculate(core::data::basic::Array2D<float> & energy_map) {

    double energy = 0.0;
    for (residue_index k = offset_; k < n_residues; ++k) {
      for (residue_index i = 0; i <= k - offset_; i++) {
        if(!energy_kernel(k, i, energy, energy_map)) return std::numeric_limits<double>::max();
      }
    }
    return energy;
  }

  /** @brief sets the sequence separation for long-range energy calculations.
   *
   * Energy between residues <code>i</code> and <code>j</code> is calculated if and only if \f$|i-j| \ge o\f$
   * where <code>o</code> is the offset.
   * @param offset - the new offset value
   */
  void residue_offset(const unsigned char offset) {
    offset_ = offset;
    correct_for_zero_offset = (offset == 0) ? 1 : 0;
  }

  /** @brief Returns the currently used value of residue offset
   * @see residue_offset(const unsigned char)
   * @return residue offset
   */
  unsigned char residue_offset() const { return offset_; }

protected:
  const systems::ResidueChain<C> & the_system; ///< the system whose energy will be evaluated
  const core::index2 n_residues; ///< the number of residues in the system
  const core::index4 n_atoms; ///< the number of atoms in the system
  unsigned char offset_ = 1;  ///< energy between residues <code>i</code> and <code>j</code> is calculated if and only if \f$|i-j| \ge o\f$ where <code>o</code> is the offset_

private:
  static const std::string name_;
  core::index2 correct_for_zero_offset = 0; ///< This is necessary when offset is 0 so self-energy is not computed twice
};

template<typename C>
const std::string LongRangeByResidues<C>::name_ = "LongRangeByResidues";

} // ~ cartesian
} // ~ ff

#endif
