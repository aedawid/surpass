#ifndef SIMULATIONS_CARTESIAN_SurpassModel_HH
#define SIMULATIONS_CARTESIAN_SurpassModel_HH

#include <string>
#include <memory>

#include <core/real.hh>

#include <simulations/representations/surpass_utils.hh>
#include <simulations/systems/surpass/SurpassAtomTyping.hh>
#include <simulations/systems/ResidueChain.hh>

#include <utils/Logger.hh>

namespace simulations {
namespace systems {
namespace surpass {

/** @brief System in SURPASS representation
 *
 * @tparam C - the type used to express coordinates; Vec3 for simple models, Vec3Cubic for modeling in periodic boundary conditions
 */
template<class C>
class SurpassModel : public ResidueChain<C> {
public:

  /** @brief Creates a new SurpassModel object based on a biomolecular structure
   *
   * @param s - a Structure object - can be in SURPASS representation or all-atom
   */
  SurpassModel(core::data::structural::Structure &s);

  /// Necessary virtual destructor
  virtual ~SurpassModel() {}

  /** @brief Returns secondary structure element assignment for each surpass atom.
   *
   * The returned vector provides, for each surpass atom, the index of secondary structure element it belongs to,
   * for example for 2GB1 this vector has size 53 and for the native SS it is:
   * 11111000000022222000333333333333330000004444000005555
   * since 0 denotes a loop (which is not a surpass secondary structure element) and all 'true' SS elements are indexed from 1
   *  (size() = N_resids = n_atoms)
   */
  const std::vector<core::index2> &ss_element_for_atoms() const { return ss_element_for_atoms_; }

  /** @brief Returns beta strand index given atom from beta (size() = N_beta)
   * The returned index refers to <code>elements_beta_</code> vector, e.g. for 2GB1:
   *
   * 000001111122223333
   */
  const std::vector<core::index2> &beta_index_for_atoms() const { return beta_index_for_atoms_; }

  /** @brief Returns a list of all atoms of type S i.e. surpass beads that are beta-strands
   *
   * for example for 2GB1 this vector has size 18 since there are 18 residues in beta strands in that protein (size() = N_beta)
   * [0,1,2,3,4,12,13,14,15,16,40,41,42,43,49,50,51,52]
   */
  const std::vector<core::index2> &atoms_in_beta() const { return atoms_in_beta_; }

  /** @brief Returns a list of all atoms of type H i.e. surpass beads that are helical
   *
   * for example for 2GB1 this vector has size 14 since there are 14 residues in a helix in that protein  (size() = N_alpha)
   * [20,21,22,23,24,25,26,27,28,29,30,31,32,33]
   */
  const std::vector<core::index2> &atoms_in_alfa() const { return atoms_in_alfa_; }

  /** @brief List of these SS element indexes which are beta
   *
   * For 2GB1 this is [1,2,4,5] since 0 is for loops and 3 is helix
   *
   * @return list of beta SSE indexes
   */
  const std::vector<core::index2> &elements_beta() const { return elements_beta_; }

  /** @brief List of these SS element indexes which are alpha
   *
   * For 2GB1 this is [3] since 0 is for loops and 1,2,4 and 5 are strands
   *
   * @return list of alpha SSE indexes
   */
  const std::vector<core::index2> &elements_alfa() const { return elements_alfa_; }

  /** @brief Provide index of the first and last index for each helix (packed into a single vector)
   *
   * @return
   */
  const std::vector<core::index2> &alfa_ranges() const { return alfa_ranges_; }

private:
  utils::Logger logger;

  std::vector<core::index2> ss_element_for_atoms_;
  std::vector<core::index2> beta_index_for_atoms_;
  std::vector<core::index2> atoms_in_beta_;
  std::vector<core::index2> atoms_in_alfa_;
  std::vector<core::index2> elements_beta_;
  std::vector<core::index2> elements_alfa_;
  std::vector<core::index2> alfa_ranges_;

  void assign_ss_elements(void);
};

template<class C>
SurpassModel<C>::SurpassModel(core::data::structural::Structure &s) :
  ResidueChain<C>(std::make_shared<SurpassAtomTyping>(),
    (representations::is_surpass_model(s)) ? representations::fix_surpass_ss_assignment(s)
                                           : *representations::surpass_representation(s)), logger("SurpassModel") {

  assign_ss_elements();
}

template<typename C>
void SurpassModel<C>::assign_ss_elements(void) {

  core::index2 last_h = 0, index_E = 0, index_H = 0;
  ss_element_for_atoms_.resize(ResidueChain<C>::n_atoms);
  if (ResidueChain<C>::coordinates[0].atom_type != 2) {
    ++last_h;
    ss_element_for_atoms_[0] = last_h;
    if (ResidueChain<C>::coordinates[0].atom_type == 1) {
      atoms_in_beta_.push_back(0);
      elements_beta_.push_back(last_h);
      beta_index_for_atoms_.push_back(index_E);
    } else {
      atoms_in_alfa_.push_back(0);
      alfa_ranges_.push_back(0);//
      elements_alfa_.push_back(last_h);
      beta_index_for_atoms_.push_back(index_H);
    }
  } else {
    ss_element_for_atoms_[0] = last_h;
    beta_index_for_atoms_.push_back(ResidueChain<C>::n_atoms);
  }
  for (atom_index i = 1; i < ResidueChain<C>::n_atoms; ++i) {
    if (ResidueChain<C>::coordinates[i].atom_type == 2) {
      ss_element_for_atoms_[i] = 0;
      beta_index_for_atoms_.push_back(ResidueChain<C>::n_atoms);
    } else {
      if (ResidueChain<C>::coordinates[i - 1].atom_type != ResidueChain<C>::coordinates[i].atom_type)
        ++last_h;
      ss_element_for_atoms_[i] = last_h;
      if (ResidueChain<C>::coordinates[i].atom_type == 1) {
        atoms_in_beta_.push_back(i);
        if (elements_beta_.size() == 0) {
          elements_beta_.push_back(last_h);
          beta_index_for_atoms_.push_back(index_E);
        } else if (elements_beta_.back() != last_h) {
          elements_beta_.push_back(last_h);
          ++index_E;
          beta_index_for_atoms_.push_back(index_E);
        } else beta_index_for_atoms_.push_back(index_E);
      } else {
        atoms_in_alfa_.push_back(i);
        if ((ResidueChain<C>::coordinates[i + 1].atom_type != ResidueChain<C>::coordinates[i].atom_type) ||
            (i == ResidueChain<C>::n_atoms - 1))
          alfa_ranges_.push_back(i);//
        if (elements_alfa_.size() == 0) {
          alfa_ranges_.push_back(i);//
          elements_alfa_.push_back(last_h);
          beta_index_for_atoms_.push_back(index_H);
        } else if (elements_alfa_.back() != last_h) {
          elements_alfa_.push_back(last_h);
          alfa_ranges_.push_back(i);//
          ++index_H;
          beta_index_for_atoms_.push_back(index_H);
        } else beta_index_for_atoms_.push_back(index_H);
      }
      continue;
    }
  }
}

}
} // ~ simulations
} // ~ cartesian

#endif
