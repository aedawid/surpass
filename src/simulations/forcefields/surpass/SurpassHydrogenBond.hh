#ifndef SIMULATIONS_CARTESIAN_FF_SurpassHydrogenBond_HH
#define SIMULATIONS_CARTESIAN_FF_SurpassHydrogenBond_HH

#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <simulations/atom_indexing.hh>
#include <core/data/basic/Array2D.hh>
#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/systems/ResidueChain.hh>
#include <simulations/systems/surpass/SurpassModel.hh>
#include <core/algorithms/UnionFind.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

using core::real;
using simulations::atom_index;
using simulations::residue_index;

template<class C>
class SurpassHydrogenBond : public ByResidueEnergy {

public:

  SurpassHydrogenBond(const systems::surpass::SurpassModel <C> &system) :
    the_system(system), n_residues(the_system.count_residues()), n_atoms(the_system.n_atoms),
    hydrogen_bonds_(the_system.atoms_in_beta().size(), std::pair<core::index4, core::index4>()),
    beta_topology_matrix_((the_system.elements_beta().size() == 0) ? 1 : the_system.elements_beta().size(),
      (the_system.elements_beta().size() == 0) ? 1 : the_system.elements_beta().size()),
    count_matrix_((the_system.elements_beta().size() == 0) ? 1 : the_system.elements_beta().size(),
      (the_system.elements_beta().size() == 0) ? 1 : the_system.elements_beta().size()) {
    for (unsigned int i = 0; i < the_system.elements_beta().size(); ++i) union_find_sheets_.add_element(i);
    find_hydrogen_bonds();
//    topology_pairs_.resize((the_system.elements_beta().size()) * (the_system.elements_beta().size()));
//    beta_pairs();
  }

  SurpassHydrogenBond(const systems::surpass::SurpassModel <C> &system, const std::vector<std::string> &parameters)
    : SurpassHydrogenBond(system) {}

  // (1) distance no longer than 6A
  // (2) cos(angle between strands) > 0.57 to avoid warkoczyki (braids)
  // from all of these which satisfy (1) and (2) find best angle for a single HB ([80,110])
  virtual void find_hydrogen_bonds() {

    beta_topology_matrix_.clear(0);
    count_matrix_.clear(0);
    union_find_sheets_.disconnect();
    if (the_system.atoms_in_beta().size() == 0) return; //skip if protein is only alfa type

    core::index2 ss = the_system.ss_element_for_atoms()[the_system.atoms_in_beta()[0]];    // index of the first beta
    core::index4 index_y = 0, ID_j, id2, id3;
    core::real dist = 6.0, H1H2 = 0.0, H1H3 = 0.0, angle1 = 0.0, diff = 0.35, r = 0.0;    //value -> angle, good if is in range 70-110; diff -> 90-20=70 & 90+20=110; r -> length of HB;
    core::real value = 0.0, angle2 = 125.0, HB_cos = 1.0;                    //the lowest acceptable value of angle between 2 hydrogen bonds (vectors); the closer to 180 the better
    bool is_good = false;
    std::vector<core::index4> hbond_partners;                            //vector with possible ACCEPTORS; the best one (by distance) of each strand (if applicable)
    std::vector<unsigned int> index, ID_remove, index4;
    core::data::basic::Vec3 H1, H2, H3;                            //vectors: y-j, y-y+2, y-j(n); to check the angles 1) between two hydrogen bonds 2) if B-sheet is linear
    unsigned int index1 = 0, index2 = 0, index3 = 0;                        //indexes of connected beta_ss_elements in the beta_topology_matrix_ matrix
//--1-- For each beta residue (DONOR HB) find all possible hydrogen bond partners (ACCEPTORS)
    // y - donor HB, j - acceptor HB (true indexes of each residue in chain)
    for (auto &y : the_system.atoms_in_beta()) {                        //the first atom in hydrogen bond (DONOR_id, core::index2)
      hbond_partners.clear();
      ID_j = n_atoms;
      ID_remove.clear();///
      index4.clear();///
//--2-- From all possible ACCEPTORS select one the best for each B-strand (if applicable)
      for (auto &j : the_system.atoms_in_beta()) {                        //the second atom in hydrogen bond (ACCEPTOR_id, core::index2)
        if (the_system.ss_element_for_atoms()[y] != the_system.ss_element_for_atoms()[j]) {    //if both residue (y and j) are laying on different B-strands...
          if (ss != the_system.ss_element_for_atoms()[j]) {                    //if we start iterating on new B-strand or j is the last residue in beta...
            dist = 6.0;                                    //the longest possible hydrogen bond (6.0A)
            diff = 20.0;                                    //the higher deviation for angle
            ss = the_system.ss_element_for_atoms()[j];
          }
//--3-- find the best ACCEPTOR in actual B-strand:
          r = the_system[y].distance_to(the_system[j]);
          if (r <= dist) {                                    //if hydrogen bond length is good (if less than 6.0 or shorter than previous)
            H1 = (the_system[j]) - (the_system[y]);
            if (y == the_system.atoms_for_chain(0).first_atom) H2 = (the_system[y + 2]) - (the_system[y]);
            else if (y == the_system.atoms_for_chain(0).last_atom) H2 = (the_system[y]) - (the_system[y - 2]);
            else H2 = (the_system[y + 1]) - (the_system[y - 1]);                    //vectors for angle are always defined only between beta residues
            if (j == the_system.atoms_for_chain(0).first_atom) H3 = (the_system[j + 2]) - (the_system[j]);
            else if (j == the_system.atoms_for_chain(0).last_atom) H3 = (the_system[j]) - (the_system[j - 2]);
            else H3 = (the_system[j + 1]) - (the_system[j - 1]);
            HB_cos = fabs(((H2.x * H3.x) + (H2.y * H3.y) + (H2.z * H3.z)) / (H2.length() * H3.length()));
            if (HB_cos > 0.57) {
              H1H2 = fabs(((H1.x * H2.x) + (H1.y * H2.y) + (H1.z * H2.z)) / (H1.length() * H2.length()));
              H1H3 = fabs(((H1.x * H3.x) + (H1.y * H3.y) + (H1.z * H3.z)) / (H1.length() * H3.length()));
              if (H1H2 < H1H3) angle1 = H1H2;
              else angle1 = H1H3;
              if (diff >= fabs(angle1)) {                        //if angle value is better
                diff = fabs(angle1);                            //the best angle value
                dist = r;                                    //the shortest hydrogen bond for actual strand
                ID_j = j;                                    //id of actual the best ACCEPTOR of this strand
                if (j == the_system.atoms_in_beta().back()) {
                  if (ID_j != n_atoms) {
                    hbond_partners.push_back(ID_j);
                    ID_j = n_atoms;
                  }
                } else if (j != the_system.atoms_for_chain(0).last_atom) {
                  if ((the_system.ss_element_for_atoms()[j + 1] != the_system.ss_element_for_atoms()[j]) &&
                      (ID_j != n_atoms)) {
                    hbond_partners.push_back(ID_j);
                    ID_j = n_atoms;
                  }
                }
              }
            }
          } else if (j == the_system.atoms_in_beta().back()) {
            if (ID_j != n_atoms) {
              hbond_partners.push_back(ID_j);
              ID_j = n_atoms;
            }
          } else if (the_system.ss_element_for_atoms()[j + 1] != the_system.ss_element_for_atoms()[j]) {
            if (ID_j != n_atoms) {
              hbond_partners.push_back(ID_j);
              ID_j = n_atoms;
            }
          }
        }
      }
//--4-- From all selected ACCEPTORS choose only two best
      is_good = false;
      angle2 = 125.0;
      if (hbond_partners.size() >= 2) {
        for (unsigned int k = 0; k < hbond_partners.size() - 1; k++) {
          for (unsigned int l = 1; l < hbond_partners.size() - k; l++) {
            H1 = (the_system[hbond_partners[k]]) - (the_system[y]);                //vector of first h-bond
            H2 = (the_system[hbond_partners[k + l]]) - (the_system[y]);            //vector of next h-bond
            H1H2 = (H1.x * H2.x) + (H1.y * H2.y) + (H1.z * H2.z);
            value = acos(H1H2 / (H1.length() * H2.length())) * 180.0 /
                    M_PI;            //value of angle between 2 hydrogen bonds
            if (value >=
                angle2) {                                //if value for planar hydrogen bonds is good (close to 180)
              is_good = true;
              angle2 = value;
              id2 = hbond_partners[k];                            //id (acceptor first h-bond)
              id3 = hbond_partners[k + l];                            //id (acceptor second h-bond)
            } else if (count_matrix_.get(the_system.beta_index_for_atoms()[y],
              the_system.beta_index_for_atoms()[hbond_partners[k]]) >
                       count_matrix_.get(the_system.beta_index_for_atoms()[y],
                         the_system.beta_index_for_atoms()[hbond_partners[k + l]])) {
//              } else if (the_system[y].distance_to(the_system[hbond_partners[k]]) <= the_system[y].distance_to(the_system[hbond_partners[k+l]])) {
              ID_j = hbond_partners[k];
              if (hbond_partners[k + l] < y) ID_remove.push_back(hbond_partners[k + l]);///
            } else {
              ID_j = hbond_partners[k + l];
              if (hbond_partners[k] < y) ID_remove.push_back(hbond_partners[k]);///
            }
          }
        }
        if (is_good == true) {
          hydrogen_bonds_[index_y] = std::make_pair(id2, id3);
          if (ID_remove.size() != 0) {///
            for (unsigned int pp = 0; pp < ID_remove.size(); ++pp) {///
              if ((ID_remove[pp] == id2) || (ID_remove[pp] == id3)) ID_remove.erase(ID_remove.begin() + pp);///
            }///
          }///
          index1 = (unsigned char) the_system.beta_index_for_atoms()[y];
          index2 = (unsigned char) the_system.beta_index_for_atoms()[id2];
          index3 = (unsigned char) the_system.beta_index_for_atoms()[id3];
          for (unsigned int pp = 0; pp < ID_remove.size(); ++pp) {///
            index4.push_back((unsigned char) the_system.beta_index_for_atoms()[ID_remove[pp]]);///
          }///
          count_matrix_.set(index1, index2, count_matrix_.get(index1, index2) + 1);
          count_matrix_.set(index1, index3, count_matrix_.get(index1, index3) + 1);
        } else {
          hydrogen_bonds_[index_y] = std::make_pair(ID_j, n_atoms);
          index1 = (unsigned char) the_system.beta_index_for_atoms()[y];
          index2 = (unsigned char) the_system.beta_index_for_atoms()[ID_j];
          for (unsigned int pp = 0; pp < ID_remove.size(); ++pp) {///
            index4.push_back((unsigned char) the_system.beta_index_for_atoms()[ID_remove[pp]]);///
          }///
          count_matrix_.set(index1, index2, count_matrix_.get(index1, index2) + 1);
        }
        if (ID_remove.size() != 0) {///
          for (unsigned int p = 0; p != the_system.atoms_in_beta().size(); ++p) {///
            for (unsigned int pp = 0; pp < ID_remove.size(); ++pp) {
              if (the_system.atoms_in_beta()[p] == ID_remove[pp]) index3 = p;
              {///
                if (hydrogen_bonds_[index3].first == y) hydrogen_bonds_[index3].first = n_atoms;///
                else if (hydrogen_bonds_[index3].second == y) hydrogen_bonds_[index3].second = n_atoms;///
              }
            }///
          }///
          for (unsigned int pp = 0; pp < index4.size(); ++pp) {///
            if (count_matrix_.get(index4[pp], index1) > 0)
              count_matrix_.set(index4[pp], index1, count_matrix_.get(index4[pp], index1) - 1);///
          }///
        }///
      } else if (hbond_partners.size() == 1) {
        hydrogen_bonds_[index_y] = std::make_pair(hbond_partners[0], n_atoms);
        index1 = (unsigned char) the_system.beta_index_for_atoms()[y];
        index2 = (unsigned char) the_system.beta_index_for_atoms()[hbond_partners[0]];
        count_matrix_.set(index1, index2, count_matrix_.get(index1, index2) + 1);
      } else hydrogen_bonds_[index_y] = std::make_pair(n_atoms, n_atoms);
      ++index_y;                                        //DONOR's index in vector: atoms_in_beta
    }
//--5-- Assign "1" for real HB and "2" for all strands (j,k, etc.) connected by hydrogen bond with the same strand (i)
    for (unsigned int ai = 0; ai < the_system.elements_beta().size(); ++ai) {
      index.clear();
      for (unsigned int aj = ai; aj < the_system.elements_beta().size(); ++aj) {
        if ((count_matrix_.get(ai, aj) > 0) && (count_matrix_.get(aj, ai) > 0)) {
          beta_topology_matrix_.set(ai, aj, 1);
          beta_topology_matrix_.set(aj, ai, 1);
          index.push_back(aj);
          union_find_sheets_.union_set(ai, aj);
        }
      }
      if (index.size() >= 2) {
        for (unsigned int ak = 0; ak < index.size() - 1; ++ak) {
          for (unsigned int al = ak + 1; al < index.size(); ++al) {
            if (beta_topology_matrix_.get(index[ak], index[al]) != 1)
              beta_topology_matrix_.set(index[ak], index[al], 2);
            if (beta_topology_matrix_.get(index[al], index[ak]) != 1)
              beta_topology_matrix_.set(index[al], index[ak], 2);
          }
        }
      }
    }
  }

  virtual const std::string &name() const { return name_; }

  virtual const std::vector<std::pair<core::index4, core::index4>> &hydrogen_bonds() const { return hydrogen_bonds_; }

  virtual const core::data::basic::Array2D<core::index1> &beta_topology_matrix() const { return beta_topology_matrix_; }

  virtual const core::data::basic::Array2D<core::index1> &count_matrix() const { return count_matrix_; }


  /** @brief Returns UnionFind object used to gather beta strands into sheets.
   *
   * @return
   */
  virtual const core::algorithms::UnionFind<unsigned int, unsigned int> &
  union_find_sheets() const { return union_find_sheets_; }

  virtual inline double calculate_by_residue(const residue_index which_residue) {

    return calculate_by_residue_rehash(which_residue, the_system[which_residue].atom_type == 1);
  }

  inline double calculate_by_residue_rehash(const residue_index which_residue, const bool if_rehash_hbonds) {

    if (if_rehash_hbonds) find_hydrogen_bonds();

    double en = 0.0;
    core::real r = 0.0;
    atom_index y = the_system.atoms_for_residue(which_residue).first_atom;
    const auto it = lower_bound(the_system.atoms_in_beta().begin(), the_system.atoms_in_beta().end(), y);
    if ((it != the_system.atoms_in_beta().end()) && (*it == y)) {
      core::index4 i = it - the_system.atoms_in_beta().begin();
      if (hydrogen_bonds_[i].first != n_atoms) {
        r = the_system[y].distance_to(the_system[hydrogen_bonds_[i].first]);
        en += -log((exp(-(r - 4.65) * (r - 4.65)) + 0.57) / 0.57);    //premium for 1 or 2 the best h-bond
      }
      if (hydrogen_bonds_[i].second != n_atoms) {
        r = the_system[y].distance_to(the_system[hydrogen_bonds_[i].second]);
        en += -log((exp(-(r - 4.65) * (r - 4.65)) + 0.57) / 0.57);    //premium for 1 or 2 the best h-bond
      }
    }

    return en;
  }

  virtual inline double calculate_by_chunk(const simulations::residue_index chunk_from,
                                           const simulations::residue_index chunk_to) {

    double en = 0.0;
    find_hydrogen_bonds();
    atom_index y_from = the_system.atoms_for_residue(chunk_from).first_atom;
    atom_index y_to = the_system.atoms_for_residue(chunk_to).first_atom;
    for (atom_index y = y_from; y <= y_to; ++y) {
      en += calculate_by_residue_rehash(y, false);
    }

    return en;
  }

  virtual inline double calculate() {

    double en = 0.0;
    find_hydrogen_bonds();
    for (auto &y : the_system.atoms_in_beta()) en += calculate_by_residue_rehash(y, false);
    return en;
  }

  static core::index2 n_ss_elements(const systems::surpass::SurpassModel <C> &system) {
    return *(std::max_element(system.ss_element_for_atoms().begin(), system.ss_element_for_atoms().end()));
  }

protected:
  const systems::surpass::SurpassModel <C> &the_system; ///< the system whose energy will be evaluated
  const core::index2 n_residues; ///< the number of residues in the system
  const core::index4 n_atoms; ///< the number of atoms in the system

private:
  static const std::string name_;
  std::vector<std::pair<core::index4, core::index4>> hydrogen_bonds_;
  core::data::basic::Array2D<core::index1> beta_topology_matrix_; ///< topology matrix for beta only (of size n_beta x n_beta), contain 0 or 1 (when two strands are H-bonded)
  core::data::basic::Array2D<core::index1> count_matrix_; ///< provides the count of hydrogen bonds between two strands (of size n_beta x n_beta)
  core::algorithms::UnionFind<unsigned int, unsigned int> union_find_sheets_; ///< binds beta strands into sheets


  inline void vec_along(core::index4 i_atom, Vec3 &output) {

    if (i_atom == the_system.atoms_for_chain(0).first_atom) {
      output.set(the_system[i_atom + 2]);
      output -= the_system[i_atom];
      return;
    }
    if (i_atom == the_system.atoms_for_chain(0).last_atom) {
      output.set(the_system[i_atom + 2]);
      output -= the_system[i_atom];
      return;
    }
    output.set(the_system[i_atom + 1]);
    output -= the_system[i_atom - 2];
  }

  void find_H_acceptors(core::index2 y, std::vector<core::index4> &hbond_partners) {

    core::real dist = 6.0, H1H2 = 0.0, H1H3 = 0.0, angle1 = 0.0, diff = 0.35, r = 0.0;
    core::index2 ss = the_system.ss_element_for_atoms()[the_system.atoms_in_beta()[0]];
    core::index4 ID_j = n_atoms;

    Vec3 H1, H2, H3;
//--2-- From all possible ACCEPTORS select one the best for each B-strand (if applicable)
    for (auto &j : the_system.atoms_in_beta()) {                        //the second atom in hydrogen bond (ACCEPTOR_id, core::index2)
      if (the_system.ss_element_for_atoms()[y] !=
          the_system.ss_element_for_atoms()[j]) {    //if both residue (y and j) are laying on different B-strands...
        if (ss !=
            the_system.ss_element_for_atoms()[j]) {                    //if we start iterating on new B-strand or j is the last residue in beta...
          dist = 6.0;                                    //the longest possible hydrogen bond (6.0A)
          diff = 20.0;                                    //the higher deviation for angle
          ss = the_system.ss_element_for_atoms()[j];
        }
//--3-- find the best ACCEPTOR in actual B-strand:
        r = the_system[y].distance_to(the_system[j]);
        if (r <= dist) {                                    //if hydrogen bond length is good (if less than 6.0 or shorter than previous)
          H1.set(the_system[j]);
          H1-=the_system[y];
          vec_along(y,H2);
          vec_along(j,H3);

          core::real HB_cos = fabs(((H2.x * H3.x) + (H2.y * H3.y) + (H2.z * H3.z)) / (H2.length() * H3.length()));
          if (HB_cos > 0.57) {
            H1H2 = fabs(((H1.x * H2.x) + (H1.y * H2.y) + (H1.z * H2.z)) / (H1.length() * H2.length()));
            H1H3 = fabs(((H1.x * H3.x) + (H1.y * H3.y) + (H1.z * H3.z)) / (H1.length() * H3.length()));
            if (H1H2 < H1H3) angle1 = H1H2;
            else angle1 = H1H3;
            if (diff >= fabs(0.0 - angle1)) {                        //if angle value is better
              diff = fabs(0.0 - angle1);                            //the best angle value
              dist = r;                                    //the shortest hydrogen bond for actual strand
              ID_j = j;                                    //id of actual the best ACCEPTOR of this strand
              if (j == the_system.atoms_in_beta().back()) {
                if (ID_j != n_atoms) {
                  hbond_partners.push_back(ID_j);
                  ID_j = n_atoms;
                }
              } else if (j != the_system.atoms_for_chain(0).last_atom) {
                if ((the_system.ss_element_for_atoms()[j + 1] != the_system.ss_element_for_atoms()[j]) &&
                    (ID_j != n_atoms)) {
                  hbond_partners.push_back(ID_j);
                  ID_j = n_atoms;
                }
              }
            }
          }
        } else if (j == the_system.atoms_in_beta().back()) {
          if (ID_j != n_atoms) {
            hbond_partners.push_back(ID_j);
            ID_j = n_atoms;
          }
        } else if (the_system.ss_element_for_atoms()[j + 1] != the_system.ss_element_for_atoms()[j]) {
          if (ID_j != n_atoms) {
            hbond_partners.push_back(ID_j);
            ID_j = n_atoms;
          }
        }
      }
    }
  }
};

template<typename C>
const std::string SurpassHydrogenBond<C>::name_ = "SurpassHydrogenBond";

} // ~ simulations
} // ~ cartesian
} // ~ ff

#endif
