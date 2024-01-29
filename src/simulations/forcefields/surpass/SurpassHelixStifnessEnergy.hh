#ifndef SIMULATIONS_GENERIC_FF_SurpassHelixStifnessEnergy_HH
#define SIMULATIONS_GENERIC_FF_SurpassHelixStifnessEnergy_HH

#include <vector>
#include <memory>
#include <iomanip>
#include <math.h>

#include <core/real.hh>
#include <core/index.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

#include <simulations/systems/surpass/SurpassModel.hh>
#include <simulations/forcefields/LongRangeByResidues.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

template<class C>
class SurpassHelixStifnessEnergy : public ByResidueEnergy {
public:

  SurpassHelixStifnessEnergy(const systems::surpass::SurpassModel<C> & system) : the_system(system) {}

  SurpassHelixStifnessEnergy(const systems::surpass::SurpassModel<C> & system, const std::vector<std::string> & parameters) : SurpassHelixStifnessEnergy(system) {}

  const std::string & name() const { return name_; }

  virtual double calculate_by_residue(const core::index2 which_residue) {

    double energy = 0.0;
    unsigned int H_len = 0;
    if (the_system[which_residue].atom_type == 0) {
      unsigned int index = (unsigned char)the_system.beta_index_for_atoms()[which_residue];
      unsigned int index1 = the_system.alfa_ranges()[2*index + 1];
      unsigned int index2 = the_system.alfa_ranges()[2*index];
      H_len = index1 - index2 + 1;
      core::real D = 1.45 * H_len - 1.8;
      core::real d = the_system[index1].distance_to(the_system[index2]);
      if ((H_len >= 10)&&(H_len <= 30)&&(d < D - 2)) energy = D - 2 - d;
      else if (d < D) energy = D - d;
    }
    return energy*H_len;
  }

  virtual double calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to) {

    double energy = 0.0;
    unsigned int index1 = (unsigned char)the_system.beta_index_for_atoms()[chunk_from];
    unsigned int index2 = (unsigned char)the_system.beta_index_for_atoms()[chunk_to];
    unsigned int index3 = 0;
    if (index1 != index2) {
      for (unsigned int i = index1; i <= index2; ++i) {
	index3 = the_system.alfa_ranges()[2*i];
	energy += calculate_by_residue(index3);
      }
    } else energy = calculate_by_residue(chunk_from);
    return energy;
  }

  virtual double calculate() {

    double energy = 0.0;
    for (unsigned int i = 0; i < the_system.elements_alfa().size(); ++i) {
      unsigned int index = the_system.alfa_ranges()[2*i];
      energy += calculate_by_residue(index);
    }
    return energy;
  }

protected:
  const systems::surpass::SurpassModel<C> & the_system; ///< the system whose energy will be evaluated

private:
  static utils::Logger logger;
  static const std::string name_;
};

template<typename C>
const std::string SurpassHelixStifnessEnergy<C>::name_ = "SurpassHelixStifnessEnergy";

template<typename C>
utils::Logger SurpassHelixStifnessEnergy<C>::logger("SurpassHelixStifnessEnergy");


} // ~ simulations
} // ~ generic
} // ~ ff
#endif
