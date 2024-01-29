#ifndef SIMULATIONS_GENERIC_FF_SurpassCentrosymetricEnergy_HH
#define SIMULATIONS_GENERIC_FF_SurpassCentrosymetricEnergy_HH

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
class SurpassCentrosymetricEnergy : public ByResidueEnergy {
public:

  SurpassCentrosymetricEnergy(const systems::surpass::SurpassModel<C> & system) : the_system(system) {}

  SurpassCentrosymetricEnergy(const systems::surpass::SurpassModel<C> & system, const std::vector<std::string> & parameters) : SurpassCentrosymetricEnergy(system) {}

  const std::string & name() const { return name_; }

  virtual double calculate_by_residue(const core::index2 which_residue) {

    core::index2 k = the_system.chain_for_atom(which_residue); // --- index of the chain having that residue
    const systems::AtomRange<C> &ch_range = the_system.atoms_for_chain(k);
    core::index4 n_atoms = ch_range.size();
    core::real R = sqrt(n_atoms + 18);
    core::data::basic::Vec3 center;
    unsigned int n = 0;
    for (atom_index i = ch_range.first_atom; i <= ch_range.last_atom; ++i) center += the_system[i];
    center /= n_atoms;

    core::real R2 = R * R;
    for (atom_index j = ch_range.first_atom; j <= ch_range.last_atom; ++j) {
      core::real x = the_system[j].x - center.x;
      core::real r2 = x * x;
      if (r2 > R2) continue;
      x = the_system[j].y - center.y;
      r2 += x * x;
      if (r2 > R2) continue;
      x = the_system[j].z - center.z;
      r2 += x * x;
      if (r2 < R2) ++n;
    }
    double energy = 0.0;
    if (n < 0.5 * n_atoms) energy = 0.2 * sqrt((int(0.5 * n_atoms)) - n) + 5;
    else if (n > 0.6 * n_atoms) energy = 5 * (n - int(0.6 * n_atoms));
    return energy;
  }

  virtual double calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to) {

    double energy = calculate_by_residue(chunk_from);
    return energy;
  }

  virtual double calculate() {
    core::index2 chains = the_system.count_chains();
    double energy = 0.0;
    for (core::index2 k = 0; k < chains; ++k) {
      core::index4 which_residue = the_system.atoms_for_chain(k).last_atom;
      double en = calculate_by_residue(which_residue);
      energy += en;
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
const std::string SurpassCentrosymetricEnergy<C>::name_ = "SurpassCentrosymetricEnergy";

template<typename C>
utils::Logger SurpassCentrosymetricEnergy<C>::logger("SurpassCentrosymetricEnergy");


} // ~ simulations
} // ~ generic
} // ~ ff
#endif
