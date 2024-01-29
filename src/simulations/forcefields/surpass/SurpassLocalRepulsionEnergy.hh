#ifndef SIMULATIONS_GENERIC_FF_SurpassLocalRepulsionEnergy_HH
#define SIMULATIONS_GENERIC_FF_SurpassLocalRepulsionEnergy_HH

#include <vector>
#include <memory>
#include <iomanip>
#include <math.h>
#include <cmath>

#include <core/index.hh>
#include <core/real.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

#include <simulations/systems/surpass/SurpassModel.hh>
#include <simulations/forcefields/LongRangeByResidues.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

template<class C>
class SurpassLocalRepulsionEnergy : public ByResidueEnergy {
public:

  /// radius of 'small ball' - the sphere around which_residue
  const static core::real SURPASS_LOCAL_REPULSION_CUTOFF;

  SurpassLocalRepulsionEnergy(const systems::surpass::SurpassModel <C> &system) : the_system(system) {}

  SurpassLocalRepulsionEnergy(const systems::surpass::SurpassModel <C> &system,
                              const std::vector<std::string> &parameters) : SurpassLocalRepulsionEnergy(system) {}

  const std::string &name() const { return name_; }

  virtual double calculate_by_residue(const simulations::residue_index which_residue) {

    const static core::real R2 = SURPASS_LOCAL_REPULSION_CUTOFF * SURPASS_LOCAL_REPULSION_CUTOFF;
    core::index2 N; // --- the maximum number of neighbours in the sphere
    if (the_system[which_residue].atom_type == 0) N = 2;    //for ss=H, N = 2
    else if (the_system[which_residue].atom_type == 2) N = 4;    //for ss=C, N = 4
    else N = 6;                         //for ss=E, N = 6
    unsigned int n = 0;   // --- actual number of neighbours in the sphere
    core::index4 max_j = (which_residue<3) ? 0 : which_residue -3;
    for (core::index4 j = 0; j < max_j; ++j) {

      core::real x = the_system[j].x - the_system[which_residue].x;
      core::real r2 = x * x;
      if (r2 > R2) continue;
      x = the_system[j].y - the_system[which_residue].y;
      r2 += x * x;
      if (r2 > R2) continue;
      x = the_system[j].z - the_system[which_residue].z;
      r2 += x * x;
      if (r2 < R2) ++n;
    }

    for (core::index4 j = which_residue + 4; j < the_system.n_atoms; ++j) {

      core::real x = the_system[j].x - the_system[which_residue].x;
      core::real r2 = x * x;
      if (r2 > R2) continue;
      x = the_system[j].y - the_system[which_residue].y;
      r2 += x * x;
      if (r2 > R2) continue;
      x = the_system[j].z - the_system[which_residue].z;
      r2 += x * x;
      if (r2 < R2) ++n;
    }
    double energy = 0.0;
    if (n > N) energy = n - N;
    return energy * energy;
  }

  virtual double calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to) {
    double en = 0.0;
    for (core::index4 y = chunk_from; y <= chunk_to; ++y) en += calculate_by_residue(y);
    return en;
  }

  virtual double calculate() {

    core::index2 chains = the_system.count_chains();
    double en = 0.0;
    for (core::index2 k = 0; k < chains; ++k) {
      for (core::index4 y = the_system.atoms_for_chain(k).first_atom; y <= the_system.atoms_for_chain(k).last_atom; ++y) en += calculate_by_residue(y);
    }
    return en;
  }

protected:
  const systems::surpass::SurpassModel<C> & the_system; ///< the system whose energy will be evaluated

private:
  static utils::Logger logger;
  static const std::string name_;
};

template<typename C>
const std::string SurpassLocalRepulsionEnergy<C>::name_ = "SurpassLocalRepulsionEnergy";

template<typename C>
const core::real SurpassLocalRepulsionEnergy<C>::SURPASS_LOCAL_REPULSION_CUTOFF = 6.0;

template<typename C>
utils::Logger SurpassLocalRepulsionEnergy<C>::logger("SurpassLocalRepulsionEnergy");


} // ~ simulations
} // ~ generic
} // ~ ff
#endif
