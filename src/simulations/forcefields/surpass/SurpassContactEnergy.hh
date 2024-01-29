#ifndef SIMULATIONS_CARTESIAN_FF_SurpassContactEnergy_HH
#define SIMULATIONS_CARTESIAN_FF_SurpassContactEnergy_HH

#include <core/real.hh>
#include <core/index.hh>

#include <core/data/io/DataTable.hh>
#include <simulations/systems/surpass/SurpassModel.hh>
#include <simulations/forcefields/LongRangeByResidues.hh>
#include <simulations/forcefields/surpass/SurpassHydrogenBond.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

using core::real;

template<typename C>
/** @brief Contact potential of a square-well shape.
 */
class SurpassContactEnergy : public LongRangeByResidues<C> {
public:

  SurpassContactEnergy(const systems::surpass::SurpassModel<C> &system, const std::vector<std::string> &parameters) :
    LongRangeByResidues<C>(system), logger("SurpassContactEnergy"), the_system(system), HB(system) {

    if (parameters.size() < 3) {
      std::runtime_error("ContactEnergy requires exactly four parameters: high_energy_level, low_energy_level, contact_shift!");
    }
    init(utils::from_string<real>(std::string(parameters[0])), utils::from_string<real>(std::string(parameters[1])),
      utils::from_string<real>(std::string(parameters[2])));
  }

  /** @brief Define the square-well shape of energy function.
   *
   */
  SurpassContactEnergy(const systems::surpass::SurpassModel<C> &system, const real high_energy_level,
                       const real low_energy_level, const real contact_shift) : LongRangeByResidues<C>(system),
                                                   logger("SurpassContactEnergy"), the_system(system), HB(system) {

    init(high_energy_level, low_energy_level, contact_shift);
  }

  virtual ~SurpassContactEnergy() {}

  bool energy_kernel(const core::index2 moved_residue, const core::index2 the_other_residue, double &energy) {

    C o_i, o_j;
    LongRangeByResidues<C>::the_system[moved_residue].wrap(o_i);
    core::index2 atom_type_i = LongRangeByResidues<C>::the_system[moved_residue].atom_type;

    LongRangeByResidues<C>::the_system[the_other_residue].wrap(o_j);
    core::index2 atom_type_j = LongRangeByResidues<C>::the_system[the_other_residue].atom_type;
    unsigned int i1 = 0, i2 = 0;
    core::index1 id = (atom_type_i << 2) + atom_type_j; // index to the array that holds contact distances
    bool is_OK = true;

    // --- Calculate only if the two residues come from different secondary structure elements (but they area allowed to be in the same loop)
    if (the_system.ss_element_for_atoms()[moved_residue] != the_system.ss_element_for_atoms()[the_other_residue]) {
      if ((the_other_residue >= moved_residue - 4) && (the_other_residue <= moved_residue + 4)) return true;
      if ((the_system[moved_residue].atom_type == 2)||(the_system[the_other_residue].atom_type == 2)) is_OK = false;
      if (the_system[moved_residue].atom_type == the_system[the_other_residue].atom_type) {
	  if ((the_system[moved_residue].atom_type == 0)&&((the_other_residue >= moved_residue - 5) && (the_other_residue <= moved_residue + 5))) {
	    return true;
	  } else if (the_system[moved_residue].atom_type == 1) {
	    i1 = (unsigned char) the_system.beta_index_for_atoms()[moved_residue];
	    i2 = (unsigned char) the_system.beta_index_for_atoms()[the_other_residue];
            if (HB.union_find_sheets().find_set(i1) == HB.union_find_sheets().find_set(i2)) is_OK = false;	// --- condition for a strand : calculate energy only if both residues are in 2 different beta sheets
	  }
      }
      real contact_shortest_distance = contact_shift_ + contact_min_distance_[id] * 0.05;
      real contact_premium_distance = contact_ave_distance_[id];
      real contact_longest_distance = contact_max_distance_[id];
      contact_shortest_distance = contact_shortest_distance * contact_shortest_distance;
      contact_premium_distance = contact_premium_distance * contact_premium_distance;
      contact_longest_distance = contact_longest_distance * contact_longest_distance;
      double d = o_i.x - o_j.x;
      double r2 = d * d;
      if (r2 > contact_longest_distance) return true;
      d = o_i.y - o_j.y;
      r2 += d * d;
      if (r2 > contact_longest_distance) return true;
      d = o_i.z - o_j.z;
      r2 += d * d;
      if (r2 < contact_longest_distance) {
        if (r2 < contact_shortest_distance) energy += high_energy_level_;
	if ((r2 > contact_premium_distance)&&(is_OK == true)) energy += low_energy_level_;
      }
    }
    return true;
  }

  virtual const std::string &name() const { return name_; }

protected:
  real high_energy_level_;
  real low_energy_level_;
  real contact_shift_;
  core::index1 contact_min_distance_[12];
  core::index1 contact_ave_distance_[12];
  core::index1 contact_max_distance_[12];

private:
  utils::Logger logger;
  const systems::surpass::SurpassModel<C> &the_system; ///< the system whose energy will be evaluated
  const SurpassHydrogenBond <C> HB;

  void init(const real high_energy_level, const real low_energy_level, const real contact_shift);

  void load_surpass_cutoffs();

  static const std::string name_;
};

template<typename C>
void SurpassContactEnergy<C>::load_surpass_cutoffs() {
 
  core::data::io::DataTable dt;
  dt.load(core::SURPASSenvironment::from_file_or_db("forcefield/surpass_contact.dat"));
  for (const auto &row : dt) {
    core::index2 i = row.get<core::index2>(0);
    core::index2 j = row.get<core::index2>(1);
    core::index1 v = core::index1(row.get<core::real>(2) * 20);
    contact_min_distance_[(i << 2) + j] = v;
    v = core::index1(row.get<core::real>(3));
    contact_ave_distance_[(i << 2) + j] = v;
    v = core::index1(row.get<core::real>(4));
    contact_max_distance_[(i << 2) + j] = v;
  }
  if (logger.is_logable(utils::LogLevel::FINE)) {
    logger << utils::LogLevel::FINE << "Contact excluded volume parameters:\n"
           << "\t     H    E    C\n" << utils::string_format("\tH %5.2f%5.2f%5.2f\n", contact_min_distance_[0] * 0.05,
      contact_min_distance_[1] * 0.05, contact_min_distance_[2] * 0.05)
           << utils::string_format("\tE %5.2f%5.2f%5.2f\n", contact_min_distance_[4] * 0.05,
             contact_min_distance_[5] * 0.05, contact_min_distance_[6] * 0.05)
           << utils::string_format("\tC %5.2f%5.2f%5.2f\n", contact_min_distance_[8] * 0.05,
             contact_min_distance_[9] * 0.05, contact_min_distance_[10] * 0.05);
  }
}

template<typename C>
void
SurpassContactEnergy<C>::init(const real high_energy_level, const real low_energy_level, const real contact_shift) {
  high_energy_level_ = high_energy_level;
  low_energy_level_ = low_energy_level;
  contact_shift_ = contact_shift;
  logger << utils::LogLevel::INFO << "Energy parameters (high_en, low_en, shift): " << high_energy_level_
         << " " << low_energy_level_ << " " << contact_shift_ << "\n";
  LongRangeByResidues<C>::offset_ = 3;

  load_surpass_cutoffs();
}

template<typename C>
const std::string SurpassContactEnergy<C>::name_ = "SurpassContactEnergy";

}
}
}
#endif
