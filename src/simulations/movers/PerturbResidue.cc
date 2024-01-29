#include <core/data/basic/Vec3.hh>
#include <core/data/basic/Vec3Cubic.hh>

#include <simulations/movers/PerturbResidue.hh>

namespace simulations {
namespace movers {

template<class C>
PerturbResidue<C>::PerturbResidue(systems::ResidueChain<C> &system,
                                  forcefields::ByResidueEnergy &energy) :
  max_step_(0.3), the_system(system), the_energy(energy),
  rand_residue_index(0, the_system.count_residues() - 1), rand_coordinate(-max_step_, max_step_),
  backup(new C[the_system.n_atoms]), logger("PerturbResidue") {}

template<class C>
core::index2
PerturbResidue<C>::move(const core::index2 n_moves, simulations::sampling::AbstractAcceptanceCriterion &mc_scheme) {
  int ret = 0;
  for (int i = 0; i < n_moves; i++)
    if (move(mc_scheme)) ret++;
  return ret;
}

template<class C>
void PerturbResidue<C>::max_move_range(core::real step) {
  max_step_ = step;
  rand_coordinate = std::uniform_real_distribution<core::real>(-max_step_, max_step_);
  logger << utils::LogLevel::INFO << "Maximum move range set to [-" << -max_step_ << "," << max_step_ << "]\n";
}

template<class C>
bool PerturbResidue<C>::move(simulations::sampling::AbstractAcceptanceCriterion &mc_scheme) {

  i_moved = rand_residue_index(generator);
  const systems::AtomRange<C> &last = the_system.atoms_for_residue(i_moved);
  core::real before = the_energy.calculate_by_residue(i_moved);
  for (core::index4 i = last.first_atom; i <= last.last_atom; ++i) {
    backup[i].set(the_system.coordinates[i]);
    the_system.coordinates[i].x += rand_coordinate(generator);
    the_system.coordinates[i].y += rand_coordinate(generator);
    the_system.coordinates[i].z += rand_coordinate(generator);
  }
  core::real after = the_energy.calculate_by_residue(i_moved);
  inc_move_counter();
  last_moved_from = last.first_atom;
  last_moved_to = last.last_atom;
  if (!mc_scheme.test(before, after)) {
    undo();
    if (logger.is_logable(utils::LogLevel::FINEST))
      logger << utils::LogLevel::FINEST <<
             utils::string_format("move accepted: residue %d; delta(Energy): %f\n", i_moved, after - before);
    return false;
  } else {
    if (logger.is_logable(utils::LogLevel::FINEST))
      logger << utils::LogLevel::FINEST
             << utils::string_format("move cancelled: residue %d; delta(Energy): %f\n", i_moved, after - before);
    return true;
  }
}

template<class C>
void PerturbResidue<C>::undo() {
  dec_move_counter();
  for (size_t i = last_moved_from; i <= last_moved_to; i++)
    the_system.coordinates[i].set(backup[i]);
}


template<class C>
const std::string PerturbResidue<C>::name_ = "PerturbResidue";

template class PerturbResidue<core::data::basic::Vec3>;
template class PerturbResidue<core::data::basic::Vec3Cubic>;

}
}
