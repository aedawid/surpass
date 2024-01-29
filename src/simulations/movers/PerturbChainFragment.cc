#include <random>
#include <memory>

#include <simulations/movers/PerturbChainFragment.hh>
#include <core/data/basic/Vec3Cubic.hh>

namespace simulations {
namespace movers {

template<class C>
PerturbChainFragment<C>::PerturbChainFragment(systems::ResidueChain<C> &system, core::index2 n_moved,
                                              forcefields::ByResidueEnergy &energy) :
  max_step_(0.5), n_moved_(n_moved), the_system(system), the_energy(energy),
  rand_bead_index(1, the_system.count_residues() - n_moved - 1), rand_coordinate(-max_step_, max_step_),
  backup(new C[the_system.n_atoms]), logger("PerturbChainFragment") {}

template<class C>
core::index2 PerturbChainFragment<C>::move(const core::index2 n_moves,
                                           simulations::sampling::AbstractAcceptanceCriterion &mc_scheme) {

  int ret = 0;
  for (int i = 0; i < n_moves; i++) if (move(mc_scheme)) ret++;
  return ret;
}

template<class C>
void PerturbChainFragment<C>::max_move_range(core::real step) {

  max_step_ = step;
  rand_coordinate = std::uniform_real_distribution<core::real>(-max_step_, max_step_);
  logger << utils::LogLevel::INFO << "Maximum move range set to [-" << -max_step_ << "," << max_step_ << "]\n";
}

template<class C>
bool PerturbChainFragment<C>::move(simulations::sampling::AbstractAcceptanceCriterion &mc_scheme) {

  last_moved_from = rand_bead_index(generator);
  last_moved_to = last_moved_from + n_moved_ - 1;
  core::real f = 2.0 / (1.0 + n_moved_);
  core::real dx = rand_coordinate(generator) * f;
  core::real dy = rand_coordinate(generator) * f;
  core::real dz = rand_coordinate(generator) * f;
  logger << utils::LogLevel::FINER << "moving the beads : " << (int) last_moved_from << " - " << (int) last_moved_to << "\n";

  core::real before = the_energy.calculate_by_chunk(last_moved_from, last_moved_to);
  for (core::index4 i = last_moved_from; i <= last_moved_to; ++i) backup[i].set(the_system.coordinates[i]);

  for (core::index4 i = 0; i < n_moved_ / 2; ++i) {
    the_system.coordinates[last_moved_from + i].x += dx * (i + 1);
    the_system.coordinates[last_moved_from + i].y += dy * (i + 1);
    the_system.coordinates[last_moved_from + i].z += dz * (i + 1);
    the_system.coordinates[last_moved_to - i].x += dx * (i + 1);
    the_system.coordinates[last_moved_to - i].y += dy * (i + 1);
    the_system.coordinates[last_moved_to - i].z += dz * (i + 1);
  }
  if (n_moved_ % 2 == 1) {
    const core::index4 ii = (last_moved_from+last_moved_to)/2;
    the_system.coordinates[ii].x += dx / f;
    the_system.coordinates[ii].y += dy / f;
    the_system.coordinates[ii].z += dz / f;
  }
  core::real after = the_energy.calculate_by_chunk(last_moved_from, last_moved_to);
  inc_move_counter();
  if (!mc_scheme.test(before, after)) {
    undo();
    if (logger.is_logable(utils::LogLevel::FINEST))
      logger << utils::LogLevel::FINEST <<
             utils::string_format("move accepted: beads %d - %d; delta(Energy): %f\n", last_moved_from, last_moved_to,
               after - before);
    return false;
  } else {
    if (logger.is_logable(utils::LogLevel::FINEST))
      logger << utils::LogLevel::FINEST
             << utils::string_format("move cancelled: beads %d - %d; delta(Energy): %f\n", last_moved_from,
               last_moved_to, after - before);
    return true;
  }
}

template<class C>
void PerturbChainFragment<C>::undo() {
  dec_move_counter();
  for (size_t i = last_moved_from; i <= last_moved_to; i++) the_system.coordinates[i].set(backup[i]);
}

template<class C>
const std::string PerturbChainFragment<C>::name_ = "PerturbChainFragment";

template class PerturbChainFragment<core::data::basic::Vec3>;
template class PerturbChainFragment<core::data::basic::Vec3Cubic>;

}
}
