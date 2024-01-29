#include <core/data/basic/Vec3.hh>
#include <core/data/basic/Vec3Cubic.hh>
#include <core/calc/statistics/Random.hh>
#include <simulations/systems/BuildPolymerChain.hh>

namespace simulations {
namespace systems {

template<class C>
BuildPolymerChain<C>::BuildPolymerChain(std::unique_ptr<C[]> &system, const core::index4 n_atoms) :
  n_atoms(n_atoms), tmp(0), system(system), rand_coordinate(-1.0, 1.0),
  generator(core::calc::statistics::Random::get()) {
  box_width = C::get_box_len(); // get the size of the simulation box
}

template<class C>
bool BuildPolymerChain<C>::generate(const core::real bond_length, const core::real cutoff,
                                    const core::index2 n_bead_attempts, const core::index2 n_chain_attempts) {

  if (box_width < std::numeric_limits<core::real>::max() / 2.1) system[0].set(box_width / 2.0);
  else system[0].set(0.0);

  for (core::index2 n_tries = 0; n_tries < n_chain_attempts; ++n_tries)
    if (try_chain(bond_length, n_bead_attempts, cutoff)) return true;

  return false;
}

template<class C>
bool  BuildPolymerChain<C>::try_chain(const core::real bond_length, const core::index2 n_attempts, const core::real cutoff) {

  const core::real cutoff2 = cutoff * cutoff;
  for (core::index4 ai = 1; ai < n_atoms; ++ai) {
    core::index2 n_attmpt = 0;
    do {
      tmp.x = rand_coordinate(generator);
      tmp.y = rand_coordinate(generator);
      tmp.z = rand_coordinate(generator);
      tmp.norm(bond_length);
      tmp += system[ai - 1];
      n_attmpt++;
    } while ((!is_good_point(ai, tmp, cutoff2)) && (n_attmpt < n_attempts));
    if (n_attmpt == n_attempts) return false;
    system[ai].set(tmp);
  }
  return true;
}

template<class C>
bool BuildPolymerChain<C>::is_good_point(const core::index4 n_atoms_so_far, const C &candidate,
                                         const core::real min_distance_square) {

  for (core::index4 i = 0; i < n_atoms_so_far - 1; ++i) {
    if (candidate.closest_distance_square_to(system[i]) < min_distance_square)
      return false;
  }

  return true;
}

template class BuildPolymerChain<core::data::basic::Vec3>;
template class BuildPolymerChain<core::data::basic::Vec3Cubic>;

}
}
