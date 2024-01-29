/** \file sampling_options.hh
 * @brief defines option objects that defines sampling parameters
 */
#ifndef UTILS_OPTIONS_sampling_options_HH
#define UTILS_OPTIONS_sampling_options_HH

#include <utils/options/Option.hh>

namespace utils {
namespace options {

static Option rnd_seed("-seed", "-sample:seed", "sets random generator seed for MC sampling");

static Option md_time("-md_time", "-sample:md_time", "the MD time to run the system");
static Option mc_inner_cycles("-i_cycles", "-sample:mc_inner_cycles", "the number of small MC cycles (inner MC loop) to perform");
static Option mc_outer_cycles("-o_cycles", "-sample:mc_outer_cycles", "the number of large MC cycles (outer MC loop) to perform");
static Option mc_cycle_factor("-cycles_x", "-sample:mc_cycle_factor", "make each MC cycle N times longer");
static Option replica_exchanges("-exchanges", "-sample:exchanges", "the number of my_sampler exchanges");

static Option begin_temperature("-t_start", "-sample:t_start", "initial temperature of the simulation");
static Option end_temperature("-t_end", "-sample:t_end", "final temperature of the simulation");
static Option temp_steps("-t_steps", "-sample:t_steps", "the number of isothermal steps to make");

static Option replicas("-replicas", "-sample:replicas", "temperatures for replicas in REMC simulation (the number of temperature values defines the number of replicas)");
static Option replica_observation_mode("-observation_mode", "-sample:replicas:observation_mode",
  "observation mode: ISOTHERMAL - same temperature (default); ISOTEMPORAL - contiguous time trajectory");

static Option n_atoms("-n_atoms", "-sample:n_atoms", "the number of atoms in the sampled system");

static Option backrub_range("-sample:backrub:range", "-sample::backrub::range", "sets the maximum rotation angle [in radians] for backrub moves");
static Option phi_psi_range("-sample:phi_psi:range", "-sample::phi_psi::range", "sets the maximum rotation angle [in radians] for phi-psi moves");
static Option chain_swing_range("-sample:chain_swing:range", "-sample::chain_swing::range", "sets the maximum rotation angle [in radians] for chain-swing moves");
static Option random_jump_range("-sample:perturb:range", "-sample::perturb::range", "sets the maximum move range for a Cartesian perturbation mover");
static Option random_n_jump_range("-sample:n_perturb:range", "-sample::n_perturb::range", "sets the maximum move range for a Cartesian N-residues perturbation mover");
static Option random_n_jump_len("-sample:n_perturb:n", "-sample::n_perturb::n", "sets the number of residues (N) for a Cartesian N-residues perturbation mover");

static utils::options::Option box_size("-c", "-sample::box_size", "dimensions of the simulation box (wx,wy,wz - separated by a comma)");

}
}

#endif
