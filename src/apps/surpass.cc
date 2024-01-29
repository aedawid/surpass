#include <cstdio>
#include <iostream>

#include <core/SURPASSenvironment.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/io/ss2_io.hh>
#include <core/data/io/Pdb.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/structural/Structure.hh>

#include <simulations/systems/surpass/SurpassModel.hh>
#include <simulations/evaluators/cartesian/CM.hh>
#include <simulations/evaluators/cartesian/RgSquare.hh>
#include <simulations/evaluators/cartesian/CrmsdEvaluator.hh>
#include <simulations/forcefields/surpass/surpass_force_field_factory.hh>
#include <simulations/movers/PerturbResidue.hh>
#include <simulations/movers/PerturbChainFragment.hh>
#include <simulations/evaluators/Timer.hh>
#include <simulations/evaluators/Evaluator.hh>
#include <simulations/evaluators/EchoEvaluator.hh>
#include <simulations/forcefields/TotalEnergyByResidue.hh>
#include <simulations/movers/Mover.hh>
#include <simulations/movers/MoversSet.hh>
#include <simulations/sampling/SimulatedAnnealing.hh>
#include <simulations/sampling/MetropolisAcceptanceCriterion.hh>
#include <simulations/observers/cartesian/PdbObserver.hh>
#include <simulations/observers/cartesian/PymolObserver.hh>
#include <simulations/observers/ObserveEnergyComponents.hh>
#include <simulations/observers/ObserveEvaluators.hh>
#include <simulations/observers/ObserveMoversAcceptance.hh>
#include <simulations/observers/TriggerLowEnergy.hh>
#include <simulations/representations/surpass_utils.hh>

#include <utils/string_utils.hh>
#include <utils/LogManager.hh>
#include <utils/options/Option.hh>
#include <utils/options/OptionParser.hh>
#include <utils/options/output_options.hh>
#include <utils/options/input_options.hh>
#include <utils/options/sampling_options.hh>
#include <utils/options/sampling_from_cmdline.hh>
#include <simulations/forcefields/ForceFieldConfig.hh>
#include <simulations/observers/ObserveReplicaFlow.hh>
#include <simulations/observers/surpass/ObserveTopologyMatrix.hh>
#include <simulations/observers/cartesian/EndVectorObserver.hh>

std::string pymol_style = R"(STYL  show spheres
show lines
select cc, name#CA
color skyblue cc
set sphere_scale 0.4 cc)";

utils::Logger logs("surpass");

using core::data::basic::Vec3;

simulations::movers::MoversSet_SP create_movers(simulations::systems::ResidueChain<Vec3> &rc,
        std::shared_ptr<simulations::forcefields::TotalEnergyByResidue> en, core::index2 which_replica) {

  using namespace utils::options;
  using namespace simulations::movers;
  using namespace simulations::systems::surpass;

  MoversSet_SP moves = std::make_shared<MoversSet>();

  std::vector<core::real> move_ranges;
  if (random_jump_range.was_used()) option_value<core::real>(random_jump_range, move_ranges);
  else move_ranges.push_back(0.5);

  std::shared_ptr<PerturbResidue<Vec3>> perturb = std::make_shared<PerturbResidue<Vec3>>(rc, *en);
  perturb->max_move_range(move_ranges[which_replica % move_ranges.size()]);
  moves->add_mover(perturb, rc.n_atoms);

  move_ranges.clear();
  if (random_n_jump_len.was_used()) {
    core::index2 n = option_value<core::index2>(random_n_jump_len);
    if (random_n_jump_range.was_used()) option_value<core::real>(random_n_jump_range, move_ranges);
    else move_ranges.push_back(0.5);
    std::shared_ptr<PerturbChainFragment<Vec3>> perturb_n = std::make_shared<PerturbChainFragment<Vec3>>(rc, n, *en);
    perturb_n->max_move_range(move_ranges[which_replica % move_ranges.size()]);
    moves->add_mover(perturb_n, rc.n_atoms / n);
  }

  return moves;
}

std::vector<core::data::structural::Structure_SP> starting_structures(
  core::data::sequence::SecondaryStructure_SP ss2_aa, core::index2 n_replicas = 1) {

  using namespace utils::options; // --- All the options are in this namespace
  using namespace core::data::structural;

  std::vector<Structure_SP> structures;
  if (input_pdb.was_used()) {
    core::data::io::Pdb reader(option_value<std::string>(input_pdb), core::data::io::is_not_alternative, true);
    reader.create_structures(structures);
  } else {
    logs << utils::LogLevel::CRITICAL << "SURPASS requires a starting conformation in PDB format\n";
    utils::exit_OK_with_message("SURPASS requires a starting conformation in PDB format\n");
  }

  while (structures.size() < n_replicas) structures.push_back(structures.back());

  for (core::index2 irepl = 0; irepl < n_replicas; ++irepl) {
    if (!simulations::representations::is_surpass_model(*structures[irepl])) {
      Structure_SP starting_structure = structures[irepl];
      core::index2 res_cnt = 0;
      for (auto res_it = starting_structure->first_residue();
           res_it != starting_structure->last_residue(); ++res_it) {
        (*res_it)->ss(ss2_aa->ss(res_cnt));
        ++res_cnt;
      }
      structures[irepl] = simulations::representations::surpass_representation(*structures[irepl]);
    }
  }
  return structures;
}

void run_annealing(core::data::structural::Structure_SP starting_structure, const simulations::forcefields::ForceFieldConfig & scoring_cfg) {

  using namespace simulations::forcefields;
  using namespace simulations::forcefields::surpass;
  using namespace simulations::movers;
  using namespace simulations::systems::surpass;
  using namespace utils::options; // --- All the options are in this namespace
  using namespace simulations::movers;
  using namespace simulations::observers;
  using namespace simulations::evaluators;

  const core::index4 n_inner_cycles = option_value<core::index2>(mc_inner_cycles, 200);
  const core::index4 n_outer_cycles = option_value<core::index2>(mc_outer_cycles, 200);
  const core::index4 cycle_size = option_value<core::index2>(mc_cycle_factor, 10);

  // --- Create the system to be sampled
  std::string input_ss2_file = option_value<std::string>(input_ss2);
  core::data::sequence::SecondaryStructure_SP ss2_aa = core::data::io::read_ss2(input_ss2_file, "");
  std::shared_ptr<SurpassModel<Vec3>> rc = std::make_shared<SurpassModel<Vec3>>(*starting_structure);

  // ---------- Prepare the scoring function ----------
  std::shared_ptr<TotalEnergyByResidue> en = create_surpass_energy<Vec3>(*rc, ss2_aa, scoring_cfg.str());

  // ---------- Movers definition ----------
  core::real move_range = (!random_jump_range.was_used()) ? 0.5 : option_value<core::real>(random_jump_range);
  logs << utils::LogLevel::INFO << "jump range: " << move_range << "\n";
  simulations::movers::MoversSet_SP movers = create_movers(*rc, en, move_range);

  // ---------- Create the sampler
  std::vector<core::real> temperatures = utils::options::annealing_temperatures_from_cmdline();
  auto sampler = simulations::sampling::SimulatedAnnealing(movers, temperatures);
  sampler.cycles(n_inner_cycles,n_outer_cycles, cycle_size);

//  auto start = std::chrono::high_resolution_clock::now();
  logs << utils::LogLevel::INFO << "Initial energy: " << en->calculate() << "\n";
  std::string out_pdb_fname = option_value<std::string>(output_pdb, "tra.pdb");
  auto tra = std::make_shared<simulations::observers::cartesian::PdbObserver<Vec3>>(*rc, *starting_structure, out_pdb_fname);

  std::shared_ptr<simulations::observers::cartesian::PdbObserver<Vec3>> min_tra = nullptr;
  if(output_pdb_min.was_used()){
    std::string fname = option_value<std::string>(output_pdb_min);
    min_tra = std::make_shared<simulations::observers::cartesian::PdbObserver<Vec3>>(*rc, *starting_structure, fname);
    core::real fraction = option_value<core::real>(output_pdb_min_fraction,0.1);
    core::real max_en = option_value<core::real>(output_pdb_min_value,en->calculate());
    std::shared_ptr<simulations::observers::TriggerLowEnergy> low_en_trigger =
        std::make_shared<simulations::observers::TriggerLowEnergy>(*en,max_en,fraction);
    min_tra->set_trigger(low_en_trigger);
  }


  // ---------- Observers & Evaluators
  ObserveEvaluators_SP stats = std::make_shared<ObserveEvaluators>("observers.dat");
  stats->add_evaluator(std::make_shared<simulations::evaluators::cartesian::RgSquare<Vec3>>(*rc));
  stats->add_evaluator(std::make_shared<simulations::evaluators::Timer>());
//    stats->add_evaluator(std::make_shared<simulations::evaluators::EchoEvaluator<core::real>>(temperature,"temperature"));
  Evaluator_SP rms = nullptr;
  if (input_pdb_native.was_used()) {
    core::data::io::Pdb native_reader(option_value<std::string>(input_pdb_native));
    core::data::structural::Structure_SP native_structure = native_reader.create_structure(0);
    rms = std::make_shared<simulations::evaluators::cartesian::CrmsdEvaluator<Vec3>>(native_structure, *rc);
  } else rms = std::make_shared<simulations::evaluators::cartesian::CrmsdEvaluator<Vec3>>(starting_structure, *rc);
  stats->add_evaluator(rms);
  stats->observe_header();
  stats->observe();

  // ---------- Create observer for energy components and movers ----------
  std::shared_ptr<ObserveEnergyComponents<ByResidueEnergy>> obs_en
    = std::make_shared<simulations::observers::ObserveEnergyComponents<ByResidueEnergy>>(*en, "energy.dat");
  obs_en->observe_header();
  ObserveMoversAcceptance_SP obs_ms = std::make_shared<simulations::observers::ObserveMoversAcceptance>(*movers, "movers.dat");
  obs_ms->observe_header();

  // --- Observer for end-to-end vector
  std::shared_ptr<simulations::observers::cartesian::EndVectorObserver<Vec3>> r_end
    = std::make_shared<simulations::observers::cartesian::EndVectorObserver<Vec3>>(*rc,"r_end.dat");

  // --- Create H-bond topology map observer
  for(core::index2 ien=0;ien<en->count_components(); ++ien) {
    std::shared_ptr<SurpassHydrogenBond<Vec3>> hb_en
      = std::dynamic_pointer_cast<SurpassHydrogenBond<Vec3>>(en->get_component(ien));
    if(hb_en!= nullptr) {
      std::shared_ptr<ObserveTopologyMatrix<Vec3>> obs_topo
        = std::make_shared<ObserveTopologyMatrix<Vec3>>(hb_en, "topology.dat");
      sampler.outer_cycle_observer(obs_topo);
    }
  }

  // --- Register all observer at the sampler
  sampler.outer_cycle_observer(stats);
  sampler.outer_cycle_observer(obs_en);
  sampler.outer_cycle_observer(obs_ms);
  sampler.outer_cycle_observer(r_end);
  sampler.outer_cycle_observer(tra);
  if (min_tra != nullptr) sampler.outer_cycle_observer(min_tra);
  sampler.run();

//  tra.finalize();
  simulations::observers::cartesian::write_pdb_conformation(*rc, *starting_structure, "final.pdb");
}

void run_replicas(std::vector<core::data::structural::Structure_SP> & starting_structures,
              const simulations::forcefields::ForceFieldConfig & scoring_cfg, std::vector<core::real> temperatures) {

  using namespace simulations::forcefields;
  using namespace simulations::forcefields::surpass;
  using namespace simulations::movers;
  using namespace simulations::systems::surpass;
  using namespace utils::options; // --- All the options are in this namespace
  using namespace simulations::movers;
  using namespace simulations::observers;
  using namespace simulations::evaluators;

  // ---------- Get parameters for the sampler and create it
  const core::index4 n_inner_cycles = option_value<core::index4>(mc_inner_cycles, 10);
  const core::index4 n_outer_cycles = option_value<core::index4>(mc_outer_cycles, 200);
  const core::index4 n_exchanges = option_value<core::index4>(replica_exchanges, 10);

  std::string input_ss2_file = option_value<std::string>(input_ss2);
  core::data::sequence::SecondaryStructure_SP ss2_aa = core::data::io::read_ss2(input_ss2_file, "");

  std::vector<std::shared_ptr<SurpassModel<Vec3>>> systems;
  std::vector<simulations::sampling::IsothermalMC_SP> replica_samplers;
  std::vector<CalculateEnergyBase_SP> energies;

  std::vector<core::real> move_ranges;
  if (random_jump_range.was_used()) option_value<core::real>(random_jump_range, move_ranges);
  else move_ranges.push_back(0.5);

  // --- Create the systems to be sampled
  for (core::index2 irepl = 0; irepl < temperatures.size(); ++irepl) {

    auto rc = std::make_shared<SurpassModel<Vec3>>(*starting_structures[irepl]);
    systems.push_back(rc);

    // ---------- Create energy function for that systems
    std::shared_ptr<TotalEnergyByResidue> en = create_surpass_energy<Vec3>(*rc, ss2_aa, scoring_cfg.str());
    energies.push_back( std::dynamic_pointer_cast<CalculateEnergyBase>(en) );

    // ---------- Movers definition ----------
    logs << utils::LogLevel::INFO << "jump range for replica " << irepl << " : "
         << move_ranges[irepl % move_ranges.size()] << "\n";
    simulations::movers::MoversSet_SP movers = create_movers(*rc, en, move_ranges[irepl % move_ranges.size()]);

    // ---------- Create the sampler
    auto sampler = std::make_shared<simulations::sampling::IsothermalMC>(movers, temperatures[irepl]);
    replica_samplers.push_back(sampler);
    sampler->cycles(n_inner_cycles,n_outer_cycles);

//    logs << utils::LogLevel::INFO << "chain length: " << rc->count_residues() << ", seq length: " << ss2_aa->length() << "\n";

//    auto start = std::chrono::high_resolution_clock::now();
    logs << utils::LogLevel::INFO << "Initial energy for replica " << irepl << " : " << en->calculate() <<
        " at temperature "<< temperatures[irepl]<<"\n";

    // ---------- Observers & Evaluators
    auto tra = std::make_shared<simulations::observers::cartesian::PdbObserver<Vec3>>(*rc, *starting_structures[irepl],
      utils::string_format("tra-%.3f.pdb",temperatures[irepl]));

    ObserveEvaluators_SP stats = std::make_shared<ObserveEvaluators>(utils::string_format("observers-%.3f.dat",temperatures[irepl]));
    stats->add_evaluator(std::make_shared<simulations::evaluators::cartesian::RgSquare<Vec3>>(*rc));
    stats->add_evaluator(std::make_shared<simulations::evaluators::Timer>());
//    stats->add_evaluator(std::make_shared<simulations::evaluators::EchoEvaluator<core::real>>(temperature,"temperature"));
    Evaluator_SP rms = nullptr;
    if (input_pdb_native.was_used()) {
      core::data::io::Pdb native_reader(option_value<std::string>(input_pdb_native));
      core::data::structural::Structure_SP native_structure = native_reader.create_structure(0);
      rms = std::make_shared<simulations::evaluators::cartesian::CrmsdEvaluator<Vec3>>(native_structure, *rc);
    } else rms = std::make_shared<simulations::evaluators::cartesian::CrmsdEvaluator<Vec3>>(starting_structures[irepl], *rc);
    stats->add_evaluator(rms);
    stats->observe_header();
//    stats->observe();

    // --- Create observer for energy components and movers
    std::shared_ptr<ObserveEnergyComponents<ByResidueEnergy>> obs_en
      = std::make_shared<simulations::observers::ObserveEnergyComponents<ByResidueEnergy>>(*en, utils::string_format("energy-%.3f.dat",temperatures[irepl]));
    obs_en->observe_header();
//    obs_en->observe();
    ObserveMoversAcceptance_SP obs_ms = std::make_shared<simulations::observers::ObserveMoversAcceptance>(*movers,
      utils::string_format("movers-%.3f.dat",temperatures[irepl]));
    obs_ms->observe_header();
//    obs_ms->observe();

    for(core::index2 ien=0;ien<en->count_components(); ++ien) {
      std::shared_ptr<SurpassHydrogenBond<Vec3>> hb_en = std::dynamic_pointer_cast<SurpassHydrogenBond<Vec3>>(en->get_component(ien));
      if(hb_en!= nullptr) {
        std::shared_ptr<ObserveTopologyMatrix<Vec3>> obs_topo = std::make_shared<ObserveTopologyMatrix<Vec3>>(hb_en,
          utils::string_format("topology-%.3f.dat",temperatures[irepl]));
        sampler->outer_cycle_observer(obs_topo);
      }
    }
    sampler->outer_cycle_observer(stats);
    sampler->outer_cycle_observer(obs_en);
    sampler->outer_cycle_observer(obs_ms);
    sampler->outer_cycle_observer(tra);
  }

  core::index2 trajectory_mode = option_value<core::index2>(utils::options::replica_observation_mode,0);
  bool replica_isothermal_observation_mode = (trajectory_mode==0);
  auto remc = std::make_shared<simulations::sampling::ReplicaExchangeMC>(replica_samplers, energies, replica_isothermal_observation_mode);
  auto remc_flow = std::make_shared<ObserveReplicaFlow>(*remc,"replica_flow.dat");
  remc->exchange_observer(remc_flow);
  remc->replica_exchanges(n_exchanges);
  remc->run();

  simulations::observers::cartesian::PdbObserver<Vec3> final(*systems[0],*starting_structures[0], "final.pdb");
  for(auto rc : systems) final.observe(*rc);
  final.finalize();
}

int main(int argc, const char *argv[]) {

  utils::LogManager::INFO();

  using namespace utils::options; // --- All the options are in this namespace

  utils::options::OptionParser &cmd = utils::options::OptionParser::get();
  cmd.register_option(utils::options::help, verbose);
  cmd.register_option(db_path, rnd_seed);
  cmd.register_option(mc_outer_cycles, mc_inner_cycles, mc_cycle_factor, random_jump_range, random_n_jump_range,
    random_n_jump_len);
  cmd.register_option(input_pdb, input_pdb_native, input_ss2);  // Input options
  cmd.register_option(output_pdb, output_pdb_min, output_pdb_min_fraction, output_pdb_min_value); // output options
  cmd.register_option(begin_temperature, end_temperature, temp_steps, replicas, replica_observation_mode, replica_exchanges);

  if (!cmd.parse_cmdline(argc, argv)) return 1;

  if (rnd_seed.was_used())
    core::calc::statistics::Random::seed(option_value<core::calc::statistics::Random::result_type>(rnd_seed));

  if (!input_ss2.was_used()) {
    logs << utils::LogLevel::SEVERE << "All-atom secondary structure must be provided with -in:ss2 command line option\n";
    return 0;
  }

  // --- Read the input secondary structure
  std::string input_ss2_file = option_value<std::string>(input_ss2);
  core::data::sequence::SecondaryStructure_SP ss2_aa = core::data::io::read_ss2(input_ss2_file);

  // --- Prepare the scoring function config
  simulations::forcefields::ForceFieldConfig scfx(utils::load_text_file(core::SURPASSenvironment::from_file_or_db("surpass.wghts", "forcefield")));
  scfx.input_ss2(input_ss2_file);

  // --- Create the sampler user requested
  if(replicas.was_used()) {
    std::vector<core::real> temperatures;
    utils::split(option_value<std::string>(replicas),temperatures,',');
    std::vector<core::data::structural::Structure_SP> starts = starting_structures(ss2_aa, temperatures.size());
    logs << utils::LogLevel::INFO << "Replica temperatures";
    for (core::real t : temperatures) logs << " " << t;
    logs << "\n";
    run_replicas(starts,scfx,temperatures);
  } else {
    core::data::structural::Structure_SP starting_structure = starting_structures(ss2_aa, 1)[0];
    run_annealing(starting_structure, scfx);
  }

}
