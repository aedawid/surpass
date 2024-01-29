#include <thread>

#include <simulations/sampling/ReplicaExchangeMC.hh>
#include <simulations/observers/ToStreamObserver.hh>


namespace simulations {
namespace sampling {

static const std::string isothermal_enum_name = "ISOTHERMAL";
static const std::string isotemporal_enum_name = "ISOTEMPORAL";

const std::string & remc_observation_mode_name(const ReplicaExchangeObservationMode mode) {

  static const std::string ret("UNKNOWN");

  if (mode == ReplicaExchangeObservationMode::ISOTEMPORAL) return isotemporal_enum_name;
  if (mode == ReplicaExchangeObservationMode::ISOTHERMAL) return isothermal_enum_name;
  return ret;
}

ReplicaExchangeObservationMode remc_observation_mode(const std::string & mode_name) {

  if(mode_name=="ISOTHERMAL") return ReplicaExchangeObservationMode::ISOTHERMAL;
  return ReplicaExchangeObservationMode::ISOTEMPORAL;
}

void ReplicaExchangeMC::run_replica(core::index2 ireplica) {

  replicas[ireplica]->my_sampler->run();
}

void ReplicaExchangeMC::run() {

  for (core::index4 iex = 0; iex < n_exchanges; ++iex) {
/* --------- Serial variant  --------- */
//    for (core::index2 ireplica = 0; ireplica < replicas.size(); ireplica++) {
//      logs << utils::LogLevel::INFO << "Running replica " << replicas[ireplica]->replica_index_ << " in T" << ireplica
//           << " = " << temperatures_[ireplica] << "\n";
//      replicas[ireplica]->my_sampler->run(ireplica);
//    }

/* --------- Concurrent variant --------- */
    std::vector<std::thread> ths;
    for (core::index2 ireplica = 0; ireplica < replicas.size(); ireplica++) {
      ths.push_back(std::thread(&ReplicaExchangeMC::run_replica,this,ireplica));
    }
    for (auto& th : ths) th.join();

    core::index2 r = random_replica(generator);
    try_exchange(r,r+1);

    call_exchange_evaluators();
    call_exchange_observers();
  }
}

/** \brief Exchange system between two parameters' sets.
 *
 * @param l1 - the index of the first parameter set, e.g. the first temperature involved in the exchange
 * @param l2 - the index of the second parameter set, e.g. the second temperature involved in the exchange
 * @return the number of successful replica exchanges; in this case either 0 or 1
 */
bool ReplicaExchangeMC::try_exchange(const core::index2 l1, core::index2  l2) {

  // ---------- The two tasks being exchanged
  std::shared_ptr<ReplicaTask> r1 = replicas[l1];
  std::shared_ptr<ReplicaTask> r2 = replicas[l2];
  core::real delta = (1.0 / temperatures_[l1] - 1.0 / temperatures_[l2]);
  core::real deltaE = (r2->energy->calculate() - r1->energy->calculate());
  delta *= deltaE;
  if((delta<0)||(rando(generator) < exp(-delta))) {
    if(logs.is_logable(utils::LogLevel::FINE))
      logs<<utils::LogLevel::FINE << utils::string_format("Exchanging replicas %d (%.2f %.2f) with %d (%.2f %.2f)\n"
        ,l1,temperatures_[l1],r1->energy->calculate(),l2,temperatures_[l2],r2->energy->calculate());
  } else {
    if(logs.is_logable(utils::LogLevel::FINE))
      logs<<utils::LogLevel::FINE << utils::string_format("Replica exchange failed %d (%.2f %.2f) with %d (%.2f %.2f)\n"
        ,l1,temperatures_[l1],r1->energy->calculate(),l2,temperatures_[l2],r2->energy->calculate());
    return false;
  }

  // ---------- Swap the two systems in the array of tasks
  replicas[l1] = r2;
  replicas[l2] = r1;
  // ---------- Update sampler indexes
  r1->temperature_index_=l2;
  r2->temperature_index_=l1;
  // ---------- Update stats for the my_sampler walk analysis
  if (l2 == 0) r1->replica_space_flag_ = 1;
  if (l1 == 0) r2->replica_space_flag_ = 1;
  if (l2 == temperatures_.size() - 1) r1->replica_space_flag_ = 2;
  if (l1 == temperatures_.size() - 1) r2->replica_space_flag_ = 2;
  n_successful_exchanges[l1]++;
  n_successful_exchanges[l2]++;

  if (isothermal_observations) {
    for (core::index1 i = 0; i < r1->my_sampler->observe_every_inner_cycle.size(); ++i) {
      auto observer1 = std::dynamic_pointer_cast<observers::ToStreamObserver>(
        r1->my_sampler->observe_every_inner_cycle[i]);
      auto observer2 = std::dynamic_pointer_cast<observers::ToStreamObserver>(
        r2->my_sampler->observe_every_inner_cycle[i]);

      if ((observer1 != nullptr) && (observer2 != nullptr)) {
        auto stream = observer1->output_stream();
        observer1->output_stream(observer2->output_stream());
        observer2->output_stream(stream);
      }
    }

    for (core::index1 i = 0; i < r1->my_sampler->observe_every_outer_cycle.size(); ++i) {
      auto observer1 = std::dynamic_pointer_cast<observers::ToStreamObserver>(
        r1->my_sampler->observe_every_outer_cycle[i]);
      auto observer2 = std::dynamic_pointer_cast<observers::ToStreamObserver>(
        r2->my_sampler->observe_every_outer_cycle[i]);

      if ((observer1 != nullptr) && (observer2 != nullptr)) {
        auto stream = observer1->output_stream();
        observer1->output_stream(observer2->output_stream());
        observer2->output_stream(stream);
      }
    }
  }
//  r1->my_sampler->evaluate_every_inner_cycle.swap(r2->my_sampler->evaluate_every_inner_cycle);
//  r1->my_sampler->evaluate_every_outer_cycle.swap(r2->my_sampler->evaluate_every_outer_cycle);

  return 1;
}

}
}

