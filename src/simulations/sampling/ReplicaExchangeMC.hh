#ifndef SIMULATIONS_SAMPLING_ReplicaExchangeMC_HH
#define SIMULATIONS_SAMPLING_ReplicaExchangeMC_HH

#include <random>
#include <vector>
#include <memory>
#include <core/real.hh>
#include <core/calc/statistics/Random.hh>

#include <simulations/forcefields/CalculateEnergyBase.hh>
#include <simulations/sampling/IsothermalMC.hh>
#include <simulations/observers/ObserveReplicaFlow.fwd.hh>

namespace simulations {
namespace sampling {

/** @brief Defines the two modes of logging in ReplicaExchangeMC sampling
 *
 * Replica Exchange Monte Carlo exchanges systems between temperatures, so a simulation of a single system
 * is ISOTEMPORAL i.e. there is always the same system being observe, but the process is not isothermal.
 * This affects all the observers registered to particular MC samplers. REMC sampler can actually
 * switch observers between temperatures so they are ISOTHERMAL, but then the reaction path is broken.
 */
enum ReplicaExchangeObservationMode {

  ISOTHERMAL, ///< Ask for isothermal observations - trajectory is not contiguous in this case
  ISOTEMPORAL ///< Ask for contiguous trajectory which is not isothermal
};

/** @brief convert from ReplicaExchangeObservationMode to its name
 *
 * @param mode - ReplicaExchangeObservationMode mode enum
 * @return observation mode name
 */
const std::string & remc_observation_mode_name(const ReplicaExchangeObservationMode mode);

/** @brief convert from ReplicaExchangeObservationMode name to its enum
 *
 * @param mode - ReplicaExchangeObservationMode string ID (name)
 * @return observation mode enum
 */
ReplicaExchangeObservationMode remc_observation_mode(const std::string & mode_name);

/** @brief Replica Exchange Monte Carlo sampler.
 *
 * The sampler uses arbitrary number of IsothermalMC samplers as replicas.
 */
class ReplicaExchangeMC {
private:

  friend class simulations::observers::ObserveReplicaFlow;

  class ReplicaTask {
  public:
    core::index2 replica_index() const { return replica_index_; }
    core::index2 temperature_index() const { return temperature_index_; }
    /// 0 - no boundary hit yet; 1 or 2 - most recently hit lowest or highest temperature, respectively
    core::index2 replica_space_flag() const { return replica_space_flag_; }

    core::index2 replica_index_;
    core::index2 temperature_index_;
    core::index2 replica_space_flag_ = 0; ///< 0 - no boundary hit yet; 1 or 2 - most recently hit lowest or highest temperature, respectively
    IsothermalMC_SP my_sampler;
    forcefields::CalculateEnergyBase_SP energy;

    ReplicaTask(core::index2 id, IsothermalMC_SP r, forcefields::CalculateEnergyBase_SP e) :
      replica_index_(id), temperature_index_(id), my_sampler(r), energy(e) {}
  };

public:

  const bool isothermal_observations; ///< True if the sampler was created in ReplicaExchangeObservationMode::ISOTHERMAL mode

  ReplicaExchangeMC(std::vector<IsothermalMC_SP> & replica_samplers,
    std::vector<forcefields::CalculateEnergyBase_SP> & total_energy, bool isothermal_observations = true) :
    isothermal_observations(isothermal_observations), n_successful_exchanges(replica_samplers.size()),
    logs("ReplicaExchangeMC"), random_replica(0,replica_samplers.size() - 2) {

    if (total_energy.size() != replica_samplers.size()) {
      // throw
    }

    for (core::index2 i = 0; i < replica_samplers.size(); ++i) {
      auto r = std::make_shared<ReplicaTask>(i, replica_samplers[i], total_energy[i]);
      replicas.push_back(r);
      temperatures_.push_back(replica_samplers[i]->temperature());
    }

  }

  /// Virtual destructor
  ~ReplicaExchangeMC() { }

  /// Returns the number of replica exchange attempts performed by <code>run()</code> call
  core::index4 replica_exchanges() const { return n_exchanges; }

  const std::vector<std::shared_ptr<ReplicaTask>> get_replicas() const { return replicas; }

  /** @brief Set the number of replica exchange attempts that will be performed by <code>run()</code> call.
   * Each exchange is attempted every \f$ N_I \times N_O \f$ Monte Carlo sweeps  where  \f$ N_I \f$ and  \f$ N_O \f$
   * is the number of inner and outer cycles, respectively.
   *  @param n_exchange - the number of outer (big) cycles
   */
  void replica_exchanges(const core::index4 n_exchange) { n_exchanges = n_exchange; }

  /** Brief Runs the sampling protocol.
   * For each temperature, the method makes  \f$ N_O \f$ = <code>outer_cycles()</code> of
   * \f$ N_I \f$ = <code>inner_cycles()</code> of Monte Carlo steps.
   * The size of each MC sweep is defined within the MoversSet instance given to constructor.
   */
  virtual void run();

  /// Call all replica exchange observers
  void call_exchange_observers() { for (const auto &e : observe_every_exchange) e->observe(); }

  /// Evaluate all replica exchange evaluators
  void call_exchange_evaluators() { for (const auto &e : evaluate_every_exchange) e->evaluate(); }

  /** @brief Adds a new ObserverInterface instance to be called after every every replica exchange
   * This observer will be called  \f$ N_{ex}\f$ times
   * @param o - shared pointer to an object inheriting Observer interface
   */
  void exchange_observer(observers::ObserverInterface_SP o) { observe_every_exchange.push_back(o); }

  /** @brief Adds a new Evaluator to be called after every every replica exchange
   * This evaluator will be called  \f$ N_{ex}\f$ times
   * @param e - shared pointer to an object inheriting Evaluator interface
   */
  void exchange_evaluator(evaluators::Evaluator_SP e) { evaluate_every_exchange.push_back(e); }

  const std::vector<core::real> & temperatures() const { return temperatures_; }

private:
  std::vector<core::real> temperatures_;
  std::vector<std::shared_ptr<ReplicaTask>> replicas;
  std::vector<core::index4> n_successful_exchanges;
  core::index4 n_exchanges;
  utils::Logger logs;
  core::calc::statistics::Random &generator = core::calc::statistics::Random::get();
  std::uniform_real_distribution<core::real> rando;
  std::uniform_int_distribution<core::index2> random_replica;

  std::vector<evaluators::Evaluator_SP> evaluate_every_exchange;
  std::vector<observers::ObserverInterface_SP> observe_every_exchange;

  void run_replica(core::index2 ireplica);

  bool try_exchange(const core::index2 l1, core::index2 l2);
};


}
}

#endif // SIMULATIONS_GENERIC_SAMPLING_ReplicaExchangeMC_HH
