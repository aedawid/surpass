#ifndef SIMULATIONS_OBSERVERS_ObserveReplicaFlow_HH
#define SIMULATIONS_OBSERVERS_ObserveReplicaFlow_HH

#include <simulations/observers/ObserverInterface.hh>
#include <simulations/sampling/ReplicaExchangeMC.hh>
#include <simulations/evaluators/Timer.hh>

namespace simulations {
namespace observers {

class ObserveReplicaFlow : public ObserverInterface {
public:

  ObserveReplicaFlow(const sampling::ReplicaExchangeMC & replicas, std::string file_name);

  virtual bool observe();

  /** @brief Does nothing, the file is always kept closed.
   */
  virtual void finalize() {}

  virtual core::index4 count_observe_calls() const { return cnt; }

private:
  std::string fname;
  const sampling::ReplicaExchangeMC & replicas_;
  evaluators::Timer timer;
  core::index4 cnt = 0;
};

} // ~ sampling
} // ~ simulations

#endif
