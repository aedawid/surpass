#include <iomanip>
#include <fstream>

#include <simulations/observers/ObserveReplicaFlow.hh>

namespace simulations {
namespace observers {

ObserveReplicaFlow::ObserveReplicaFlow(const sampling::ReplicaExchangeMC & replicas, std::string file_name) :
  fname(file_name), replicas_(replicas) {
  std::ofstream out(fname);
  out.close();
}

/// This method is called to take observations
bool ObserveReplicaFlow::observe() {

  if(!ObserverInterface::trigger->operator()()) return false;

  core::index1 sw = log10(double(replicas_.temperatures().size())) + 1;
  std::ofstream out(fname, std::fstream::out | std::fstream::app);
  out << std::fixed << std::showpoint << std::setw(timer.min_width())<< std::setprecision(int(timer.precision())) << timer.evaluate()<<"   ";
  for(core::index2 i=0;i<replicas_.temperatures().size();++i)
    out << std::setw(sw) << replicas_.replicas[i]->replica_index_ << " ";
  out << "  ";

  for(core::index2 i=0;i<replicas_.temperatures().size();++i)
    out << std::setw(1) << replicas_.replicas[i]->replica_space_flag_ << " ";
  out << "  ";

  for(core::index2 i=0;i<replicas_.temperatures().size();++i)
    out << std::setw(4) << replicas_.n_successful_exchanges[i] << " ";

  out << "\n";
  out.close();

  return true;
}

} // ~ observers
} // ~ simulations
