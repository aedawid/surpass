#include <simulations/movers/Mover.hh>
#include <simulations/sampling/SimulatedAnnealing.hh>
#include <simulations/sampling/MetropolisAcceptanceCriterion.hh>

namespace simulations {
namespace sampling {

void SimulatedAnnealing::run() {

  for (core::index2 itemp = 0; itemp < temperatures.size(); itemp++) {
    logger << utils::LogLevel::INFO << "Temperature set to " << temperatures[itemp] << "\n";
    IsothermalMC::run(temperatures[itemp]);
  }
}

}
}
