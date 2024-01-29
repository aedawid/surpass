#include <simulations/movers/Mover.hh>
#include <simulations/sampling/IsothermalMC.hh>
#include <simulations/sampling/MetropolisAcceptanceCriterion.hh>

namespace simulations {
namespace sampling {

void IsothermalMC::run() {

  MetropolisAcceptanceCriterion mc(temperature_);

  for (core::index4 i = 0; i < n_outer_cycles; i++) {
    for (core::index2 j = 0; j < n_inner_cycles; j++) {
      for (core::index4 k = 0; k < n_cycle_size; ++k) {
        for (movers::MoversIterator m_it = movers->begin(); m_it != movers->end(); ++m_it) (*m_it)->move(mc);
      }
      call_inner_cycle_evaluators();
      call_inner_cycle_observers();
    }
    call_outer_cycle_evaluators();
    call_outer_cycle_observers();
  }
}

}
}
