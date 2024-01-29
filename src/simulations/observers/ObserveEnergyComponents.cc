#include <iostream>
#include <iomanip>

#include <simulations/observers/ObserveEnergyComponents.hh>
#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/forcefields/CalculateEnergyBase.hh>

namespace simulations {
namespace observers {

template <typename E>
ObserveEnergyComponents<E>::ObserveEnergyComponents(forcefields::TotalEnergy<E> & total_energy, std::shared_ptr<std::ostream> out) :
logger("ObserveEnergyComponents"), total_energy_(total_energy), outstream(out), is_file_(false) {}

template <typename E>
bool ObserveEnergyComponents<E>::observe() {

  ++cnt;
  if(!ObserverInterface::trigger->operator()()) return false;

  *(outstream) << std::setw(5) << cnt << " ";
  double en = 0.0;
  const std::vector<core::real> & factors =  total_energy_.get_factors();
  *(outstream) << std::fixed;
  for (core::index2 i = 0; i < total_energy_.count_components(); ++i) {
    double ee = total_energy_.calculate_component(i);
    *(outstream) << ' ' << std::setw(total_energy_.get_sw()[i]) << std::setprecision(total_energy_.precision()) << ee;
    en += ee * factors[i];
  }
  *(outstream) << ' ' << std::setprecision(total_energy_.precision()) << std::setw(total_energy_.name().size()) << en << "\n";
  outstream->flush();

  return true;
}

template <typename E>
void ObserveEnergyComponents<E>::finalize() {

  if(is_file_) {
    std::shared_ptr<std::ofstream> of = std::static_pointer_cast<std::ofstream>(outstream);
    of->close();
  } else outstream->flush();
}

template <typename E>
std::string ObserveEnergyComponents<E>::header_string() const {

  return total_energy_.header_string();
}

template <typename E>
void ObserveEnergyComponents<E>::observe_header(const std::string & prefix) {
  (*outstream) << prefix << header_string() << "\n";
}

template
class ObserveEnergyComponents<simulations::forcefields::CalculateEnergyBase>;

template
class ObserveEnergyComponents<simulations::forcefields::ByResidueEnergy>;

} // ~ simulations
} // ~ observers
