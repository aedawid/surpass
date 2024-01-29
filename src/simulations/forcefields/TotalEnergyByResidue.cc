#include <core/index.hh>

#include <simulations/forcefields/TotalEnergyByResidue.hh>

namespace simulations {
namespace forcefields {

double TotalEnergyByResidue::calculate_by_residue(const core::index2 which_residue) {
  double en = 0.0;
  for (core::index2 i = 0; i < components.size(); ++i)
    en += components[i]->calculate_by_residue(which_residue) * factors[i];
  return en;
}

double TotalEnergyByResidue::calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to) {
  double en = 0.0;
  for (core::index2 i = 0; i < components.size(); ++i)
    en += components[i]->calculate_by_chunk(chunk_from, chunk_to) * factors[i];
  return en;
}

const std::string TotalEnergyByResidue::name_ = "TotalEnergyByResidue";

std::string TotalEnergyByResidue::header_string() const {

  std::stringstream ss;
  ss << std::setw(sw[0]) << components[0]->name();
  for (size_t i = 1; i < components.size(); ++i)
    ss << ' ' << std::setw(sw[i]) << components[i]->name();
  ss << ' ' << std::setw(name_.size()) << name_;
  return ss.str();
}

std::ostream &operator<<(std::ostream &out, const TotalEnergyByResidue &e) {

  double en = 0.0;
  out << std::fixed;
  for (core::index2 i = 0; i < e.count_components(); ++i) {
    double ee = e.calculate_component(i);
    out << ' ' << std::setw(e.get_sw()[i]) << std::setprecision(e.precision()) << ee;
    en += ee * e.get_factors()[i];
  }
  out << ' ' << std::setprecision(e.precision()) << std::setw(e.name_.size()) << en;

  return out;
}

} // ~ generic
} // ~ ff
