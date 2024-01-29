#include <vector>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <random>

#include <core/index.hh>
#include <core/calc/statistics/Random.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

#include <simulations/movers/Mover.hh>
#include <simulations/movers/MoversSet.hh>

namespace simulations {
namespace movers {

void MoversSet::add_mover(Mover_SP mover, const size_t moves_each_step) {

  factors.push_back(moves_each_step);
  movers.push_back(mover);
  for (size_t i = 0; i < moves_each_step; ++i) sweep.push_back(mover);
  sw.push_back(std::max(min_width, core::index2(mover->name().size())));
  logger << utils::LogLevel::INFO << "added " << mover->name() << " which attempts " << moves_each_step
         << " moves at each MC step\n";
}

core::calc::statistics::RandomSequenceIterator<Mover_SP> MoversSet::begin() {

  return core::calc::statistics::RandomSequenceIterator<Mover_SP>::begin(sweep,sweep.size());
}

const std::string MoversSet::header_string() const {

  std::stringstream ss;
  ss << std::setw(sw[0]) << movers[0]->name();
  for (size_t i = 1; i < movers.size(); ++i)
    ss << ' ' << std::setw(sw[i]) << movers[i]->name();

  return ss.str();
}

utils::Logger MoversSet::logger("MoversSet");

const core::index1 MoversSet::precision = 4;
const core::index2 MoversSet::min_width = 7;

std::ostream &operator<<(std::ostream &out, const MoversSet &e) {

  out << std::fixed;
  for (size_t i = 0; i < e.movers.size(); ++i) {
    out << ' ' << std::setw(e.sw[i]) << std::setprecision(e.precision) << e.movers[i]->get_and_clear_success_rate();
  }
  return out;
}

} // ~ simulations
} // ~ movers
