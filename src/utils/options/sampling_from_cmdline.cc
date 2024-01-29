#include <core/calc/numeric/basic_math.hh>

#include <utils/options/sampling_from_cmdline.hh>
#include <utils/options/sampling_options.hh>
#include <utils/options/OptionParser.hh>

namespace utils {
namespace options {

std::vector<core::real> annealing_temperatures_from_cmdline() {

  std::vector<core::real> temperatures;
  return annealing_temperatures_from_cmdline(temperatures);
}

std::vector<core::real> &annealing_temperatures_from_cmdline(std::vector<core::real> &temperatures) {

  using namespace utils::options;

  if (!begin_temperature.was_used() && temperatures.size() > 0) return temperatures;

  temperatures.clear();
  temperatures.push_back( option_value<core::real>(begin_temperature, 1.0)); // --- we need at least one temperature value
  if (end_temperature.was_used()) {
    temperatures.push_back(option_value<core::real>(end_temperature));
    if (temp_steps.was_used()) {
      const core::index2 t_steps = option_value<core::index2>(temp_steps, 2);
      core::calc::numeric::evenly_spaced_values(temperatures[0], temperatures[1], t_steps, temperatures);
    }
  }

  return temperatures;
}

}
}
