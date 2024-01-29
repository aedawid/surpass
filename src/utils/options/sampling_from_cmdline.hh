#ifndef UTILS_OPTIONS_sampling_from_cmdline_HH
#define UTILS_OPTIONS_sampling_from_cmdline_HH

#include <vector>
#include <core/real.hh>

namespace utils {
namespace options {

/** @brief Returns a set of temperature values to be used in simulated annealing
 *
 * This method tests the following options:
 *     - <code>begin_temperature</code> which corresponds to <code>-sample:t_start</code> command line flag
 *     - <code>end_temperature</code> which corresponds to <code>-sample:t_end</code> command line flag
 *     - <code>temp_steps</code> which corresponds to <code>-sample:t_steps</code> command line flag
 *
 * @return a newly created vector of temperature values
 */
std::vector<core::real> annealing_temperatures_from_cmdline();

/** @brief Returns a set of temperature values to be used in simulated annealing
 *
 * This method tests the following options:
 *     - <code>begin_temperature</code> which corresponds to <code>-sample:t_start</code> command line flag
 *     - <code>end_temperature</code> which corresponds to <code>-sample:t_end</code> command line flag
 *     - <code>temp_steps</code> which corresponds to <code>-sample:t_steps</code> command line flag
 *
 * @param temperatures - a vector where the temperatures are inserted; the vector will be cleared prior filling
 * @return a reference to a vector of temperature values
 */
std::vector<core::real> &annealing_temperatures_from_cmdline(std::vector<core::real> &temperatures);

}
}

#endif
