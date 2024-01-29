/** @file MoversSet.hh
 *  @brief Provides MoversSet class
 */
#ifndef SIMULATIONS_GENERIC_MOVERS_MoversSet_HH
#define SIMULATIONS_GENERIC_MOVERS_MoversSet_HH

#include <vector>
#include <memory>
#include <algorithm>
#include <iomanip>

#include <core/index.hh>
#include <core/calc/statistics/RandomSequenceIterator.hh>
#include <utils/Logger.hh>

#include <simulations/movers/Mover.hh>

namespace simulations {
namespace movers {

typedef typename core::calc::statistics::RandomSequenceIterator<Mover_SP> MoversIterator;

/** @brief A set of movers which is required by any Monte Carlo sampling scheme
 *
 */
class MoversSet {
public:

  /** @brief Prints statistics about movers success rates.
   *
   * This operator prints a nice looking table with statistics for all the mover in the given set.
   * Every call of this operator results in a single row of the table, providing most recent success rates
   * The success rates are cleared by this call, so statistics for the next call are being collected from scratch (independently)
   *
   * @param out - outout stream
   * @param e - a set ov movers
   * @return reference to the stream
   */
  friend std::ostream & operator<<(std::ostream &out,const MoversSet & e);

  /** @brief Adds a new mover to the set.
   *
   * @param mover - points to a mover object
   * @param moves_each_step - how many times this mover should be called within a single MC sweep
   */
  void add_mover(Mover_SP mover, const size_t moves_each_step);

  /** @brief Starts an iteration over MC sweep
   *
   * @return iterator pointing to the very first mover (<strong>randomly selected!</strong>) of a new MC sweep
   */
  MoversIterator begin();

  /** @brief End of a Monte Carlo sweep
   *
   * @return pass-the-end iterator for a sequence of movers comprising a single MC sweep
   */
  MoversIterator end() { return core::calc::statistics::RandomSequenceIterator<Mover_SP>::end(); }

  /** @brief A nicely looking header for the movers table.
   *
   * The string returned by this method should be printed as a header for the table produced by
   * <code>operator<<(std::ostream &,const MoversSet &)</code> operator
   * @return
   */
  const std::string header_string() const;

  /// Returns the size of each sweep i.e. how many movers are called
  core::index2 sweep_size() const { return sweep.size(); }

private:
  std::vector<size_t> factors;
  std::vector<Mover_SP> sweep;
  std::vector<Mover_SP> movers;
  static utils::Logger logger;
  static const core::index1 precision;
  std::vector<core::index2> sw;
  static const core::index2 min_width;
};

std::ostream & operator<<(std::ostream &out,const MoversSet & e);

/// Defines a new type representing a shared pointer to a MoversSet object
typedef std::shared_ptr<MoversSet> MoversSet_SP;

} // ~ simulations
} // ~ movers
#endif
