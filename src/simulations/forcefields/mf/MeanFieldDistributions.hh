/** @file MeanFieldDistributions.hh
 *  @brief Defines a set of knowledge-based potentials (defined by a spline interpolation of probability).
 *
 *  This class will turn a probability distribution into a mean field potential.
 */
#ifndef SIMULATIONS_CARTESIAN_FF_MF_MeanFieldDistributions_HH
#define SIMULATIONS_CARTESIAN_FF_MF_MeanFieldDistributions_HH

#include <cmath>

#include <iostream>
#include <memory>
#include <vector>
#include <map>

#include <core/index.hh>
#include <core/real.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/calc/numeric/Function1D.hh>

#include <utils/Logger.hh>
#include <utils/string_utils.hh>

namespace simulations {
namespace forcefields {
namespace mf {

using namespace core::calc::numeric;
using core::data::sequence::Sequence;

/// Shorter name for the type of each interpolated mean-field distribution.
typedef Function1D<core::real> EnergyComponent;
typedef std::shared_ptr<Function1D<core::real>> EnergyComponent_SP;

/** @brief Holds knowledge-based energy terms in a dictionary assigned to string keys.
 *
 * The functions held by this class are to be used to evaluate knowledge-based sequence-dependent energy.
 * Each energy term has a string key assigned. All this data is loaded from a config file which provides
 * (1) the string key to store the energy component, and (2) spline function parameters.
 *
 * Derived class should repack the energy terms stored here into by-residue indexed vector
 *
 * The example tabulates requested distribution (e.g. to make a plot of it) :
 * @include ex_MeanFieldDistributions.cc
 */
class MeanFieldDistributions {
public:

  /** @brief Creates an empty container for distributions.
   *
   * The energy functions should be loaded with <code>load_1D_distribitions()</code> method
   * @param energy_name - name of this mean-field potential; basically - an ID string
   */
  MeanFieldDistributions() : logger("MeanFieldDistributions") { }

  /// const-iterator points to the beginning of energy components stored in this energy function object
  std::map<std::string, EnergyComponent_SP>::const_iterator components_begin() const { return ff.cbegin(); }

  /// const-iterator points to the end of energy components stored in this energy function object
  std::map<std::string, EnergyComponent_SP>::const_iterator components_end() const { return ff.cend(); }

  /// Name of this energy function, as stored in the header of the input file
  const std::string & name() const { return name_; }

  /// Sets the name of this energy function
  void name(const std::string & new_name) { name_ = new_name; }

  /// Registers a new energy component; it will be assigned to a given key string
  void add_component(const std::string &key, EnergyComponent_SP distribution) {
    ff.insert(std::pair<std::string, EnergyComponent_SP>(key, distribution));
  }

  /** @brief Returns a requested energy component.
   *
   * If the component has not been registered in this container, an exception is thrown
   * @param component_key - a string identifier pointing to the energy component
   * @return a requested energy component
   */
  EnergyComponent_SP at(const std::string & component_key) { return ff.at(component_key); }

  /** @brief Returns true if a certain energy component was registered in this container.
   *
   * @param component_key - a string identifier pointing to the energy component
   * @return true if it is stored in this container; false otherwise
   */
  bool contains_distribution(const std::string & component_key) const { return ff.find(component_key) != ff.end(); }

  /** @brief Returns a vector of all keys denoting energy components (distributions) known to this container.
   *
   * @return vector of all keys mapping to distributions
   */
  const std::vector<std::string> known_distributions() const;

private:
  std::string name_; ///< Name of this energy function, as stored in the header of the input file
  std::map<std::string, EnergyComponent_SP> ff; ///< Spline functions stored for each component type (string tag)
  utils::Logger logger;
};

/** @brief Loads the components from a file.
 *
 * This method can also convert probabilities to energies. The conversion is done always when the given
 * pseudocounts fraction is non-negative. By default it <strong>is</strong> set to <code>-1</code>, so no conversion is applied
 *
 * @param ff_file - the name of file with energy functions
 * @param pseudocounts - pseudocounts fraction used to convert probabilities into energy values
 */
std::shared_ptr<MeanFieldDistributions> load_1D_distributions(const std::string & ff_file, const core::real pseudocounts_fraction = -1.0);

}
}
}

#endif

/**
 * @example ex_MeanFieldDistributions.cc
 */

