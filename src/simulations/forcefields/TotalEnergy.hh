#ifndef SIMULATIONS_FORCEFIELDS_TotalEnergy_HH
#define SIMULATIONS_FORCEFIELDS_TotalEnergy_HH

#include <vector>
#include <memory>
#include <algorithm>
#include <iomanip>

#include <core/real.hh>
#include <core/data/basic/Array2D.hh>

#include <utils/string_utils.hh>
#include <utils/Logger.hh>

#include <simulations/forcefields/CalculateEnergyBase.hh>
#include <simulations/evaluators/Evaluator.hh>

namespace simulations {
namespace forcefields {

/** @brief Calculates energy as a weighted combination of components.
 *
 * @tparam E - the type of energy components. If want to use <code>ByResidueEnergy</code> here,
 * take just <code>TotalEnergyByResidue</code> which is the concrete instantiation with some useful extra stuff
 */
template<class E>
class TotalEnergy : public virtual CalculateEnergyBase, public evaluators::Evaluator {
public:

  /// Bare constructor
  TotalEnergy() : logger("TotalEnergy") { }

  /// Virtual destructor
  virtual ~TotalEnergy() {}

  /** @brief Adds a new energy component that will be evaluated every time this total energy is evaluated.
   *
   * Energy returned by the new term will be multiplied by the given factor.
   * @param term - energy term
   * @param factor - factor used to scale the energy value
   */
  void add_component(std::shared_ptr<E> term, const core::real factor) {

    factors.push_back(factor);
    components.push_back(term);
    sw.push_back(std::max(min_width(), core::index2(term->name().size())));
    logger << utils::LogLevel::INFO << "added  new energy component with weight = " << factor << "\n";
  }

  /// Return the name of this energy class
  virtual const std::string &name() const { return name_; }

  virtual core::index1 precision() const { return precision_; }

  virtual core::index2 min_width() const { return min_width_; }

  /** @brief Calculates the total energy of a system.
   *
   * The energy is calculated as a weighted sum of energy terms stored in this container
   */
  virtual double calculate() {
    double en = 0.0;
    for (core::index2 i = 0; i < components.size(); ++i)
      en += components[i]->calculate() * factors[i];
    return en;
  }

  /** @brief Calculates energy of a given component only
   *
   * @return <strong>unweighted</strong> value of a given energy component that participates to this total energy
   */
  virtual double calculate_component(const core::index2 which_component) const { return components[which_component]->calculate(); }

  /// This method has been inherited from Evaluator class is just a synonym for <code>calculate()</code>
  virtual core::real evaluate() { return calculate(); }

  /// Returns the number of energy terms stored in this container
  core::index2 count_components() const { return components.size(); }

  /// Returns a requested energy term
  const std::shared_ptr<E> get_component(core::index2 id) const { return components[id]; }

  /// The returned string is the header line for energy components table (printed by ostream operator)
  virtual std::string header_string() const {

    std::stringstream ss;
    ss << std::setw(sw[0]) << components[0]->name();
    for (size_t i = 1; i < components.size(); ++i)
      ss << ' ' << std::setw(sw[i]) << components[i]->name();
    ss << ' ' << std::setw(name_.size()) << name_;
    return ss.str();
  }

  /// Returns the weights used to scale energy components
  const std::vector<core::real> & get_factors() const { return TotalEnergy::factors; }

  /** @brief Provides the text field width for each energy component.
   * i-th element of the returned vector may be used to print nicely i-th energy component; just say:
   * @code
   *  std::cout <<  std::setw(get_sw()[i]) << calculate_component(i);
   * @endcode
   * @return vector of field widths
   */
  const std::vector<core::index2> & get_sw() const { return TotalEnergy::sw; }

protected:
  std::vector<core::real> factors; ///< Weights used to scale energy components
  std::vector<std::shared_ptr<E>> components; ///< Objects that evaluates particular energy types
  std::vector<core::index2> sw; ///< Width of each energy value when converted to string (used for printing energy table)

private:
  utils::Logger logger;
  static const std::string name_;
  static const size_t min_width_;
  static const core::index2 precision_;
};

template<class E>
const std::string TotalEnergy<E>::name_ = "TotalEnergy";

template<class E>
const core::index2 TotalEnergy<E>::precision_ = 2;

template<class E>
const size_t TotalEnergy<E>::min_width_ = 7;

}
}

#endif //SIMULATIONS_GENERIC_FF_TotalEnergy_HH
