#ifndef SIMULATIONS_CARTESIAN_EVALUATORS_RgSquare_HH
#define SIMULATIONS_CARTESIAN_EVALUATORS_RgSquare_HH

#include <memory>

#include <core/real.hh>
#include <core/calc/structural/calculate_from_structure.hh>

#include <utils/string_utils.hh>
#include <simulations/evaluators/Evaluator.hh>
#include <simulations/systems/CartesianAtomsSimple.hh>
#include <core/calc/structural/calculate_from_structure.hh>

namespace simulations {
namespace evaluators {
namespace cartesian {

/** @brief Evaluates square of the radius of gyration of a given system.
 */
template <typename C>
class RgSquare : public simulations::evaluators::Evaluator {
public:

  friend std::ostream & operator<<(std::ostream &out,const Evaluator & e);

  /** @brief set up \f$R_g^2\f$ evaluator for a given system
   * @param system - \f$R_g^2\f$ of this system will be evaluated at every <code>evaluate()</code> call
   */
  RgSquare(const systems::CartesianAtomsSimple<C> & system) : xyz(system) {}

  /// Evaluates \f$R_g^2\f$ of the system of interest
  virtual core::real evaluate() { return sqrt(core::calc::structural::calculate_Rg_square(&xyz[0], (&xyz[0]) + xyz.n_atoms)); }

  /// Returns the name of this evaluator which is "RgSquare"
  virtual const std::string & name() const { return name_; }

  /// Returns precision used to print Rg values (2 digits)
  virtual core::index1 precision() const { return 2; }

  /// Returns the width of the numeric value when printed (7 characters in total)
  virtual core::index2 min_width() const { return 7; }

  /// virtual destructor
  virtual ~RgSquare() {}

private:
  const systems::CartesianAtomsSimple<C> & xyz;
  static const std::string name_;
};

template <typename C>
const std::string RgSquare<C>::name_ = "RgSquare";

} // ~ cartesian
} // ~ evaluators
} // ~ simulations

#endif
