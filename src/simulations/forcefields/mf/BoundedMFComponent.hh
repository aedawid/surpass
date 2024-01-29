#ifndef SIMULATIONS_CARTESIAN_FF_MF_BoundedMFComponent_HH
#define SIMULATIONS_CARTESIAN_FF_MF_BoundedMFComponent_HH

#include <memory>

#include <core/real.hh>
#include <core/calc/numeric/Function1D.hh>

namespace simulations {
namespace forcefields {
namespace mf {

using namespace core::calc::numeric;

/** @brief Decorates an energy term with a harmonic penalty for arguments that reach outside a certain range.
 *
 * If argument is within the given range <code>[left_bound,right_bound]</code>, the energy is calculated from a given 1D function object
 * Otherwise, a quadratic penalty is applied.
 */
class BoundedMFComponent : public Function1D<core::real> {
public:

  /** Defines an energy function based on the given function object.
   * @param energy_function - a function used to calculate energy within the range
   * @param left_bound - left bound of the range
   * @param right_bound - right bound of the range
   */
  BoundedMFComponent(const std::shared_ptr<Function1D<core::real>> energy_function, const core::real left_bound,const core::real right_bound)
    :  x_b(left_bound), x_e(right_bound), energy_function_(energy_function) { }


  /** @brief Calculates the energy, applying quadratic penalty is required
   *
   * @param x - energy argument (distance, angle, dihedral angle, etc)
   */
  virtual core::real operator()(const core::real x) const {

    if (x < x_b) return (x - x_b) * (x - x_b);
    if (x > x_e) return (x_e - x) * (x_e - x);
    return (*energy_function_)(x);
  }

private:
  core::real x_b, x_e;
  const std::shared_ptr<Function1D<core::real>> energy_function_;
};

}
}
}

#endif