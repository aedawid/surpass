#ifndef CORE_CALC_NUMERIC_Interpolate1D_HH
#define CORE_CALC_NUMERIC_Interpolate1D_HH

#include <exception>
#include <core/index.hh>
#include <utils/string_utils.hh>
#include <core/calc/numeric/Function1D.hh>

namespace core {
namespace calc {
namespace numeric {

/** @brief Interpolates a function based on given points placed on a regular grid in 1D.
 *
 * @tparam D - the type of the data container that holds the points for interpolations, e.g. std::vector<E>
 * @tparam E - the type of data being interpolated
 * @tparam I - the type of interpolator to be used
 */
template<typename D, typename E, typename I>
class Interpolate1D : public Function1D<E> {
public:

  const core::index4 n_points; ///< the number of points used for interpolation
  const E x_min; ///< the smallest argument value
  const E x_max; ///< the smallest argument value

  /** @brief Spacing between the points used to create this interpolation.
   *
   * The spacing must be equal for the whole set of points. The maximum argument value covered by this interpolation
   * may be easily calculated as: <code>n_points * step + x_min</code>
   */
  const E step;

  /// const access to the data this interpolation is based on: independent argument (X)
  const D &get_x() const { return x; }

  /// const access to the data this interpolation is based on: dependent values (Y)
  const D &get_y() const { return y; }

  /** @brief Creates the interpolator based on the given points.
   *
   * @param x - a vector of x-coordinates of the points; this interpolator will make a deep copy of this data
   * @param y - a vector of y-coordinates of the points; this interpolator will make a deep copy of this data
   * @param interpolator - instance of the low-level interpolator engine to be used
   * @see CubicInterpolator LinearInterpolator CatmullRomInterpolator : currently implemented interpolators
   * defined on \f$[0,1]\f$ range that may be used by Interpolate1D
   */
  Interpolate1D(const D &x, const D &y, const I &interpolator) :
    n_points(x.size()), x_min(x[0]), x_max(x.back()), step(x[1] - x[0]), x(x), y(y), interpolator(interpolator),
    div(1.0 / step), min_allowed_value(x_min + step + step), max_allowed_value((n_points - 2) * step + x_min) {
  }

  /** @brief Calculate an interpolated value at x = arg.
   *
   * @param arg - a point for which the interpolated value will be computed; in the range from Interpolate1D::x_min
   * to x_max
   * @return interpolated value
   */
  virtual E operator()(const E arg) const {
    if (arg < min_allowed_value) return y[0];
    if (arg > max_allowed_value) return y[n_points - 1];

    const int klo = (core::index4) ((arg - x_min) * div);

    const E mu = ((arg) - x[klo]) * div;
    return interpolator(y[klo - 1], y[klo], y[klo + 1], y[klo + 2], mu);
  }

protected:
  const D x;
  const D y;
  const I &interpolator;
  const E div;
private:
  const E min_allowed_value;
  const E max_allowed_value;
};

}
}
}

#endif
