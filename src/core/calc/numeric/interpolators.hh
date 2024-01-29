/** \file interpolators.hh
 * @brief Provides low-level routines for data interpolation
 */
#ifndef CORE_CALC_NUMERIC_interpolators_HH
#define CORE_CALC_NUMERIC_interpolators_HH

#include <core/index.hh>

namespace core {
namespace calc {
namespace numeric {

/** @brief linear interpolation on irregularly spaced data points.
 *
 * @tparam D - the type of  data container, e.g. std::vector
 * @tparam E - the type of  data elements, e.g. core::real, double, etc.
 * @param x_sorted - x values, sorted
 * @param y_values - y values, evaluated for the given <code>x_sorted</code> points
 * @param arg - point where the function will be interpolated
 * @return interpolated value
 */
template<typename D, typename E>
E linear_irregular_interpolator(D x_sorted, D y_values, E arg) {

  if (x_sorted[0] >= arg) return y_values[0];

  if (x_sorted.back() <= arg) return y_values.back();

  core::index4 klo = 0;
  core::index4 khi = x_sorted.size() - 1;
  while (khi - klo > 1) {
    core::index4 k = (khi + klo) >> 1;
    if (x_sorted[k] > arg) khi = k;
    else klo = k;
  }

  E a = (y_values[khi] - y_values[klo]) / (x_sorted[khi] - x_sorted[klo]);
  E b = -a * x_sorted[klo] + y_values[klo];

  return a * arg + b;
}

/** @brief linear interpolation on the range \f$[0,1]\f$
 *
 * Returns the value corresponding to a point lying on the range \f$[0,1]\f$. The value is computed from
 * a linear interpolation based on the two points: \f$(0,y_1)\f$ and \f$(0,y_2)\f$.
 * User has to properly convert his x value to the <code>mu</code> argument.
 *
 * @tparam E - the type of  data elements, e.g. core::real, double, etc.
 * @param y0 - unused by the method but necessary for consistence with other interpolators
 * @param y1 - the value at the left-hand side (one that corresponds to x=0 of the interpolation range) point
 * @param y1 - the value at the right-hand side (one that corresponds to x=1 of the interpolation range) point
 * @param y3 - unused by the method but necessary for consistence with other interpolators
 * @param mu - argument as a fraction of the range over which the function is interpolated,
 * must be in the range of \f$[0,1]\f$
 * @return interpolated value
 */
template<typename E>
struct LinearInterpolator {
  E operator()(const E y0, const E y1, const E y2, const E y3, const E mu) const {

    return (y1 * (1 - mu) + y2 * mu);
  }
};

/** @brief Cubic interpolation on the range \f$[0,1]\f$
 *
 * Returns the value corresponding to a point lying on the range \f$[0,1]\f$. The value is computed from
 * a cubic interpolation based on the four points: \f$(-1,y_0)\f$, \f$(0,y_1)\f$, \f$(1,y_2)\f$ and \f$(2,y_3)\f$.
 * User has to properly convert his x value to the <code>mu</code> argument.
 *
 * @tparam E - the type of  data elements, e.g. core::real, double, etc.
 * @param y0 - the value at the point preceding \f$(0,y_1)\f$
 * @param y1 - the value at the left-hand side (one that corresponds to x=0 of the interpolation range) point
 * @param y2 - the value at the right-hand side (one that corresponds to x=1 of the interpolation range) point
 * @param y3 - the value at the point following \f$(1,y_2)\f$
 * @param mu - argument as a fraction of the range over which the function is interpolated,
 * must be in the range of \f$[0,1]\f$
 * @return interpolated value
 */
template<typename E>
struct CubicInterpolator {
  E operator()(const E y0, const E y1, const E y2, const E y3, const E mu) const {

    const E mu2 = mu * mu;
    const E a0 = y3 - y2 - y0 + y1;
    const E a1 = y0 - y1 - a0;
    const E a2 = y2 - y0;
    const E a3 = y1;

    return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
  }
};

/** @brief Cubic Catmull-Rom interpolation on the range \f$[0,1]\f$ provides much smoother curve than the regular cubic variant.
 *
 * Returns the value corresponding to a point lying on the range \f$[0,1]\f$. The value is computed from
 * a Catmull-Rom cubic interpolation based on the four points: \f$(-1,y_0)\f$, \f$(0,y_1)\f$, \f$(1,y_2)\f$ and \f$(2,y_3)\f$.
 * User has to properly convert his x value to the <code>mu</code> argument.
 *
 * @tparam E - the type of  data elements, e.g. core::real, double, etc.
 * @param y0 - the value at the point preceding \f$(0,y_1)\f$
 * @param y1 - the value at the left-hand side (one that corresponds to x=0 of the interpolation range) point
 * @param y2 - the value at the right-hand side (one that corresponds to x=1 of the interpolation range) point
 * @param y3 - the value at the point following \f$(1,y_2)\f$
 * @param mu - argument as a fraction of the range over which the function is interpolated,
 * must be in the range of \f$[0,1]\f$
 * @return interpolated value
 */
template<typename E>
struct CatmullRomInterpolator {

  E operator()(const E y0, const E y1, const E y2, const E y3, const E mu) const {

    const E mu2 = mu * mu;
    const E a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
    const E a1 = y0 - 2.5 * y1 + 2 * y2 - 0.5 * y3;
    const E a2 = -0.5 * y0 + 0.5 * y2;
    const E a3 = y1;

    return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
  }
};

}
}
}

#endif
