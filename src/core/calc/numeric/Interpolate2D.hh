#ifndef CORE_CALC_NUMERIC_Interpolate2D_HH
#define CORE_CALC_NUMERIC_Interpolate2D_HH

#include <vector>
#include <algorithm>
#include <stdexcept>

#include <core/real.hh>
#include <core/index.hh>
#include <core/data/basic/Array2D.hh>
#include <utils/Logger.hh>

namespace core {
namespace calc {
namespace numeric {

using core::real;

/** @brief Interpolates a function based on given points placed on a regular grid in 2D.
 *
 * @tparam E - the type of data being interpolated
 * @tparam I - the type of interpolator to be used
 */
template<typename E, typename I>
class Interpolate2D {
public:

  const E x_min; ///< the smallest X-argument value
  const E x_max; ///< the largest X-argument value
  const E step_x; ///< the separation between points of the interpolation grid in the X dimension - must be constant

  const E y_min; ///< the smallest Y-argument value
  const E y_max; ///< the largest Y-argument value
  const E step_y; ///< the separation between points of the interpolation grid in the Y dimension - must be constant

  Interpolate2D(const E min_x, const E step_x, const index2 nx, const E min_y, const E step_y, const index2 ny,
                const std::shared_ptr<core::data::basic::Array2D<E>> &z, const I &interpolator) :
    x_min(min_x), x_max(min_x + step_x * nx), step_x(step_x), y_min(min_y), y_max(min_y + step_y * ny), step_y(step_y),
    n_points_x_(nx), n_points_y_(ny),
    interpolator(interpolator), div_x(1.0 / step_x), div_y(1.0 / step_y), data_values(z), logger("Interpolate2D") {

    if (n_points_x_ != data_values->count_columns()) {
      logger << utils::LogLevel::CRITICAL << " the number of X points = " << n_points_x_
          << " differs from the number of rows in the given array2d = " << data_values->count_columns() << "\n";
      throw std::length_error("the number of X points differs from the number of rows in the given Array2D");
    }
    if (n_points_y_ != data_values->count_rows()) {
      logger << utils::LogLevel::CRITICAL << " the number of X points = " << n_points_y_
          << " differs from the number of rows in the given array2d = " << data_values->count_rows() << "\n";
      throw std::length_error("the number of Y points differs from the number of rows in the given Array2D");
    }
  }

  inline core::index2 n_points_x() { return Interpolate2D< E, I>::data_values->count_rows(); }

  inline core::index2 n_points_y() { return Interpolate2D< E, I>::data_values->count_columns(); }


  inline E operator()(E tx, E ty) {

    const core::data::basic::Array2D<E> & z = *data_values;
    if (tx < x_min) return z(0, 0);
    if (tx > x_max) return z(n_points_x_ - 1, 0);
    if (ty < y_min) return z(0, 0);
    if (ty > y_max) return z(0, n_points_y_ - 1);

    const int klo_x = (core::index4) ((tx - x_min) * div_x);
    const int klo_y = (core::index4) ((ty - y_min) * div_y);
    const int khi_y = klo_y + 1;
    E mu = (tx - (klo_x * step_x + x_min)) / step_x;
    core::real vx0 = interpolator(z(klo_x - 1, klo_y - 1), z(klo_x, klo_y - 1),
        z(klo_x + 1, klo_y - 1), z(klo_x + 2, klo_y - 1), mu);
    core::real vx1 = interpolator(z(klo_x - 1, klo_y), z(klo_x, klo_y),
        z(klo_x + 1, klo_y), z(klo_x + 2, klo_y), mu);
    core::real vx2 = interpolator(z(klo_x - 1, khi_y), z(klo_x, khi_y),
        z(klo_x + 1, khi_y), z(klo_x + 2, khi_y), mu);
    core::real vx3 = interpolator(z(klo_x - 1, khi_y + 1), z(klo_x, khi_y + 1),
        z(klo_x + 1, khi_y + 1), z(klo_x + 2, khi_y + 1), mu);

    mu = (ty - (klo_y * step_y + y_min)) / step_y;
    return interpolator(vx0, vx1, vx2, vx3, mu);
  }

protected:
  const core::index2 n_points_x_; ///< the number of points used for interpolation in the X direction
  const core::index2 n_points_y_; ///< the number of points used for interpolation in the Y direction
  const I & interpolator; ///< Object used to calculate interpolation in 1D
  const E div_x, div_y; ///< Inverse of bin size
  const std::shared_ptr<core::data::basic::Array2D<E>> data_values; ///< Data matrix

private:
  utils::Logger logger;
};

}
}
}

#endif
