#ifndef CORE_CALC_NUMERIC_InterpolatePeriodic2D_HH
#define CORE_CALC_NUMERIC_InterpolatePeriodic2D_HH

#include <vector>
#include <algorithm>

#include <core/real.hh>
#include <core/index.hh>
#include <core/data/basic/Array2D.hh>
#include <core/calc/numeric/Interpolate2D.hh>
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
class InterpolatePeriodic2D : public Interpolate2D<E, I> {
public:

  InterpolatePeriodic2D(const E min_x, const E step_x, const index2 nx, const E min_y, const E step_y, const index2 ny,
      const std::shared_ptr<core::data::basic::Array2D<E>> & z, const I & interpolator) :
      Interpolate2D<E, I>(min_x, step_x, nx, min_y, step_y, ny, z, interpolator), period_x(step_x * nx), period_y(
          step_y * ny), logger("InterpolatePeriodic2D") {
  }

  inline E operator()(E tx, E ty) {

    while (tx < Interpolate2D<E, I>::x_min) tx += period_x;
    while (tx > Interpolate2D<E, I>::x_max) tx -= period_x;
    while (ty < Interpolate2D<E, I>::y_min) ty += period_y;
    while (ty > Interpolate2D<E, I>::y_max) ty -= period_y;

    int klo_x_p = ((int) ((tx - Interpolate2D<E, I>::x_min) * Interpolate2D<E, I>::div_x)) - 1;
    const int klo_x = (klo_x_p + 1) % Interpolate2D<E, I>::n_points_x();
    const int klo_x_n = (klo_x + 1) % Interpolate2D<E, I>::n_points_x();
    const int klo_x_nn = (klo_x_n + 1) % Interpolate2D<E, I>::n_points_x();
    if (klo_x_p < 0) klo_x_p += Interpolate2D<E, I>::n_points_x();

    int klo_y_p = ((int) ((ty - Interpolate2D<E, I>::y_min) * Interpolate2D<E, I>::div_y)) - 1;
    const int klo_y = (klo_y_p + 1) % Interpolate2D<E, I>::n_points_y();
    const int klo_y_n = (klo_y + 1) % Interpolate2D<E, I>::n_points_y();
    const int klo_y_nn = (klo_y_n + 1) % Interpolate2D<E, I>::n_points_y();
    if (klo_y_p < 0) klo_y_p += Interpolate2D<E, I>::n_points_y();

    const E x_klo_x = klo_x * Interpolate2D<E, I>::step_x + Interpolate2D<E, I>::x_min;

    E mu = (tx - x_klo_x) * Interpolate2D<E, I>::div_x;
    const core::data::basic::Array2D<E> & z = *Interpolate2D<E, I>::data_values;
    core::real vx0 = Interpolate2D<E, I>::interpolator(z(klo_x_p, klo_y_p), z(klo_x, klo_y_p), z(klo_x_n, klo_y_p),
        z(klo_x_nn, klo_y_p), mu);
    core::real vx1 = Interpolate2D<E, I>::interpolator(z(klo_x_p, klo_y), z(klo_x, klo_y), z(klo_x_n, klo_y),
        z(klo_x_nn, klo_y), mu);
    core::real vx2 = Interpolate2D<E, I>::interpolator(z(klo_x_p, klo_y_n), z(klo_x, klo_y_n), z(klo_x_n, klo_y_n),
        z(klo_x_nn, klo_y_n), mu);
    core::real vx3 = Interpolate2D<E, I>::interpolator(z(klo_x_p, klo_y_nn), z(klo_x, klo_y_nn), z(klo_x_n, klo_y_nn),
        z(klo_x_nn, klo_y_nn), mu);

    const E y_klo_y = klo_y * Interpolate2D<E, I>::step_y + Interpolate2D<E, I>::y_min;
    mu = (ty - y_klo_y) * Interpolate2D<E, I>::div_y;

    return Interpolate2D<E, I>::interpolator(vx0, vx1, vx2, vx3, mu);
  }

private:
  std::string s;
  const E period_x, period_y;
  utils::Logger logger;
};

}
}
}

#endif
