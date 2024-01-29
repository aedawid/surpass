/** \file calculate_from_structure.hh
 * @brief Simple calculations from Cartesian coordinates: \f$ R_g \f$, center of mass, etc.
 */
#ifndef CORE_CALC_STRUCTURE_calculate_from_structure_H
#define CORE_CALC_STRUCTURE_calculate_from_structure_H

#include <core/real.hh>
#include <core/index.hh>

namespace core {
namespace calc {
namespace structural {

/** @brief Calculates the square of radius of gyration from given coordinates
 *
 * @param from - begin iterator pointing to the first atom
 * @param to - end iterator pointing behind the last atom
 * @tparam It - iterator type
 * @return  \f$ R_g^2 \f$ value
 */
template<typename It>
double calculate_Rg_square(const It from,const It to) {

  double cx = 0, cy = 0, cz = 0, n = 0;
  for (auto ic = from; ic != to; ++ic) {
    ++n;
    cx += ic->x;
    cy += ic->y;
    cz += ic->z;
  }
  cx /= n;
  cy /= n;
  cz /= n;

  double s=0;
  double cc = 0;
  for (auto ic = from; ic != to; ++ic) {
    const auto &c = *ic;
    cc = c.x - cx;
    s += cc * cc;
    cc = c.y - cy;
    s += cc * cc;
    cc = c.z - cz;
    s += cc * cc;
  }

  return s / n;
}

}
}
}

#endif
