/** \file angles.hh
 * @brief evaluates planar and dihedral angles (such as \f$\Phi\f$, \f$\Psi\f$ or \f$\omega\f$) based on given coordinates
 */
#ifndef CORE_CALC_STRUCTURAL_angles_H
#define CORE_CALC_STRUCTURAL_angles_H

#include <cmath>
#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/structural/Residue.fwd.hh>

namespace core {
namespace calc {
namespace structural {

/** @brief Converts from radians to degrees.
 * @param radians - angle in radians
 */
inline static core::real to_degrees(const core::real radians) {
  return radians * 180.0 / 3.14159265;
}

/** @brief Converts from degrees to radians.
 * @param degrees - angle in degrees
 */
inline static core::real to_radians(const core::real degrees) {
  return degrees / 180.0 * 3.14159265;
}


/** @brief Converts a dihedral angle value from \f$[-180.0,180.0]\f$ degrees to \f$[0.0,360.0]\f$.
 *
 * Dihedral angles should be defined in the range \f$[-\pi,\pi]\f$ as required by standard C library. This method
 * should be used to make the angle nice-looking e.g. when plotting rotamers.
 * @param degrees - angle in degrees
 * @return transformed angle value
 */
inline static core::real start_deg_angle_from_zero(const core::real degrees) {

  return (degrees<0) ? 360.0 + degrees : degrees;
}

/** @brief Evaluates a planar angle between two vectors.
 * @param v1 - the first vector
 * @param v2 - the second vector
 * @return planar angle value
 */
template<typename T>
inline static core::real evaluate_planar_angle(const T & v1, const T &v2) {

  return acos(v1.dot_product(v2) / sqrt(v1.length() * v2.length()));
}

/** @brief Evaluates a planar angle between three points.
 * @param v1 - the first point
 * @param v2 - the second point
 * @param v3 - the third point
 * @return planar angle value
 */
template<typename T>
inline static core::real evaluate_planar_angle(const T & v1, const T &v2, const T &v3) {

  // Vector V3 --> V2
  core::real dx1 = v2.x - v3.x;
  core::real dy1 = v2.y - v3.y;
  core::real dz1 = v2.z - v3.z;
  // Vector V1 --> V2
  core::real dx2 = v2.x - v1.x;
  core::real dy2 = v2.y - v1.y;
  core::real dz2 = v2.z - v1.z;

  core::real a = dx2 * dx1 + dy2 * dy1 + dz2 * dz1;
  a /= sqrt((dx2 * dx2 + dy2 * dy2 + dz2 * dz2) * (dx1 * dx1 + dy1 * dy1 + dz1 * dz1));

  return acos(a);
}

/** @brief Evaluates a dihedral angle between four points.
 * @param v1 - the first point
 * @param v2 - the second point
 * @param v3 - the third point
 * @param v4 - the fourth point
 * @return torsion angle value
 */
template<typename T>
inline static core::real evaluate_dihedral_angle(const T & v1, const T &v2, const T &v3, const T &v4) {

  using namespace core::data::basic;

  T t1(v2); // t1 = v1 -> v2
  t1 -= v1;
  T t2(v3); // t2 = v2 -> v3
  t2 -= v2;
  T t3(v4); // t3 = v3 -> v4
  t3 -= v3;
  T n1;
  core::data::basic::cross_product(t1, t2, n1); // the first normal
  n1.norm();
  T n2;
  core::data::basic::cross_product(t2, t3, n2); // the second normal; also the u1 versor
  n2.norm();
  t2.norm(); // t2 is the u3 versor
  T u2;
  core::data::basic::cross_product(t2, n2, u2); // u2 = u3 x u1

  return -atan2(n1.dot_product(u2),n1.dot_product(n2) );
}

}
}
}
#endif
