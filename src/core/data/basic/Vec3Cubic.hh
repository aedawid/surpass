#ifndef CORE_DATA_BASIC_VEC3CUBIC_H
#define CORE_DATA_BASIC_VEC3CUBIC_H

#include <cmath>
#include <cstdlib>

#include <core/real.hh>
#include <core/data/basic/Vec3.hh>

namespace core {
namespace data {
namespace basic {

using core::real;

/** @brief A vector in 3D Cartesian space with cubic Periodic Boundary Conditions (PBC) built in
 *
 * This class behaves mostly as its base (which is Vec3). All the vector's coordinates may take any real value
 * (so they are not restrained to the simulation box) The only difference is in
 * <code>closest_distance_square_to</code> method, which takes into account  PBC while evaluating the distance.
 * The  extra features of this class are <code>wrap()</code> methods which calculate the in-the-box coordinate
 * methods used to.
 */
class Vec3Cubic: public Vec3 {
public:

  /// Return the with of the periodic box
  static inline real get_box_len() { return box_len; }

  /// Sets the new the with of the periodic box
  static inline void set_box_len(const real new_box_len) {
    box_len = new_box_len;
    box_len_half = box_len / 2.0; 
  }

  Vec3Cubic() : Vec3() {}

  Vec3Cubic(real x, real y, real z) : Vec3(x,y,z) {}

  Vec3Cubic(real v)  : Vec3(v) {}

  /// Returns the X coordinate wrapped to the box
  inline real wrap_x() const { return (x - box_len * (floor(x / box_len))); }

  /// Returns the Y coordinate wrapped to the box
  inline real wrap_y() const { return (y - box_len * (floor(y / box_len))); }

  /// Returns the Z coordinate wrapped to the box
  inline real wrap_z() const { return (z - box_len * (floor(z / box_len))); }

  /** @brief Wraps coordinates of this vector to the cubic box and stores the resulting coordinates in the given vector.
   *
   * @include apply_pbc.cc
   */
  inline void wrap(Vec3Cubic & out) const {
    out.x = wrap_x();
    out.y = wrap_y();
    out.z = wrap_z();
  }

  /** Returns the closest value of <code>this.x - v.x</code>.
   * Note, that the result might be negative!
   */
  inline real closest_delta_x(const Vec3Cubic &v) const {
    const real dx = wrap_x() - v.wrap_x();
    if (dx > 0)
      return (dx > box_len_half) ? box_len - dx : dx;
    else
      return (dx < -box_len_half) ? -box_len - dx : dx;
  }

  /** Returns the closest value of <code>this.x - v.x</code>.
   * Note, that the result might be negative!
   */
  inline real closest_delta_y(const Vec3Cubic &v) const {

    const real dy = wrap_y() - v.wrap_y();
    if (dy > 0)
      return (dy > box_len_half) ? box_len - dy : dy;
    else
      return (dy < -box_len_half) ? -box_len - dy : dy;
  }

  /** Returns the closest value of <code>this.x - v.x</code>.
   * Note, that the result might be negative!
   */
  inline real closest_delta_z(const Vec3Cubic &v) const {
    const real dz = wrap_z() - v.wrap_z();
    if (dz > 0)
      return (dz > box_len_half) ? box_len - dz : dz;
    else
      return (dz < -box_len_half) ? -box_len - dz : dz;
  }

  /** @brief Calculates the closest possible (according to PBC) square distance between a given point and this point.
   *
   * @param v - the square distance will be measured between these images of <code>v</code> and this point
   * that are located within the original simulation box
   * @return the square distance
   */
  inline real closest_distance_square_to(const Vec3Cubic & v) const {
    real r = closest_delta_x(v);
    real r2 = r * r;
    r = closest_delta_y(v);
    r2 += r * r;
    r = closest_delta_z(v);
    r2 += r * r;

    return r2;
  }

private:
  static real box_len;
  static real box_len_half;
};

} // ~ basic
} // ~ data
} // ~ core

#endif
/**
 */
