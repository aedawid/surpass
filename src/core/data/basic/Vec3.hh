/** @file Provides Vec3 data type - a vector in 3D
 */
#ifndef CORE_DATA_BASIC_VEC3_H
#define CORE_DATA_BASIC_VEC3_H

#include <cmath>
#include <limits>
#include <functional>
#include <xmmintrin.h>

#include <iostream>
#include <vector>
#include <memory>

#include <core/data/basic/Vec3.fwd.hh>
#include <core/calc/numeric/basic_math.hh>
#include <core/real.hh>
#include <core/index.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace basic {

/** \brief Represents a vector in the 3D space.
 */
class Vec3 {
public:

  /** @name X, Y, Z coordinates of the vector
   */
  ///@{
  core::real x, y, z;
  ///@}

  char chain_id = 'A';    ///< Character denoting a chain this atom belongs to
  core::index1 residue_type = 0; ///< A place to keep information about residue type
  /// General purpose register; but primarily padding to ensure even placement in memory.
  core::index2 register_ = 0;
  core::index2 atom_type = 0; ///< A place to keep information about atom type
  core::index2 residue_index = 0; ///< A place to keep the order number of a  residue this atom belongs to

  /// Constructor creates a new vector placed at (0,0,0)
  Vec3() : x(0), y(0), z(0) { }

  /// Constructor creates a new vector placed at (x,y,z)
  Vec3(real x, real y, real z) : x(x), y(y), z(z) { }

  /// Constructor creates a new vector placed at (v[0],v[1],v[2])
  Vec3(const std::vector<float> & v) : x(v[0]), y(v[1]), z(v[2]) { }

  /// Constructor creates a new vector placed at (v[0],v[1],v[2])
  Vec3(const std::vector<double> & v) : x(v[0]), y(v[1]), z(v[2]) { }

  /// Constructor creates a new vector placed at (v,v,v)
  Vec3(real v) : x(v), y(v), z(v) { }

  /** @brief Indexing operator provides X, Y or Z coordinate based on the index.
   * @param axis - index of the coordinate: 0, 1 or 2 to get X, Y or Z, respectively
   */
  real operator[](const index1 axis) const {

    if(axis==0) return x;
    if(axis==1) return y;
    return z;
  }

  /// Add a constant to the coordinates of this vector
  inline Vec3& operator+=(const Vec3& r);

  /// Subtract a constant from the coordinates of this vector
  inline Vec3& operator-=(const Vec3& r);

  /// Divide the coordinates of this vector by a core::real constant
  inline Vec3& operator/=(const real r);

  /// Multiply the coordinates of this vector by a core::real constant
  inline Vec3& operator*=(const real r);

  inline bool operator==(const Vec3& v) const;

  void set(const real v) { x = y = z = v; }

  void set(const real ax, const real ay, const real az) {
    x = ax;
    y = ay;
    z = az;
  }

  /// Copy coordinates of a given Vec3 object to this vector
  void set(const Vec3& v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }

  /// Copy the first three coordinates of a given std::vector to this vector
  void set(const std::vector<core::real>& v) {
    x = v[0];
    y = v[1];
    z = v[2];
  }

  /** @brief Returns the size of a simulation box.
   *
   * Vec3 class uses open boundary conditions so there is no box at all. This method returns
   * <code>std::numeric_limits<core::real>::max()</code> value
   * @return <code>std::numeric_limits<core::real>::max()</code> value
   */
  static inline real get_box_len() { return std::numeric_limits<real>::max(); }

  /** @brief Sets the new value of the size of a simulation box.
   *
   * Vec3 class uses open boundary conditions so there is no box at all. This method does nothing
   */
  static inline void set_box_len(real new_box_len) {}

  inline void norm() { *this /= length(); }

  inline void norm(real new_length) { *this /= (length() / new_length); }

  inline real dot_product(const Vec3& r) const { return x * r.x + y * r.y + z * r.z; }

  inline real length_squared() const { return x * x + y * y + z * z; }

  inline real length() const { return sqrt(length_squared()); }

//  void round_self() {
//    x = round(x);
//    y = round(y);
//    z = round(z);
//  }

  /** @brief Calculates the square distance between a given point and this point.
   *
   * This method first calculates the squared distance of X coordinates. If this is already higher
   * than the given cutoff value, this method returns the cutoff without further distance evaluation.
   * @param v - the square distance will be measured between <code>v</code> and this point
   * @param cutoff2 - cutoff for the squared distance value.
   * @return the square distance
   */
  inline real distance_square_to(const Vec3& v,real cutoff2) const {
    real r = v.x - x;
    real r2 = r * r;
    if (r2 >= cutoff2) return cutoff2;
    r = v.y - y;
    r2 += r * r;
    if (r2 >= cutoff2) return cutoff2;
    r = v.z - z;
    r2 += r * r;

    return r2;
  }

  /** @brief Calculates the square distance between a given point and this point without applying periodic boundary conditions.
   * @param v - the square distance will be measured between <code>v</code> and this point
   * @return the square distance
   */
  inline real distance_square_to(const Vec3& v) const {
    real r = v.x - x;
    real r2 = r * r;
    r = v.y - y;
    r2 += r * r;
    r = v.z - z;
    r2 += r * r;

    return r2;
  }

  /** @brief Calculates the closest possible (according to PBC) square distance between a given point and this point.
   *
   * This method is identical to <code>distance_square_to()</code> as Vec3 class does not implement periodic boundary conditions.
   * It is however overridden by the derived class Vec3Cubic which provides PBC.
   * @param v - the square distance will be measured between <code>v</code> and this point
   * @return the square distance
   */
  inline real closest_distance_square_to(const Vec3 & v) const { return distance_square_to(v); }

  /** @brief Calculates the distance between a given point and this point without applying periodic boundary conditions.
   * @param v - the distance will be measured between <code>v</code> and this point
   * @return the distance
   */
  inline real distance_to(const Vec3& v) const {
    real r = v.x - x;
    real r2 = r * r;
    r = v.y - y;
    r2 += r * r;
    r = v.z - z;
    r2 += r * r;

    return sqrt(r2);
  }

  inline real wrap_x() const { return x; }

  inline real wrap_y() const { return y; }

  inline real wrap_z() const { return z; }

  /** Returns the closest value of <code>this.x - v.x</code>.
   * Note, that the result might be negative!
   */
  inline real closest_delta_x(const Vec3 &v) const { return x - v.x; }

  /** Returns the closest value of <code>this.y - v.y</code>.
   * Note, that the result might be negative!
   */
  inline real closest_delta_y(const Vec3 &v) const { return y - v.y; }

  /** Returns the closest value of <code>this.z - v.z</code>.
   * Note, that the result might be negative!
   */
  inline real closest_delta_z(const Vec3 &v) const { return z - v.z; }

  inline void wrap(Vec3 & out) const {
    out.x = wrap_x();
    out.y = wrap_y();
    out.z = wrap_z();
  }

  friend std::ostream& operator<<(std::ostream &out, const Vec3 &v);
  friend utils::Logger &operator <<(utils::Logger &logger, const Vec3 & coordinates);
};

inline Vec3& Vec3::operator+=(const Vec3& r) {

  this->x += r.x;
  this->y += r.y;
  this->z += r.z;
  return *this;
}

inline bool Vec3::operator==(const Vec3& v) const {

  return ((this->x == v.x) && (this->y == v.y) && (this->z == v.z));
}

inline Vec3& Vec3::operator-=(const Vec3& r) {

  this->x -= r.x;
  this->y -= r.y;
  this->z -= r.z;
  return *this;
}

inline Vec3& Vec3::operator/=(const real r) {

  this->x /= r;
  this->y /= r;
  this->z /= r;
  return *this;
}

inline Vec3& Vec3::operator*=(const real r) {

  this->x *= r;
  this->y *= r;
  this->z *= r;
  return *this;
}

inline Vec3 operator*(real r, Vec3 v) {
  v *= r;
  return v;
}

inline Vec3 operator*(Vec3 v, real r) {
  v *= r;
  return v;
}

inline Vec3 operator+(Vec3 v, const Vec3& r) {
  v += r;
  return v;
}

inline Vec3 operator-(Vec3 v, const Vec3& r) {
  v -= r;
  return v;
}

inline real operator*(const Vec3& l, const Vec3& r) {
  real dot_product = l.x*r.x + l.y*r.y + l.z*r.z;
  return dot_product;    
}

static inline void cross_product(const Vec3& l, const Vec3& r, Vec3 &result) {

  result.x = l.y * r.z - l.z * r.y;
  result.y = l.z * r.x - l.x * r.z;
  result.z = l.x * r.y - l.y * r.x;
}

/** @brief Calculates determinant of a matrix defined by its row vectors.
 *
 * @param v1 - the first row of the input matrix
 * @param v2 - the second row of the input matrix
 * @param v3 - the third row of the input matrix
 * @return determinant of a 3x3 matrix
 */
static inline real det(const Vec3& v1, const Vec3& v2, const Vec3 &v3) {

  return v1.x*((v2.y*v3.z) - (v3.y*v2.z)) -v1.y*(v2.x*v3.z - v3.x*v2.z) + v1.z*(v2.x*v3.y - v3.x*v2.y);
}

/** @brief Calculates signed distance between the first and fourth point.
 *
 * The signed distance is defined as:
 * \f[
 * r_{14}^{*} =  \textrm{sign}(
 * \det \left( \begin{array}{c}
 * p_1 - p_0 \\
 * p_2 - p_1 \\
 * p_3 - p_2 \end{array} \right)
 * ) \times \|p_3-p_0\|
 * \f]
 * @param p0 - the first atom (point)
 * @param p1 - the second atom (point)
 * @param p2 - the third atom
 * @param p3 - the fourth atom (point)
 * @return determinant of a 3x3 matrix
 */
static inline real r14x(const Vec3& p0,const Vec3& p1, const Vec3& p2, const Vec3 &p3) {

  Vec3 v1(p1);
  v1-=p0;
  Vec3 v2(p2);
  v2-=p1;
  Vec3 v3(p3);
  v3-=p2;
  return core::calc::numeric::sgn(det(v1, v2, v3)) * p3.distance_to(p0);
}

static inline __m128 cross_product_sse(__m128 & a, __m128 & b) {

  __m128 result = _mm_sub_ps(_mm_mul_ps(b, _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1))),
      _mm_mul_ps(a, _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1))));
  return _mm_shuffle_ps(result, result, _MM_SHUFFLE(3, 0, 2, 1));
}

static inline void cross_product_sse(const Vec3& l, const Vec3& r, Vec3 &result) {

  __m128 a = _mm_set_ps(0.0, l.x, l.y, l.z);
  __m128 b = _mm_set_ps(0.0, r.x, r.y, r.z);
  __m128 o = cross_product_sse(a, b);

  float res[4] = { 0, 0, 0, 0 };
  _mm_storer_ps(res, o);
  result.x = res[1];
  result.y = res[2];
  result.z = res[3];
}

/** @brief Less-than comparison between Vec3 objects is implemented as comparing their X coordinates
 *  @param ai -  the first vector being compared
 *  @param aj -  the second vector being compared
 */
static inline bool operator<(const Vec3 & ai,const Vec3 & aj)  {

  if (ai.x != aj.x) return ai.x < aj.x;
  if (ai.y != aj.y) return ai.y < aj.y;
  if (ai.z != aj.z) return ai.z < aj.z;
  return false;
}

} // ~ basic
} // ~ data
} // ~ core

namespace std {

/** @brief Calculates a hash of a Vec3 object from its coordinates
 */
template<>
struct hash<core::data::basic::Vec3> {
  std::size_t operator()(const core::data::basic::Vec3 &k) const {
    return (int(k.x * 10) + (int(k.y * 10) << 10) + (int(k.x * 10) << 20));
  }
};
}

#endif
