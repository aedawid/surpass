#ifndef CORE_CALC_STRUCTURAL_TRANSFORMATIONS_Rototranslation_H
#define CORE_CALC_STRUCTURAL_TRANSFORMATIONS_Rototranslation_H

#include <xmmintrin.h>

#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <core/calc/sse_math.hh>

#include <core/calc/structural/transformations/Rototranslation.fwd.hh>
#include <core/calc/structural/transformations/transformation_utils.hh>

using core::data::basic::Vec3;

namespace core {
namespace calc {
namespace structural {
namespace transformations {

/** \brief Rigid body transformation for vectors or atoms.
 *
 * Rototranslation can rotate and translate a point in 3D space. It may be used to move a point from one coordinate system to another, to
 * rotate a point around an axis, etc. Formally, the <code>apply(Vec3 & v)</code> method performs the following operation:
 * \f[ \vec{v} \leftarrow \mathcal{R} (\vec{v} - \vec{v}_B) + \vec{v}_A
 * \f]
 * where \f$\mathcal{R}\f$ is the rotation matrix, \f$\vec{v}_B\f$ and \f$\vec{v}_A\f$ are vectors of translation before and after rotation, respectively
 * Rotation matrix is stored <strong>row-wise</strong>
 */
class Rototranslation {
public:

  /// Default constructor initializes this object with the unit transformation (i.e. non-modifying)
  Rototranslation();

  /** \brief Applies this transformation to a given point
   *
   * Coordinates of the given vector will be changed according to this transformation.
   * @param v - vector to be translated and rotated
   */
	inline void apply(Vec3 & v) const {

		v.x -= tr_before_.x;
		v.y -= tr_before_.y;
		v.z -= tr_before_.z;
		const core::real tmpx = v.x * rot_x_.x + v.y * rot_x_.y + v.z * rot_x_.z;
		const core::real tmpy = v.x * rot_y_.x + v.y * rot_y_.y + v.z * rot_y_.z;
		const core::real tmpz = v.x * rot_z_.x + v.y * rot_z_.y + v.z * rot_z_.z;
		v.x = tmpx + tr_after_.x;
		v.y = tmpy + tr_after_.y;
		v.z = tmpz + tr_after_.z;
	}

  /** \brief Applies the inverse transformation to a given point
 *
 * Coordinates of the given vector will be changed according to this transformation.
 * @param v - vector to be translated and rotated
 */
  inline void apply_inverse(Vec3 & v) const {

	  v.x -= tr_after_.x;
	  v.y -= tr_after_.y;
	  v.z -= tr_after_.z;
	  const real tmpx = v.x * rot_x_.x + v.y * rot_y_.x + v.z * rot_z_.x;
	  const real tmpy = v.x * rot_x_.y + v.y * rot_y_.y + v.z * rot_z_.y;
	  const real tmpz = v.x * rot_x_.z + v.y * rot_y_.z + v.z * rot_z_.z;
	  v.x = tmpx + tr_before_.x;
	  v.y = tmpy + tr_before_.y;
	  v.z = tmpz + tr_before_.z;
  }

  /** \brief Applies this transformation to a given point
   *
   * Coordinates of a result will be stored in a given vector
   * @param v - vector to be translated and rotated
   * @param result - where the output coordinates are stored
   */
	inline void apply(const Vec3 & v, Vec3 & result) const  {

#ifdef __SSE2__
    __m128 vq = _mm_set_ps(0.0,v.z,v.y,v.x);
    vq = _mm_sub_ps(vq,m_tr_before);
    __m128 va = _mm_set1_ps(vq[0]);
    __m128 vc = _mm_mul_ps(va,m_rot_x);
    va = _mm_set1_ps(vq[1]);
    vc = _mm_add_ps(vc,_mm_mul_ps(va,m_rot_y));
    va = _mm_set1_ps(vq[2]);
    vc = _mm_add_ps(vc,_mm_mul_ps(va,m_rot_z));
    vc = _mm_add_ps(vc,m_tr_after);
    result.set(vc[0],vc[1],vc[2]);
#else
		result.x = v.x - tr_before_.x;
		result.y = v.y - tr_before_.y;
		result.z = v.z - tr_before_.z;
		const real tmpx = result.x * rot_x_.x + result.y * rot_x_.y + result.z * rot_x_.z;
		const real tmpy = result.x * rot_y_.x + result.y * rot_y_.y + result.z * rot_y_.z;
		const real tmpz = result.x * rot_z_.x + result.y * rot_z_.y + result.z * rot_z_.z;
		result.x = tmpx + tr_after_.x;
		result.y = tmpy + tr_after_.y;
		result.z = tmpz + tr_after_.z;
#endif
	}

  /** \brief Applies the inverse transformation to a given point
   *
   * Coordinates of the given vector will be changed according to this transformation.
   * @param v - vector to be translated and rotated
   */
  inline void apply_inverse(const Vec3 & v, Vec3 & result) const {

    result.x = v.x - tr_after_.x;
    result.y = v.y - tr_after_.y;
    result.z = v.z - tr_after_.z;
    const real tmpx = result.x * rot_x_.x + result.y * rot_y_.x + result.z * rot_z_.x;
    const real tmpy = result.x * rot_x_.y + result.y * rot_y_.y + result.z * rot_z_.y;
    const real tmpz = result.x * rot_x_.z + result.y * rot_y_.z + result.z * rot_z_.z;
    result.x = tmpx + tr_before_.x;
    result.y = tmpy + tr_before_.y;
    result.z = tmpz + tr_before_.z;
  }

  /** \brief Calculates distance between transformed query point and a template point.
   *
   * Coordinates of the <code>query</code> atom (vector) are transformed but not stored; only the distance
   * between the result and a <code>tmplt</code> point is evaluated and returned.
   * When __SSE2__ instructions are turned on, this is much faster than just to apply the transformation
   * and the to calculate the distance.
   * @param query - vector to be translated and rotated
   * @param tmplt - reference point to compute the distance
   */
	inline real distance_squared(const Vec3 & query,const Vec3 & tmplt) const {

#ifdef __SSE2__
    __m128 vq = _mm_sub_ps(_mm_set_ps(0.0,query.z,query.y,query.x),m_tr_before);
    __m128 va = _mm_set1_ps(vq[0]);
    __m128 vc = _mm_mul_ps(va,m_rot_x);
    va = _mm_set1_ps(vq[1]);
    vc = _mm_add_ps(vc,_mm_mul_ps(va,m_rot_y));
    va = _mm_set1_ps(vq[2]);
    vc = _mm_add_ps(vc,_mm_mul_ps(va,m_rot_z));
    vc = _mm_add_ps(vc,m_tr_after);
    vc = _mm_sub_ps(vc,_mm_set_ps(0.0,tmplt.z,tmplt.y,tmplt.x));
    vc = _mm_mul_ps(vc,vc);
	  return vc[0]+vc[1]+vc[2];
#else
	  Vec3 tmp;
	  apply(query,tmp);
	  return tmp.distance_square_to(tmplt);
#endif
	}

	/** \brief Replaces data for this rototranslation with data taken from the given <code>rt</code> object
	 * @param rt - source rototranslation
	 */
  inline void set(const Rototranslation & rt) {

    rot_x_.set(rt.rot_x_);
    rot_y_.set(rt.rot_y_);
    rot_z_.set(rt.rot_z_);
    tr_before_.set(rt.tr_before_);
    tr_after_.set(rt.tr_after_);
    update_mm();
  }

  /// Returns the translation vector \f$\vec{v}_A\f$ added after rotation
  inline const Vec3 & tr_after() const { return tr_after_; }

  /// Returns the translation vector \f$\vec{v}_B\f$ subtracted after rotation
  inline const Vec3 & tr_before() const { return tr_before_; }

  /// Returns the X versor of rotation (the first row of the rotation matrix)
  inline const Vec3 & rot_x() const { return rot_x_; }

  /// Returns the Y versor of rotation (the second row of the rotation matrix)
  inline const Vec3 & rot_y() const { return rot_y_; }

  /// Returns the Z versor of rotation (the third row of the rotation matrix)
  inline const Vec3 & rot_z() const { return rot_z_; }

  /** \brief Sets the new translation vector \f$\vec{v}_A\f$ that will be added after each rotation
   * @param v - the new vector \f$\vec{v}_A\f$ to be used in the transformation
   */
  inline void tr_after(const Vec3 & v) {
    tr_after_.set(v);
    update_mm();
  }

  /** \brief Sets the new translation vector \f$\vec{v}_B\f$ that will be subtracted before each rotation
   * @param v - the new vector \f$\vec{v}_B\f$ to be used in the transformation
   */
  inline void tr_before(const Vec3 & v) {
    tr_before_.set(v);
    update_mm();
  }

  /** \brief Sets the new X versor of the rotation matrix.
   * @param v - the new versor for the first row in the rotation matrix
   */
  inline void rot_x(const Vec3 & v) {
    rot_x_.set(v);
    update_mm();
  }

  /** \brief Sets the new Y versor of the rotation matrix.
   * @param v - the new versor for the second row in the rotation matrix
   */
  inline void rot_y(const Vec3 & v) {
    rot_y_.set(v);
    update_mm();
  }

  /** \brief Sets the new Z versor of the rotation matrix.
   * @param v - the new versor for the third row in the rotation matrix
   */
  inline void rot_z(const Vec3 & v) {
    rot_z_.set(v);
    update_mm();
  }

  /** \brief Sets the new translation vector \f$\vec{v}_A\f$ that will be added after each rotation
   * @param v - the new vector \f$\vec{v}_A\f$ to be used in the transformation
   */
  inline void tr_after(const real x, const real y, const real z) {
    tr_after_.set(x, y, z);
    update_mm();
  }

  /** \brief Sets the new translation vector \f$\vec{v}_B\f$ that will be subtracted before each rotation
   * @param v - the new vector \f$\vec{v}_B\f$ to be used in the transformation
   */
  inline void tr_before(const real x, const real y, const real z) {
    tr_before_.set(x, y, z);
    update_mm();
  }

  /** \brief Sets the new X versor of the rotation matrix.
   * @param x - x coordinate the new versor for the first row in the rotation matrix
   * @param y - y coordinate the new versor for the first row in the rotation matrix
   * @param z - z coordinate the new versor for the first row in the rotation matrix
   */
  inline void rot_x(const real x, const real y, const real z) {
    rot_x_.set(x, y, z);
    update_mm();
  }

  /** \brief Sets the new Y versor of the rotation matrix.
   * @param x - x coordinate the new versor for the second row in the rotation matrix
   * @param y - y coordinate the new versor for the second row in the rotation matrix
   * @param z - z coordinate the new versor for the second row in the rotation matrix
   */
  inline void rot_y(const real x, const real y, const real z) {
    rot_y_.set(x, y, z);
    update_mm();
  }

  /** \brief Sets the new Z versor of the rotation matrix.
   * @param x - x coordinate the new versor for the third row in the rotation matrix
   * @param y - y coordinate the new versor for the third row in the rotation matrix
   * @param z - z coordinate the new versor for the third row in the rotation matrix
   */
  inline void rot_z(const real x, const real y, const real z) {
    rot_z_.set(x, y, z);
    update_mm();
  }

  /** @brief Prepares a transformation that rotates points (e.g. atoms) around a given axis.
   *
   * The transformation may be used to rotate atoms around a bond. In the following example atoms are rotated
   * so to alter a \f$Phi\f$, \f$\Psi\f$ or \f$\Omega\f$ angle.
   * \include set_dihedral.cc
   *
   * @params axis_from - the first endpoint of the rotation axis, e.g. one of the two bonded atoms
   * @params axis_to - the second endpoint of the rotation axis, e.g. the second of the two bonded atoms
   * @params angle - the angle of rotation
   * @params center - the center point of rotation
   * @params destination - Rototranslation object where the transformation parameters (rotation matrix and translation vectors) will be stored
   */
  static void around_axis(const Vec3 &axis_from,
			const Vec3 &axis_to, const core::real angle, const Vec3 &center,Rototranslation & destination) {


		core::real cosi = cos(angle / 2.0);
		core::real sine = sin(angle / 2.0);
		core::real ax = axis_to.x - axis_from.x;
		core::real ay = axis_to.y - axis_from.y;
		core::real az = axis_to.z - axis_from.z;
		core::real an = sqrt(ax * ax + ay * ay + az * az);
		ax /= an;
		ay /= an;
		az /= an;
		core::real b = sine * ax;
		core::real c = sine * ay;
		core::real d = sine * az;
		core::real a2 = cosi * cosi;
		core::real b2 = b * b;
		core::real c2 = c * c;
		core::real d2 = d * d;

		destination.tr_before_.set(center);
		destination.tr_after_.set(center);
		destination.rot_x_.set(a2 + b2 - c2 - d2, 2 * b * c - 2 * cosi * d,
				2 * cosi * c + 2 * b * d);
		destination.rot_y_.set(2 * cosi * d + 2 * b * c, a2 - b2 + c2 - d2,
				2 * c * d - 2 * cosi * b);
		destination.rot_z_.set(2 * b * d - 2 * cosi * c, 2 * cosi * b + 2 * c * d,
				a2 - b2 - c2 + d2);
	}

	static Rototranslation around_axis(const Vec3 &axis_from,
			const Vec3 &axis_to, const core::real angle, const Vec3 &center) {

		Rototranslation r;
		around_axis(axis_from,axis_to,angle,center,r);

		return r;
	}

  /** @brief Prepares a transformation that rotates points (e.g. atoms) around a given axis.
   *
   * @params axis -  the rotation axis, e.g. a bond
   * @params angle - the angle of rotation
   * @params center - the center point of rotation
   * @params destination - Rototranslation object where the transformation parameters (rotation matrix and translation vectors) will be stored
   * @see around_axis(const Vec3 &, const Vec3 &, const, const Vec3 &, Rototranslation &)
   */
	static void around_axis(const Vec3 &axis,
			 const core::real angle, const Vec3 &center,Rototranslation & destination) {


		core::real cosi = cos(angle / 2.0);
		core::real sine = sin(angle / 2.0);
		core::real b = sine * axis.x;
		core::real c = sine * axis.y;
		core::real d = sine * axis.z;
		core::real a2 = cosi * cosi;
		core::real b2 = b * b;
		core::real c2 = c * c;
		core::real d2 = d * d;

		destination.tr_before_.set(center);
		destination.tr_after_.set(center);
		destination.rot_x_.set(a2 + b2 - c2 - d2, 2 * b * c - 2 * cosi * d,
				2 * cosi * c + 2 * b * d);
		destination.rot_y_.set(2 * cosi * d + 2 * b * c, a2 - b2 + c2 - d2,
				2 * c * d - 2 * cosi * b);
		destination.rot_z_.set(2 * b * d - 2 * cosi * c, 2 * cosi * b + 2 * c * d,
				a2 - b2 - c2 + d2);
	}


  /** @brief Prepares a transformation that rotates points (e.g. atoms) around a given axis.
   *
   * @params axis -  the rotation axis, e.g. a bond
   * @params angle - the angle of rotation
   * @params center - the center point of rotation
   * @see around_axis(const Vec3 &, const Vec3 &, const, const Vec3 &, Rototranslation &)
   */
  static Rototranslation around_axis(const Vec3 &axis, const core::real angle, const Vec3 &center) {

    Rototranslation r;
    around_axis(axis, angle, center, r);

    return r;
  }

  friend std::ostream &operator<<(std::ostream &out, const Rototranslation &r);

  friend void local_coordinates_three_atoms(const Vec3 &a1, const Vec3 &a2, const Vec3 &a3, Rototranslation &r);

private:
	Vec3 rot_x_, rot_y_, rot_z_;
	Vec3 tr_before_, tr_after_;
#ifdef __SSE2__
	 __m128 m_rot_x, m_rot_y,m_rot_z; ///< rotation matrix packed COLUMNWISE !!!
	 __m128 m_tr_before, m_tr_after;
#endif
protected:
	  inline void update_mm() {
	#ifdef __SSE2__
	    m_tr_before = _mm_set_ps(0.0,tr_before_.z,tr_before_.y,tr_before_.x);
	    m_rot_x = _mm_set_ps(0.0,rot_z_.x,rot_y_.x,rot_x_.x);
	    m_rot_y = _mm_set_ps(0.0,rot_z_.y,rot_y_.y,rot_x_.y);
	    m_rot_z = _mm_set_ps(0.0,rot_z_.z,rot_y_.z,rot_x_.z);
	    m_tr_after = _mm_set_ps(0.0,tr_after_.z,tr_after_.y,tr_after_.x);
	#endif
	  }
};

}
}
}
}

#endif
/**
 * \example set_dihedral.cc
 */
