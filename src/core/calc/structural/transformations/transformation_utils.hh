/** \file transformation_utils.hh
 * @brief Utility methods that create various rototranslation transformations
 */
#ifndef CORE_CALC_STRUCTURAL_TRANSFORMATIONS_transformation_utils_H
#define CORE_CALC_STRUCTURAL_TRANSFORMATIONS_transformation_utils_H

#include <core/real.hh>
#include <core/data/basic/Vec3.hh>

#include <core/calc/structural/transformations/Rototranslation.fwd.hh>
#include <core/data/structural/Residue.hh>

using core::data::basic::Vec3;

namespace core {
namespace calc {
namespace structural {
namespace transformations {

/** \brief Creates a rototranslation that moves points into a local coordinate system.
 *
 * The local coordinate system is defined by three atoms (points, vectors) as follows:
 *   - \f$ \vec{v}_z  = \vec{v}_{2,3} -  \vec{v}_{1,2}\f$
 *   - \f$ \vec{v}_x  = \vec{v}_{2,3} +  \vec{v}_{1,2}\f$
 *   - \f$ \vec{v}_y  = \vec{v}_z \times   \vec{v}_x\f$
 * where:
 *   - \f$ \vec{v}_{1,2} = |\vec{2} - \vec{1}| \f$
 *   - \f$ \vec{v}_{2,3} = |\vec{3} - \vec{2}| \f$
 *
 * The example program prints local coordinates of all atoms from a residue defined in its local coordinate system
 *
 * @param a1 - the first point
 * @param a2 - the second point
 * @param a3 - the third point
 * \include ex_local_coordinates_three_atoms.cc
 */
Rototranslation_SP local_coordinates_three_atoms(const Vec3 & a1,const Vec3 & a2,const Vec3 & a3);

/** \brief Creates a rototranslation that moves points into a local coordinate system.
 *
 * Looks for three atoms of given names in the given residue and calls <code>local_coordinates_three_atoms(const Vec3 & ,const Vec3 & ,const Vec3 & );</code>
 * @param r - a residue containing the three atoms of interest
 * @param atom_names - vector of three names of atoms
 * @returns a rototranslation that transform from global to local coordinate system
 */
Rototranslation_SP local_coordinates_three_atoms(const core::data::structural::Residue &r,const std::vector<std::string> & atom_names);

/** \brief Creates a rototranslation that moves points into a local coordinate system.
 *
 * Works like <code>local_coordinates_three_atoms(const Vec3 & a1,const Vec3 & a2,const Vec3 & a3)</code>, but result is stored in the
 * given Rototranslation object and nothing new  is allocated.
 */
void local_coordinates_three_atoms(const Vec3 &a1, const Vec3 &a2, const Vec3 &a3, Rototranslation & r);

/** \brief Creates a rototranslation that moves points into a local coordinate system defined for a base pair.
 *
 * Local reference frame for a pair of paired bases according to <br/>
 * Olson et al. <em>A Standard Reference Frame for the Description of Nucleic Acid Base-pair Geometry</em>
 * J. Mol. Biol. (2001) <bf>313</bf>, 229-237 <br/>
 * is defined as:
\begin{align*}
 y &= C_{1B}^{'} - C_{1B}^{'} \\
 y &= C_{1B}^{'} - C_{1B}^{'} \\
 z &= x \times \y \\
\end{align*}

 */
void local_coordinates_three_atoms(const Vec3 &a1, const Vec3 &a2, const Vec3 &a3, Rototranslation & r);

/** @brief Calculates Euler angles corresponding to a given rotation matrix.
 *
 * The rotation matrix should be provided row-wise.
 * @param row_x - the first row of the rotation matrix
 * @param row_y - the second row of the rotation matrix
 * @param row_z - the third row of the rotation matrix
 * @return three Euler angles stored in a 3D vector object
 */
Vec3 euler_angles(const Vec3 & row_x,const Vec3 & row_y,const Vec3 & row_z);

/** @brief Calculates Euler angles corresponding to a given rotation matrix.
 *
 * The rotation matrix is obtained from the provided instance of Rototranslation
 * @param r - rototranslation transformation
 * @return three Euler angles stored in a 3D vector object
 */
Vec3 euler_angles(const Rototranslation & r);

/** @brief Calculates Euler angles between two local reference frames
 *
 * @param r1 - the first rototranslation transformation
 * @param r2 - the second rototranslation transformation
 * @return three Euler angles stored in a 3D vector object
 */
Vec3 euler_angles(const Rototranslation & r1,const Rototranslation & r2);

}
}
}
}

#endif

/**
 * \example ex_local_coordinates_three_atoms.cc
 */
