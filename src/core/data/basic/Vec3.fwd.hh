#ifndef CORE_DATA_BASIC_VEC3_FWD_HH
#define CORE_DATA_BASIC_VEC3_FWD_HH

namespace core {
namespace data {
namespace basic {

#include <memory>

class Vec3;

/** @brief Type that represents a vector of points in 3D.
 *
 * Such coordinates may be conveniently loaded from a PDB file by using
 * utils::options::models_from_cmdline() declared in input_utils.hh
 */
typedef std::vector<core::data::basic::Vec3> Coordinates;

/** @brief Type that represents a pointer to a vector of points in 3D
 *
 * Such coordinates may be conveniently loaded at the command line level from a PDB file by using
 * utils::options::models_from_cmdline() declared in input_utils.hh
 */
typedef std::shared_ptr<std::vector<core::data::basic::Vec3>> Coordinates_SP;

} // ~ basic
} // ~ data
} // ~ core

#endif
