#include <iostream>

#include <core/data/basic/Vec3.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace basic {

using core::real;

std::ostream& operator<<(std::ostream &out, const Vec3 &v) {

  out << utils::string_format("%8.3f %8.3f %8.3f", v.x, v.y, v.z);
	return out;
}

utils::Logger &operator <<(utils::Logger &logger, const Vec3 & coordinates) {

	logger << utils::string_format("%8.3f %8.3f %8.3f",coordinates.x,coordinates.y,coordinates.z);

	return logger;
}
} // ~ basic
} // ~ data
} // ~ core
