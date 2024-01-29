#include <iostream>

#include <core/data/basic/Vec3.hh>
#include <utils/string_utils.hh>
#include <core/calc/structural/transformations/Rototranslation.hh>

using core::data::basic::Vec3;

namespace core {
namespace calc {
namespace structural {
namespace transformations {

Rototranslation::Rototranslation() :
    rot_x_(1, 0, 0), rot_y_(0, 1, 0), rot_z_(0, 0, 1), tr_before_(0.0), tr_after_(0.0) { update_mm(); }

std::ostream& operator<<(std::ostream &out, const Rototranslation &r) {

	out << utils::string_format("x`  =  |%8.3f %8.3f %8.3f|     (x - %8.3f)     %8.3f\n",r.rot_x_.x,r.rot_x_.y,r.rot_x_.z,r.tr_before_.x,r.tr_after_.x);
	out << utils::string_format("y`  =  |%8.3f %8.3f %8.3f|  *  (y - %8.3f)  +  %8.3f\n",r.rot_y_.x,r.rot_y_.y,r.rot_y_.z,r.tr_before_.y,r.tr_after_.y);
	out << utils::string_format("z`  =  |%8.3f %8.3f %8.3f|     (z - %8.3f)     %8.3f\n",r.rot_z_.x,r.rot_z_.y,r.rot_z_.z,r.tr_before_.z,r.tr_after_.z);
	return out;
}

}
}
}
}
