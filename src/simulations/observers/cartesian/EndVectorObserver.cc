#include <core/data/basic/Vec3.hh>
#include <core/data/basic/Vec3Cubic.hh>
#include <simulations/observers/cartesian/EndVectorObserver.hh>

namespace simulations {
namespace observers {
namespace cartesian {

template<typename C>
bool EndVectorObserver<C>::observe() {

  ++cnt;
  if(!ObserverInterface::trigger->operator()()) return false;

  core::real cx = observed_object_[0].closest_delta_x(observed_object_[observed_object_.n_atoms-1]);
  core::real cy = observed_object_[0].closest_delta_y(observed_object_[observed_object_.n_atoms-1]);
  core::real cz = observed_object_[0].closest_delta_z(observed_object_[observed_object_.n_atoms-1]);

  *outstream << utils::string_format("%6d %8.3f %8.3f %8.3f\n", cnt, cx, cy, cz);
  outstream->flush();

  return true;
}

template class EndVectorObserver<core::data::basic::Vec3>;
template class EndVectorObserver<core::data::basic::Vec3Cubic>;

}
}
}
