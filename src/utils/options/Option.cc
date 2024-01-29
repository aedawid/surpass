#include <stdexcept>

#include <utils/options/Option.hh>
#include <utils/options/OptionParser.hh>

namespace utils {
namespace options {

bool Option::was_used() const {
  return OptionParser::get().was_used(*this);
}

}
}
