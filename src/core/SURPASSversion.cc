#include <sstream>

#include <core/SURPASSversion.hh>

namespace core {

const std::string SURPASSversion::GIT_HASH = "d81726d2747d5446bbe3b48f464f5b199f0e3ed1";
const std::string SURPASSversion::GIT_TIMESTAMP = "2024-01-29 06:10:25";
const std::string SURPASSversion::GIT_BRANCH = "surpass";

std::string SURPASSversion::to_string() {

  std::stringstream str;
  str << *this;
  return str.str();
}

std::ostream & operator<<(std::ostream & out, const SURPASSversion & ver) {

  out << " Git SHA       : " << ver.GIT_HASH << "\n"
      << " Git branch    : " << ver.GIT_BRANCH << "\n"
      << " Git timestamp : " << ver.GIT_TIMESTAMP << "\n";
  return out;
}

}
