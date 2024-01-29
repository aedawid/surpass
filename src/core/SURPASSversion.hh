#ifndef CORE_SURPASSversion_HH
#define CORE_SURPASSversion_HH

#include <string>

#include <utils/Logger.hh>

namespace core {

struct SURPASSversion {
  static const std::string GIT_HASH;
  static const std::string GIT_TIMESTAMP;
  static const std::string GIT_BRANCH;

  std::string to_string();
};

std::ostream & operator<<(std::ostream & ostream, const SURPASSversion & ver);

}

#endif
