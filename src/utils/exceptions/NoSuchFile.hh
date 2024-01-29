#ifndef UTILS_EXCEPTIONS_NoSuchFile_H
#define UTILS_EXCEPTIONS_NoSuchFile_H

#include <stdexcept>
#include <string>

namespace utils {
namespace exceptions {

class NoSuchFile : public std::runtime_error {
public:
  const std::string missing_file_name;
  NoSuchFile(const std::string & missing_file_name) :
      std::runtime_error("Can't locate the file: " + missing_file_name), missing_file_name(missing_file_name) {
  }
};

}
}

#endif
