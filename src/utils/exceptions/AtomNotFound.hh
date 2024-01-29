#ifndef UTILS_EXCEPTIONS_AtomNotFound_H
#define UTILS_EXCEPTIONS_AtomNotFound_H

#include <stdexcept>
#include <string>

namespace utils {
namespace exceptions {

/** @brief Exception thrown when a residue lacks one of its atoms
 *
 */
class AtomNotFound : public std::runtime_error {
public:
  /// The name of the missing atom
  const std::string missing_atom_name;
  /// The name of the residue the  atom could not be found
  const std::string source_residue_name;

  /// Constructor creates an exception instance
  AtomNotFound(const std::string & missing_atom_name, const std::string & residue) :
      std::runtime_error("Can't locate the atom: " + missing_atom_name+" in residue "+residue),
      missing_atom_name(missing_atom_name), source_residue_name(residue) {
  }
};

}
}

#endif
