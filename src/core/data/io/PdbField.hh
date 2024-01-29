/** \file PdbField.hh
 * @brief A base class that represents a field found in a PDB file, e.g. SEQRES, HEADER or ATOM
 */
#ifndef CORE_DATA_IO_PdbField_H
#define CORE_DATA_IO_PdbField_H

#include <string>

namespace core {
namespace data {
namespace io {

/** @brief Base class for entries that may be found in a PDB file.
 */
class PdbField {
public:
  /// Object of any derived class may be printed in PDB format
  virtual std::string to_pdb_line() const = 0;
  /// virtual destructor to satisfy a compiler
  inline virtual ~PdbField() {}
};

}
}
}

#endif
