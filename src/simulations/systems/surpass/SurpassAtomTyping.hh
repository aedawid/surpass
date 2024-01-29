#ifndef SIMULATIONS_CARTESIAN_SurpassAtomTyping_HH
#define SIMULATIONS_CARTESIAN_SurpassAtomTyping_HH

#include <vector>
#include <string>

#include <core/index.hh>

#include <simulations/systems/AtomTypingBase.hh>

namespace simulations {
namespace systems {
namespace surpass {

/** @brief Atom typing for <code>SURPASS</code> model
 */
class SurpassAtomTyping : public AtomTypingBase {
public:
  /// Names of all distinct atoms defined in the CABS++ model
  static const std::vector<std::string> surpass_atom_names;

  /// Constructor will initialize the object based on its static data
  SurpassAtomTyping() : AtomTypingBase(surpass_atom_names) {}

  /// Virtual destructor to satisfy a compiler
  virtual ~SurpassAtomTyping() {}
};

}
}
}

#endif

