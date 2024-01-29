/** @file AtomTypingVariants.hh
 * @brief Provides AtomTypingVariants enumeration and ostream operators to print it nicely
 */
#ifndef SIMULATIONS_SYSTEMS_AtomTypingVariants_HH
#define SIMULATIONS_SYSTEMS_AtomTypingVariants_HH

#include <iostream>
#include <utils/Logger.hh>

namespace simulations {
namespace systems {

/** @brief AtomTypingVariant tells AtomTyping classes how to modify a standard residue type.
 *
 * For example, terminal residues in Molecular Mechanics force fields should be modified to include amine or carboxyl group.
 * This enum is only a flag denoting that a modification is necessary.
 */
enum class AtomTypingVariants  {
  STANDARD, // just as without any variant
  N_TERMINAL, // make H2N group
  C_TERMINAL  // make COO- group, add OXT atom
};

/// Operator to print AtomTypingVariants name to a std stream
std::ostream & operator<<(std::ostream & out,const AtomTypingVariants v);

/// Operator to print AtomTypingVariants name to logs
utils::Logger & operator<<(utils::Logger & logs,const AtomTypingVariants v);

}
}
#endif
