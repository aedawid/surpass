/** \file structure_selectors.fwd.hh
 * @brief Provides forward declarations for types declared in structure_selectors.hh
 */

#ifndef CORE_DATA_STRUCTURAL_structure_selectors_FWD_HH
#define CORE_DATA_STRUCTURAL_structure_selectors_FWD_HH

#include <memory>

namespace core {
namespace data {
namespace structural {

class AtomSelector;

typedef std::shared_ptr<core::data::structural::AtomSelector> AtomSelector_SP;

class IsCA;
class IsCB;
class IsBB;
class IsBBCB;
class IsNamedAtom;

class ResidueSelector;
typedef std::shared_ptr<ResidueSelector> ResidueSelector_SP;
class SelectResidueRange;

class ChainSelector;
typedef std::shared_ptr<ChainSelector> ChainSelector_SP;

class SelectChainResidues;
class SelectChainResidueAtom;

class LogicalANDSelector;
}
}
}

#endif
