/** @file PairwiseSequenceAlignment.fwd.hh
 * Provides declaration of PairwiseSequenceAlignment_SP type
 */
#ifndef CORE_ALIGNMENT_PairwiseSequenceAlignment_FWD_H
#define CORE_ALIGNMENT_PairwiseSequenceAlignment_FWD_H

#include <memory>

namespace core {
namespace alignment {

class PairwiseSequenceAlignment;

/// PairwiseSequenceAlignment_SP is a shared pointer to PairwiseSequenceAlignment type
typedef std::shared_ptr<PairwiseSequenceAlignment> PairwiseSequenceAlignment_SP;

}
}

#endif
