/** @file PairwiseStructureAlignment.hh
 *  @brief defines PairwiseStructureAlignment and PairwiseStructureAlignment_SP types
 */
#ifndef CORE_ALIGNMENT_PairwiseStructureAlignment_H
#define CORE_ALIGNMENT_PairwiseStructureAlignment_H

#include <string>

#include <core/alignment/PairwiseSequenceAlignment.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/alignment/CoordinatesSet.hh>

#include <core/protocols/StructureAlignmentProtocol.fwd.hh>

namespace core {
namespace alignment {

using namespace core::data::sequence;

/** @brief Represents a pairwise structure alignment.
 *
 * Pairwise structure alignment may be created from an abstract alignment object and the two structures
 */
class PairwiseStructureAlignment : public PairwiseSequenceAlignment {
  friend core::protocols::StructureAlignmentProtocol; // StructureAlignmentProtocol must be able to access job_id
public:

  const CoordinatesSet & query_structure; ///< Holds information about the query structure
  const CoordinatesSet & template_structure; ///< Holds information about the template structure

  /** @brief Create a PairwiseStructureAlignment object.
   *
   * This object will hold the given pointers as a shallow copies, i.e. any change to the alignment
   * pointed by <code>alignment</code> pointer will affect this object.
   * @param alignment - the alignment
   * @param query_sequence - query sequence data
   * @param template_sequence - template sequence data
   */
  PairwiseStructureAlignment(const PairwiseAlignment_SP alignment,
                            const CoordinatesSet & query_structure, const CoordinatesSet & template_structure) :
    PairwiseSequenceAlignment(alignment,query_structure.create_superimposed_sequence(),
      template_structure.create_superimposed_sequence()),query_structure(query_structure),template_structure(template_structure) {

  }

  /// Returns the total alignment score value stored in the alignment
  inline real alignment_score() const { return alignment->alignment_score; }

  /** @brief Returns the query sequence as it is aligned to the target.
   *
   * @param gap - character used to denote a gap in an alignment; by default '-' is used
   */
  inline std::string get_aligned_query(char gap = '-') const {

    return alignment->get_aligned_query(query_sequence->sequence,gap);
  }

  /** @brief Returns the template sequence as it is aligned to the target.
   *
   * @param gap - character used to denote a gap in an alignment; by default '-' is used
   */
  inline std::string get_aligned_template(const char gap='-') const {

    return alignment->get_aligned_template(template_sequence->sequence,gap);
  }

private:
  index4 job_id; ///< Used by StructureAlignmentProtocol to store index of an alignment that is computed in a concurrent thread
};

/// Define PairwiseStructureAlignment_SP type which is a shared pointer to PairwiseStructureAlignment
typedef std::shared_ptr<PairwiseStructureAlignment> PairwiseStructureAlignment_SP;

} // ~alignment
} // ~core

#endif
