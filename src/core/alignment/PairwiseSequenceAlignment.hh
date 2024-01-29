#ifndef CORE_ALIGNMENT_PairwiseSequenceAlignment_H
#define CORE_ALIGNMENT_PairwiseSequenceAlignment_H

#include <string>

#include <core/real.hh>
#include <core/alignment/PairwiseSequenceAlignment.fwd.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/io/fasta_io.hh>
#include <core/alignment/PairwiseAlignment.hh>

namespace core {
namespace alignment {

using namespace core::data::sequence;

/** @brief Represents a pairwise sequence alignment
 *
 * Pairwise sequence alignment may be created from an abstract alignment object and the two aligned sequences
 * (sequence profiles, secondary profiles, etc) as shown below:
 *
 * \include ex_PairwiseSequenceAlignment.cc
 */
class PairwiseSequenceAlignment {
public:

  PairwiseAlignment_SP alignment; ///< Alignment between the two sequences
  Sequence_SP query_sequence;     ///< The query sequence
  Sequence_SP template_sequence;  ///< The template sequence

  /** @brief Create a PairwiseSequenceAlignment object.
   *
   * This object will hold the given pointers as a shallow copies, i.e. any change to the alignment
   * pointed by <code>alignment</code> pointer will affect this object.
   * @param alignment - the alignment
   * @param query_sequence - query sequence data
   * @param template_sequence - template sequence data
   */
  PairwiseSequenceAlignment(const PairwiseAlignment_SP alignment, const Sequence_SP query_sequence, const Sequence_SP template_sequence) :
      alignment(alignment), query_sequence(query_sequence), template_sequence(template_sequence) {

    if (alignment->query_length() != query_sequence->length()) {
      throw std::runtime_error(
        utils::string_format("Wrong length of the query sequence %s! %d in alignment, %d in sequence!",
          query_sequence->header().c_str(), alignment->query_length(), query_sequence->length()));
    }
    if (alignment->template_length() != template_sequence->length()) {
      throw std::runtime_error(
        utils::string_format("Wrong length of the template sequence %s! %d in alignment, %d in sequence!",
          template_sequence->header().c_str(), alignment->template_length(), template_sequence->length()));
    }
  }

  /** @brief Create a PairwiseSequenceAlignment object along with its all data.
   *
   * This constructor creates a new PairwiseAlignment object based on the two aligned sequences. Respective Sequence
   * objects are created as well. Pointers to the newly created objects are stored within this PairwiseSequenceAlignment instance.
   *
   * @param alignment - the alignment
   * @param query_name - name of the query sequence; will be used to display the alignment and will be stored in a sequence header when printed in FASTA format
   * @param aligned_query - query sequence - already aligned with the template
   * @param template_name - name of the template sequence, will work as the name for the query
   * @param aligned_template - template sequence - as it is aligned with the query
   * @param score - alignment score, passed to a constructor of the  PairwiseAlignment object contained by this object
   * @param expand_unaligned - parameter passed  to a constructor of the PairwiseAlignment object contained by this object
   */
  PairwiseSequenceAlignment(const std::string & query_name, const std::string & aligned_query, const core::index2 first_query_residue,
      const std::string & template_name, const std::string &  aligned_template, const core::index2 first_template_residue, const core::real score,
      bool expand_unaligned = true);

  /** @brief Create a PairwiseSequenceAlignment object along with its all data.
   *
   * This constructor creates a new PairwiseAlignment object based on the two aligned sequences. Respective Sequence
   * objects are created from the given gapped sequence objects by removing all gaps. Pointers to the newly created objects
   * are stored within this PairwiseSequenceAlignment instance. Original data is not affected.
   *
   * @param aligned_query - query sequence - already aligned with the template
   * @param aligned_template - template sequence - as it is aligned with the query
   * @param score - alignment score
   */
  PairwiseSequenceAlignment(const Sequence_SP aligned_query, const Sequence_SP aligned_template, const core::real score);

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
};

std::ostream & operator<<(std::ostream & out, const PairwiseSequenceAlignment & seq_ali );

} // ~alignment
} // ~core

#endif
/**
 * \example ex_PairwiseSequenceAlignment.cc
 */
