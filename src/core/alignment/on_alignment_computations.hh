/** \file on_alignment_computations.hh
 * @brief Performs some basic calculation based on a sequence alignment (multiple or pairwise)
 *
 */
#ifndef CORE_ALIGNMENT_on_alignment_computations_HH
#define CORE_ALIGNMENT_on_alignment_computations_HH

#include <string>

#include <core/index.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/alignment/PairwiseAlignment.hh>
#include <core/alignment/PairwiseSequenceAlignment.hh>
#include <core/alignment/scoring/SimilarityMatrix.hh>

namespace core {
namespace alignment {

/// Compares every i-th position in the first aligned sequence with i-th position in the second one and counts how many of them are identical
/**
 * @param s1 - reference to the first of the two aligned sequences
 * @param s2 - reference to the second of the two aligned sequences
 */
core::index2 sum_identical(const std::string &s1, const std::string &s2);

/// Compares every i-th position in the first sequence with i-th position in the second one and counts how many of them are identical
/**
 * @param s1 - reference to the first of the two aligned sequences
 * @param s2 - reference to the second of the two aligned sequences
 */
core::index2 sum_identical(const core::data::sequence::Sequence &s1, const core::data::sequence::Sequence &s2);

/// Compares every i-th position in the first sequence with i-th position in the second one and counts how many of them are identical
/**
 * @param alignment - alignment saying how the two sequences are aligned to one another
 * @param query_sequence - reference to the first of the two compared sequences (query sequence)
 * @param tmplate_sequence - reference to the second of the two compared sequences (template sequence)
 */
core::index2 sum_identical(const PairwiseAlignment &alignment, const std::string &query_sequence,
                           const std::string &tmplate_sequence);

/// Compares every i-th position in the first sequence with i-th position in the second one and counts how many of them are identical
/**
 * @param alignment - pairwise sequence alignment
 */
core::index2 sum_identical(const PairwiseSequenceAlignment &alignment);

/// Calculates the number of aligned positions in a given sequence alignment
/**
 * @param s1 - reference to the first of the two aligned sequences
 * @param s2 - reference to the second of the two aligned sequences
 */
core::index2 sum_aligned(const std::string &s1, const std::string &s2);

/// Calculates the number of aligned positions in a given sequence alignment
/**
 * @param s1 - reference to the first of the two aligned sequences
 * @param s2 - reference to the second of the two aligned sequences
 */
core::index2 sum_aligned(const core::data::sequence::Sequence &s1, const core::data::sequence::Sequence &s2);

/// Compares every i-th position in the first sequence with i-th position in the second one and counts how many of them are identical
/**
 * @param alignment - alignment saying how the two sequences are aligned to one another
 */
core::index2 sum_aligned(const PairwiseAlignment &alignment);

/// Computes Aln_k measure between two alignments.
/**
 * @param reference_alignment - reference alignment
 * @param k - maximum alignment shift allowed
 * @param another_alignment - a score alignment
 */
core::index2 compute_aln_n(const PairwiseAlignment &reference_alignment, const core::index2 k,
                           const PairwiseAlignment &another_alignment);

/** @brief Calculates an alignment score for two sequences that are already aligned
 * @param s1 - reference to the first of the two aligned sequences
 * @param s2 - reference to the second of the two aligned sequences
 */
short int calculate_score(const core::data::sequence::Sequence &s1, const core::data::sequence::Sequence &s2,
                             const scoring::NcbiSimilarityMatrix &scoring, short int gap_open, short int gap_extend);

/** @brief Calculates an alignment score for two sequences that are already aligned
 * @param s1 - the first of the two aligned sequences
 * @param s2 - the second of the two aligned sequences
 */
short int calculate_score(const std::string &s1, const std::string &s2,
                             const scoring::NcbiSimilarityMatrix &scoring, short int gap_open, short int gap_extend);

}
}

#endif
