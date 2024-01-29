#ifndef CORE_ALIGNMENT_AlignmentRow_H
#define CORE_ALIGNMENT_AlignmentRow_H

#include <vector>
#include <string>
#include <iostream>

#include <core/index.hh>

namespace core {
namespace alignment {

class PairwiseAlignment;

/** @brief Represents a single sequence aligned with others i.e one row in a pairwise or multiple alignment.
 *
 * Note, that some of these aligned positions may actually be aligned with gaps from another alignment row
 */
class AlignmentRow {
public:

  friend PairwiseAlignment;

  /// Index of the first position in the aligned part of a sequence
  core::index2 first_pos;

  /** @brief Creates an empty object.
   *
   * The data will be filled by relevant methods from the friend classes, e.g. PairwiseAlignment class
   * @param first_pos - the index of the first residue position covered by this alignment row
   */
  AlignmentRow(const core::index2 first_pos) : first_pos(first_pos) { }

  /** @brief Creates an alignment row based on an aligned (i.e. gapped) sequence.
   *
   * The data will be filled by relevant methods from the friend classes, e.g. PairwiseAlignment class
   * @param first_pos - the index of the first residue position covered by this alignment row
   */
  AlignmentRow(const core::index2 first_pos, const std::string &aligned_seq) : first_pos(first_pos) {
    for (const char c : aligned_seq)
      if ((c == '_') || (c == '-')) append_gapped();
      else append_aligned();
  }

  /** @brief Copy constructor.
   * @param src - the source object
   */
  AlignmentRow(const AlignmentRow & src) : first_pos(src.first_pos) {

    for(auto e : src.alignment_to_sequence) alignment_to_sequence.push_back(e);
    for(auto e : src.sequence_to_alignment) sequence_to_alignment.push_back(e);
    last_pos_aligned = src.last_pos_aligned;
  }

  /// Returns the length of this Alignment row (including gaps)
  core::index2 length() const { return alignment_to_sequence.size(); }

  /** @brief Returns the number of aligned positions (which is the full length of this row minus the number of gaps).
   *
   * Note, that some of these aligned positions may actually be aligned with gaps from another alignment row
   */
  core::index2 n_aligned() const { return last_pos_aligned + 1; }

  /// Returns the number of gaps in this alignment row
  core::index2 n_gaps() const { return alignment_to_sequence.size() - last_pos_aligned; }

  /** @brief Says which residue of the relevant sequence appears at <code>alignment_pos</code> position of this AlignmentRow
 *
 * If the requested position in the alignment row is a gap, a negative number is returned.
 * @param alignment_pos - position within this alignment row (zero-related)
 */
  int position_in_sequence(const core::index2 alignment_pos) const { return alignment_to_sequence[alignment_pos]; }

  /** @brief Returns a sequence of arbitrary elements as it appears in an alignment aligned according to this AlignmentRow.
   *
   * This method takes a vector of elements of the aligned sequence (e.g. chars for a sequence string, Residue_SP, PdbAtom_SP or other objects)
   * and the gap symbol (e.g. '-', nullptr, etc) and produces the aligned sequence, stored in the given vector
   *
   * <strong>Note: </strong>This method clears the given <code>aligned_elements</code> vector prior populating it with new content
   *
   * @param sequence_elements - the elements of the sequence that will be aligned, e.g. characters, Residue objects or \f$C\alpha\f$ atoms
   * @param gap_element - the element used to denote a gap in the alignment
   * @param aligned_elements - the sequence aligned according to this alignment
   * @returns the reference to the given aligned_elements vector
   */
  template<typename T>
  std::vector<T> &aligned_sequence(const std::vector<T> &sequence_elements, const T gap_element,
                                   std::vector<T> &aligned_elements) const {

    aligned_elements.clear();
    for (const short int t:alignment_to_sequence) {
      if (t >= 0) aligned_elements.push_back(sequence_elements[t]);
      else aligned_elements.push_back(gap_element);
    }
    return aligned_elements;
  }

  /** @brief Returns a sequence as it appears in an alignment aligned according to this AlignmentRow.
   *
   * @param a_sequence - the elements of the sequence that will be aligned, e.g. characters, Residue objects or \f$C\alpha\f$ atoms
   * @param gap - the element used to denote a gap in the alignment
   * @returns the aligned sequence as a string
   */
  std::string aligned_sequence(const std::string &a_sequence, const char gap = '-') const {

    std::vector<char> t_chars(a_sequence.begin(), a_sequence.end());
    std::vector<char> t_gapped_chars;
    aligned_sequence(t_chars, gap, t_gapped_chars);
    return std::string(t_gapped_chars.begin(), t_gapped_chars.end());
  }


  /** @brief Returns a sequence of arbitrary elements as it appears in an alignment aligned according to this AlignmentRow.
   *
   * This method takes a vector of elements of the aligned sequence (e.g. chars for a sequence string, Residue_SP, PdbAtom_SP or other objects)
   * and produces the aligned sequence, stored in the given vector. Gaps are omitted by this method.
   *
   * <strong>Note: </strong>This method clears the given <code>aligned_elements</code> vector prior populating it with new content
   *
   * @param sequence_elements - the elements of the sequence that will be aligned, e.g. characters, Residue objects or \f$C\alpha\f$ atoms
   * @param aligned_elements - the sequence aligned according to this alignment
   * @returns the reference to the given aligned_elements vector
   */
  template<typename T>
  std::vector<T> &aligned_sequence(const std::vector<T> &sequence_elements, std::vector<T> &aligned_elements) const {

    aligned_elements.clear();
    for (const short int t:alignment_to_sequence) {
      std::cerr << t<<" "<<sequence_elements[t]<<" "<<aligned_elements.size()<<"\n";

      if (t >= 0) aligned_elements.push_back(sequence_elements[t]);
    }
    std::cerr << alignment_to_sequence.size()<<" "<<sequence_elements.size()<<" -> "<<aligned_elements.size()<<"\n";

    return aligned_elements;
  }

private:
  short last_pos_aligned = -1;
  /// Binds position in a sequence to an appropriate in this alignment row.
  std::vector<short int> sequence_to_alignment;
  /** @brief Binds alignment indexes to sequence indexes.
   *
   * <code>alignment_to_sequence[i]</code> holds an index in the sequence (or structure) that is related to this alignment row.
   * If position <code>i</code> of this alignment row is aligned to a gap, <code>alignment_to_sequence[i]</code> holds <code>-1</code>
   */
  std::vector<short int> alignment_to_sequence;

  inline void append_aligned() {
    alignment_to_sequence.push_back((++last_pos_aligned));
    sequence_to_alignment.push_back(alignment_to_sequence.size() - 1);
  }

  inline void append_gapped() { alignment_to_sequence.push_back(-1); }
};


} // ~alignment
} // ~core

#endif
