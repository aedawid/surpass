#ifndef CORE_ALIGNMENT_PairwiseAlignment_H
#define CORE_ALIGNMENT_PairwiseAlignment_H

#include <vector>
#include <string>
#include <algorithm>
#include <tuple>
#include <exception>
#include <memory>   // for shared_ptr

#include <core/index.hh>
#include <core/real.hh>
#include <core/alignment/AlignmentRow.hh>
#include <core/alignment/AlignmentBlock.hh>
#include <core/data/sequence/Sequence.hh>

#include <utils/Logger.hh>

namespace core {
namespace alignment {

typedef std::shared_ptr<PairwiseAlignment> PairwiseAlignment_SP;

/** @brief Represents a generic pairwise alignment.
 *
 * <strong>Note</strong> that this class doesn't know anything about objects that have been aligned. It just
 * holds a path on an alignment matrix and allow simple operations on it. Once you have a vector of
 * objects that correspond to the aligned residues, you may reorder them as they would appear in the alignment, e.g.:
 * \include ex_PairwiseAlignment.cc
 *
 * In the example above, the query "objects" are just characters denoting the residues, but in general
 * any object may be used.
 *
 * If you need a concrete sequence alignment class, use PairwiseSequenceAlignment.
 */
class PairwiseAlignment {
public:

  /// Total alignment score
  const real alignment_score;

  /** @brief Create a PairwiseAlignment object.
   *
   * @param first_query_residue - index of the first residue in the query sequence (important e.g. in the case of a local alignment)
   * @param first_template_residue - index of the first residue in the template sequence
   * @param score - total alignment score
   */
  PairwiseAlignment(const core::index2 first_query_residue, const core::index2 first_template_residue, const core::real score) :
    alignment_score(score), the_query(first_query_residue), the_template(first_template_residue), logs("PairwiseAlignment") { }

  /** @brief Create a PairwiseAlignment object from two strings representing an alignment
   *
   * @param aligned_query - query sequence aligned to a template
   * @param first_query_residue - index of the first residue in the query sequence (important e.g. in the case of a local alignment)
   * @param aligned_template - template sequence aligned to a query
   * @param first_template_residue - index of the first residue in the template sequence
   * @param score - total alignment score
   * @param expand_unaligned - whether unaligned pairs (which are marked by lower case letters) should be expanded, i.e.
   * <code><pre>
   *  LWagc
   *  IWdef
   * </pre></code>
   *
   * will be transformed to:
   * <code><pre>
   *  LW---AGC
   *  IWDEF---
   *  </pre></code>
   */
  PairwiseAlignment(const std::string &aligned_query, const core::index2 first_query_residue, const std::string &aligned_template,
                    const core::index2 first_template_residue, const core::real score = 0.0, bool expand_unaligned = true);

  /** @brief Create a PairwiseAlignment instance based on two aligned sequences.
   *
   * The given sequences are used only to extract an alignment path from them.
   *
   * @param aligned_query - a query sequence aligned to its template
   * @param aligned_template - a template sequence aligned to the query
   * @param score - alignment score value (mainly used when the alignment is printed on a screen or to a file)
   */
  PairwiseAlignment(const core::data::sequence::Sequence &aligned_query,
                    const core::data::sequence::Sequence &aligned_template, const core::real score);

  /** @brief Create a PairwiseAlignment instance based on an alignment path
   *
   *
   * @param first_query_residue - a query sequence aligned to its template
   * @param first_template_residue - a template sequence aligned to the query
   * @param score - alignment score value (mainly used when the alignment is printed on a screen or to a file)
   * @param alignment_path - alignment path is a string that encodes an alignment. It may contain only three different characters:
   *    - '*' means a match (two positions are aligned)
   *    - '-' means gap in the first (query) sequence
   *    - '|' means gap in the second (template) sequence
   */
  PairwiseAlignment(const core::index2 first_query_residue, const core::index2 first_template_residue,
                    const core::real score, const std::string &alignment_path);

  /** @brief Creates a new alignment from an alignment path.
   *
   * @see PairwiseAlignment(const core::index2,const core::index2,  const core::real, const std::string &) docs
   *    for description of the alignment path format
   */
  PairwiseAlignment(const std::string &alignment_path) : PairwiseAlignment(0, 0, 0.0, alignment_path) { }

  /** @brief Creates a pairwise alignment from two AlignmentRow objects.
   *  @param score - score value of this alignment
   *  @param query - the query alignment row
   *  @param tmplt - the template alignment row
   */
  PairwiseAlignment(const real score, const AlignmentRow &query, const AlignmentRow &tmplt);

  /// Says which position in the query is aligned with the given position in the template
  inline int which_template_for_query(const core::index2 query_index) const {
    return the_template.alignment_to_sequence[the_query.sequence_to_alignment[query_index]];
  }

  /// Says which position in the template is aligned with the given position in the query
  inline int which_query_for_template(const core::index2 template_index) const {
    return the_query.alignment_to_sequence[the_template.sequence_to_alignment[template_index]];
  }

  /// Returns the total length of this alignment
  inline core::index2 length() const { return the_query.alignment_to_sequence.size(); }

  /// Returns the number of aligned positions from the query sequence (gaps are excluded)
  inline core::index2 query_length() const { return the_query.n_aligned(); }

  /// Returns the number of residue-to-residue alignment positions
  inline core::index2 n_aligned() const { return n_aligned_; }

  /// Returns the number of aligned positions from the template sequence (gaps are excluded)
  inline core::index2 template_length() const { return the_template.n_aligned(); }

  /// Appends a match to this alignment path
  inline void append_match() {
    the_query.append_aligned();
    the_template.append_aligned();
    n_aligned_++;
  }

  inline void append_query_gap() {
    the_query.append_gapped();
    the_template.append_aligned();
  }

  inline void append_template_gap() {
    the_template.append_gapped();
    the_query.append_aligned();
  }

  inline const std::pair<short int, short int> at(const core::index2 alignment_position) const {

    return std::pair<short int, short int>(the_query.alignment_to_sequence[alignment_position],
      the_template.alignment_to_sequence[alignment_position]);
  }

  /// Returns true if this sequence contains at least one gap
  inline bool is_gapped(core::index2 position) const {

    return ((the_query.alignment_to_sequence[position] < 0) || (the_template.alignment_to_sequence[position] < 0));
  }

  /** @brief Returns a vector of objects as they would be a query sequence aligned to the template.
   *
   * This is a generic method that takes a vector of objects of type <code>T</code>
   * and distributes them as they would be the query sequence in this alignment.
   * The given <code>gap</code> object is used to represent gaps in the alignment
   * @param query_residues - a sequence of query object that will be aligned
   * @param gap - object used to represent a gap.
   * @param aligned_query_data - vector where the result is stored
   * @tparam T - type of aligned objects, e.g. Residue, PdbAtom, char etc.
   */
  template<typename T>
  inline std::vector<T> &get_aligned_query(const std::vector<T> &query_residues, const T gap,
                                           std::vector<T> &aligned_query_data) const {

    return the_query.aligned_sequence(query_residues, gap, aligned_query_data);
  }

  /** @brief Formats an aligned template sequence given a raw (gapless) sequence.
   *
   * @param template_sequence - the sequence aligned as the template
   * @param gap - the gap symbol
   */
  inline std::string get_aligned_query(const std::string &query_sequence, const char gap = '-') const {

    return the_query.aligned_sequence(query_sequence, gap);
  }

  /** @brief Formats an aligned template sequence given a raw (gapless) sequence.
   *
   * @param template_sequence - the sequence aligned as the template
   * @param gap - the gap symbol
   */
  inline std::string get_aligned_template(const std::string &template_sequence, const char gap = '-') const {

    return the_template.aligned_sequence(template_sequence, gap);
  }

  /** @brief Returns a vector of objects as they would be a template sequence aligned to the query.
   *
   * @param template_residues - a sequence of template object that will be aligned
   * @param gap - object used to represent a gap.
   * @param aligned_template_data - vector where the result is stored
   * @see get_aligned_query<T>() - a similar method that retrieves an aligned query
   */
  template<typename T>
  inline std::vector<T> &get_aligned_template(std::vector<T> &template_residues, T gap,
                                              std::vector<T> &aligned_template_data) const {

    return the_template.aligned_sequence(template_residues, gap, aligned_template_data);
  }

  /** @brief Creates a vector of query objects aand a vector of template objects as they are aligned with themselves.
   *
   * @param query_residues - a sequence of query object that will be aligned
   * @param tmplt_residues - a sequence of template object that will be aligned
   * @param gap - object used to represent a gap.
   * @param aligned_query_data - vector where the resulting query is stored
   * @param aligned_tmplt_data - vector where the resulting template is stored
   * @see get_aligned_query<T>(std::vector<T> &,std::vector<T> &) - a similar method that retrieves an aligned query
   */
  template<typename T>
  void get_aligned_query_template(const std::vector<T> &query_residues, const std::vector<T> &tmplt_residues,
                                  const T &gap, std::vector<T> &aligned_query_data,
                                  std::vector<T> &aligned_tmplt_data) const {

    the_query.aligned_sequence(query_residues, gap, aligned_query_data);
    the_template.aligned_sequence(tmplt_residues, gap, aligned_tmplt_data);
    typename std::vector<T>::iterator q_it = aligned_query_data.begin();
    typename std::vector<T>::iterator t_it = aligned_tmplt_data.begin();
    while (q_it != aligned_query_data.end()) {
      if ((*q_it == gap) || (*t_it == gap)) {
        t_it = aligned_tmplt_data.erase(t_it);
        q_it = aligned_query_data.erase(q_it);
      } else {
        ++q_it;
        ++t_it;
      }
    }
  }

  /** @brief Creates a vector of query objects and a vector of template objects as they are aligned with themselves.
   *
   * This method does not put gaps into the aligned sequence.
   * @param query_residues - a sequence of query object that will be aligned
   * @param tmplt_residues - a sequence of template object that will be aligned
   * @param aligned_query_data - vector where the resulting query is stored
   * @param aligned_tmplt_data - vector where the resulting template is stored
   */
  template<typename T>
  void get_aligned_query_template(const std::vector<T> &query_residues, const std::vector<T> &tmplt_residues,
                                  std::vector<T> &aligned_query_data, std::vector<T> &aligned_tmplt_data) const {

    aligned_query_data.clear();
    aligned_tmplt_data.clear();

    core::index2 iq = 0, it = 0;
    for (core::index2 i = 0; i < the_query.alignment_to_sequence.size(); ++i) {
      if (the_query.position_in_sequence(i) < 0) ++it;
      if (the_template.position_in_sequence(i) < 0) ++iq;
      if ((the_query.position_in_sequence(i) >= 0) && (the_template.position_in_sequence(i) >= 0)) {
        aligned_query_data.push_back(query_residues[iq]);
        aligned_tmplt_data.push_back(tmplt_residues[it]);
        ++it;
        ++iq;
      }
    }
  }

  /** \brief Creates a path string this encodes the alignment.
   */
  std::string to_path() const;

  /// Returns a vector of all solid (i.e. gapless) blocks of this alignment
  /**
   * A solid block is the longest possible fragment extracted from an alignment without any gap; neither in the query
   * nor in the template sequence.
   *
   * @returns a vector of aligned blocks; each block is returned as a four-tuple, providing:
   * 		- index of the first residue of the block in the first (i.e. query) sequence
   * 		- index of the last residue of the block in the first (i.e. query) sequence
   * 		- index of the first residue of the block in the second (i.e. template) sequence
   * 		- index of the last residue of the block in the second (i.e. template) sequence
   */
  std::vector<AlignmentBlock> aligned_blocks() const;

private:
  AlignmentRow the_query;
  AlignmentRow the_template;
  core::index2 n_aligned_ = 0;
  utils::Logger logs;
  void clear_unaligned(core::index2 &n, const bool flag);

  /** @brief increments the first residue index in the_row.
   * The method is used when a new alignment is createdfrom two alignment rows. The value of <code>first_index</code>
   * must be increased if both rows start with a gap. Such gaps aligned to other gaps will be removed by
   * PairwiseAlignment constructor.
   */
  static int first_index(int index, const AlignmentRow &the_row, const AlignmentRow &the_other_row);
};

} // ~alignment
} // ~core

/**
 * \example ex_PairwiseAlignment.hh
 */
#endif
