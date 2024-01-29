/** @file AlignmentBlock.hh
 * Provides data type representing a single contiguous block of a sequence alignment
 */
#ifndef CORE_ALIGNMENT_AlignmentBlock_H
#define CORE_ALIGNMENT_AlignmentBlock_H

#include <memory>   // for shared_ptr
#include <iostream>

#include <core/index.hh>
#include <utils/string_utils.hh>

namespace core {
namespace alignment {

/** @brief Represents a contiguous (i.e. gapless) block of an alignment.
 *
 * Any alignment may be written as a combination of gapless aligned blocks linked by stretches of sequence fragments aligned with gaps
 */
struct AlignmentBlock {

  core::index2 first_query_pos; ///< Index of the first residue of a query sequence included in this block
  core::index2 last_query_pos; ///< Index of the last residue of a query sequence included in this block
  core::index2 first_tmplt_pos; ///< Index of the first residue of a template sequence included in this block
  core::index2 last_tmplt_pos; ///< Index of the last residue of a template sequence included in this block

  /** @brief Creates a new block.
   *
   * @param first_q - index of the first query residue
   * @param last_q - index of the last query residue
   * @param first_t - index of the first query template
   * @param last_t - index of the last query template
   */
  AlignmentBlock(const core::index2 first_q,const core::index2 last_q,const core::index2 first_t,const core::index2 last_t) {
    first_query_pos = first_q;
    first_tmplt_pos = first_t;
    last_query_pos = last_q;
    last_tmplt_pos = last_t;
  }

  /// Calculate the number of aligned residues in this block
  core::index2 size() const { return last_tmplt_pos - first_tmplt_pos + 1; }
};

typedef std::shared_ptr<AlignmentBlock> AlignmentBlock_SP;

std::ostream & operator<<(std::ostream & out, const AlignmentBlock & block);

std::ostream & operator<<(std::ostream & out, const AlignmentBlock_SP block);


}
}

#endif
