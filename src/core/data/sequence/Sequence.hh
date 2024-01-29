#ifndef CORE_DATA_SEQUENCE_Sequence_H
#define CORE_DATA_SEQUENCE_Sequence_H

#include <string>
#include <memory>
#include <vector>
#include <algorithm>

#include <core/data/sequence/Sequence.fwd.hh>

#include <core/data/sequence/sequence_utils.hh>
#include <core/chemical/Monomer.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace sequence {

/** @brief Represents amino acid or nucleic acid sequence.
 *
 * The following example reads sequences from a file in FASTA format and filters them by a few criteria.
 * \include ex_filter_fasta.cc
 */
class Sequence {
public:
  static const char gap_symbol; ///< The symbol used by BioShell to mark gaps in alignment; other commonly used gap symbols will be converted to this one
  const std::string sequence; ///< the sequence itself - as a string
  const bool is_nucleic_sequence; ///< True if this is a nucleotide sequence
  const bool is_protein_sequence; ///< True if this is a protein sequence

  /** @brief  Basic constructor fills the data fields of this instance.
   * @param header - a string that will show up as a header when this sequence will be printed in FASTA or PIR format
   * @param seq - the sequence itself (single-character code)
   * @param first_pos - if this sequence is actually a subsequence (e.g. comes from a local alignment), use a non-zero value
   */
  Sequence(const std::string & header, const std::string & seq, const core::index2 first_pos = 0) :
    sequence(fix_gaps(seq)), is_nucleic_sequence(is_nucleic(sequence)), is_protein_sequence((!is_nucleic_sequence) && (is_protein(sequence))), first_pos_(first_pos), header_(header) {}

  /** @brief  Basic constructor fills the data fields of this instance.
   * @param header - a string that will show up as a header when this sequence will be printed in FASTA or PIR format
   * @param seq - the sequence itself (single-character code)
   * @param first_pos - if this sequence is actually a subsequence (e.g. comes from a local alignment), use a non-zero value
   */
  Sequence(const char* header, const char* seq, const core::index2 first_pos = 0) :
    sequence(fix_gaps(seq)), is_nucleic_sequence(is_nucleic(sequence)), is_protein_sequence((!is_nucleic_sequence) && (is_protein(sequence))), first_pos_(first_pos), header_(header) {}

  /** @brief Creates a sequence from a vector of residue types
   * @param header - a string that will show up as a header when this sequence will be printed in FASTA or PIR format
   * @param seq - the sequence itself (vector of residue types)
   * @param first_pos - if this sequence is actually a subsequence (e.g. comes from a local alignment), use a non-zero value
   */
  Sequence(const char* header, const std::vector<core::chemical::Monomer> seq, const core::index2 first_pos = 0) :
    sequence(create_seq(seq)), is_nucleic_sequence(is_nucleic(sequence)), is_protein_sequence((!is_nucleic_sequence) && (is_protein(sequence))), first_pos_(first_pos), header_(header) {}

  /** @brief Creates a new sequence based on a contigus fragment of a source sequence
   * @param sequence - original sequence
   * @param start_pos - identifies the first residue of the copied fragment.
   * @param end_pos - identifies the last residue of the copied fragment. <strong>inclusive!</strong>
   */
  Sequence(const Sequence & seq, const core::index2 start_pos, const core::index2 end_pos) :
    sequence(seq.sequence.substr(start_pos-seq.first_pos(), end_pos-seq.first_pos()+1)), is_nucleic_sequence(seq.is_nucleic_sequence),
    is_protein_sequence(seq.is_protein_sequence), first_pos_(start_pos), header_(seq.header()) {}

  /// Bare virtual destructor to satisfy a compiler
  virtual ~Sequence() {}

  /// Returns the index of the first residue in this sequence
  core::index2 first_pos() const { return first_pos_; }

  /// Returns the header of this sequence
  const std::string &header() const { return header_; }

  /// Defines the new staring position for this sequence
  void first_pos(core::index2 new_index) { first_pos_ = new_index; }

  /** @brief  Defines a new header string for this sequence.
   *
   * This method is virtual so derived class can alter the header visible for this base class.
   * <code>const std::string &header() const</code> method however is not virtual, so overriding changes the way how
   * a header is formatted
   * @param new_header - the new header for this sequence
   */
  virtual void header(const std::string &new_header) { header_ = new_header; }

  /// Returns the number of residues in this sequence (including gaps!)
  size_t length() const { return sequence.length(); }

  /// true if and only if there is no single gap in this sequence
  inline bool is_gapped() const { return (sequence.find('-') != std::string::npos); }

  /// Returns the number of residues in this sequence that are not gaps
  inline core::index2 length_without_gaps() const { return std::count_if(sequence.begin(), sequence.end(), [](char c) {return (c!= '-')&&(c !='_');}); }

  /// Returns the number of residues in this sequence that are not gaps
  inline core::index2 count_gaps() const { return std::count_if(sequence.begin(), sequence.end(), [](char c) {return (c== '-')||(c =='_');}); }

  /** @brief Returns the residue type at a given position.
   * @param pos - position in this sequence
   * @return a monomer type defining the residue at the requested position
   */
  const core::chemical::Monomer & get_monomer(core::index4 pos) const;

  /** @brief Returns a character denoting the residue type at a given position.
   * @param pos - position in this sequence
   * @return a character defining the residue at the requested position
   */
  inline char operator[](const index4 pos) const {
    return sequence[pos];
  }

  /// Creates and returns a new sequence by removing gaps from this sequence
  virtual std::shared_ptr<Sequence> create_ungapped_sequence() const {

    std::string tmp_seq = sequence;
    tmp_seq.erase(std::remove(tmp_seq.begin(), tmp_seq.end(), '-'), tmp_seq.end());
    tmp_seq.erase(std::remove(tmp_seq.begin(), tmp_seq.end(), '_'), tmp_seq.end());
    return std::make_shared<Sequence>(header_, tmp_seq, first_pos());
  }

  /** @brief Says whether this object bears information about secondary structure for this sequence.
   * @return this class always returns false, but that may be changed in a derived type
   */
  virtual bool has_ss() const { return false; }

  /** @brief Returns true if the given string represents an amino acid sequence
   * IsProteinSequence filter  instance is used to make the decision
   * @param seq - a putative amino acid sequence
   * @return true if this sequence comprises at least 90% amino acid
   */
  static bool is_protein(const std::string &seq);

  /** @brief Returns true if the given string represents a nucleic acid sequence.
   * IsNucleicSequence filter  instance is used to make the decision
   * @param seq - a putative nucleic acid sequence
   * @return true if this is a nucleic acid sequence;
   */
  static bool is_nucleic(const std::string &seq);

private:
  core::index2 first_pos_;    ///< index of the first residue in the sequence; indexes for all the other residues will be assigned consecutively
  std::string header_; ///< a string identifying this sequence

  static std::string fix_gaps(const std::string & input_seq);
  static std::string fix_gaps(const char* & input_seq);

  static std::string create_seq(const std::vector<core::chemical::Monomer> seq);
};

}
}
}

/**
 * \example ex_filter_fasta.cc
 */

#endif
