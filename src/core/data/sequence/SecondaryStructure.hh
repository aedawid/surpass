#ifndef CORE_DATA_SEQUENCE_SecondaryStructure_H
#define CORE_DATA_SEQUENCE_SecondaryStructure_H

#include <string>
#include <vector>
#include <memory>

#include <core/index.hh>
#include <core/real.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SecondaryStructureAnnotation.hh>
#include <core/data/sequence/SecondaryStructure.fwd.hh>

namespace core {
namespace data {
namespace sequence {

/// Represents a sequence annotated with secondary structure probabilities.
class SecondaryStructure : public Sequence, public SecondaryStructureAnnotation {
public:

  /** @brief Creates a secondary structure of a coiled chain.
   *
   * This secondary structure will have all secondary structure set to coil (C)
   * @param header - sequence header, passed to the relevant Sequence base class constructor
   * @param seq - amino acid or nucleic sequence (required by Sequence base class constructor)
   * @param first_pos - index of the very first residue in the sequence (again, required by Sequence base class constructor)
   */
  SecondaryStructure(const std::string & header, const std::string & seq, core::index2 first_pos) :
      Sequence(header, seq, first_pos), SecondaryStructureAnnotation(seq.size()) { }

  /** @brief Creates a secondary structure from a string-type definition.
   *
   * This secondary structure will have all probabilities equal to 1.0 or 0.0, according to the given string
   * @param header - sequence header, passed to the relevant Sequence base class constructor
   * @param seq - amino acid or nucleic sequence (required by Sequence base class constructor)
   * @param first_pos - index of the very first residue in the sequence (again, required by Sequence base class constructor)
   * @param ss_string - secondary structure as a string
   */
  SecondaryStructure(const std::string & header, const std::string & seq, core::index2 first_pos,
      const std::string & ss_string) :
      Sequence(header, seq, first_pos), SecondaryStructureAnnotation(ss_string) {}

  /** @brief Creates a new secondary structure  based on a contigus fragment of a source data
   * @param source - original  secondary structure object
   * @param start_pos - zero-related index indicates the first position of the new sequence
   * @param end_pos - zero-related index indicates the last position of the new sequence <strong>inclusive!</strong>
   */
  SecondaryStructure(const SecondaryStructure & source, const core::index2 start_pos, const core::index2 end_pos) :
    SecondaryStructure(source.header(),source.sequence.substr(start_pos-source.first_pos(), end_pos-source.first_pos()+1),
      start_pos, source.str().substr(start_pos-source.first_pos(), end_pos-source.first_pos()+1)) {}

  /// Bare virtual destructor to satisfy a compiler
  virtual ~SecondaryStructure() {}

  /** @brief Says whether this object bears information about secondary structure for this sequence.
   * @return always true
   */
  virtual bool has_ss() const { return true; }

  /// Creates and returns a new SecondaryStructure object by removing gaps from this sequence
  virtual std::shared_ptr<Sequence> create_ungapped_sequence() const {

    std::string tmp_seq = sequence;
    tmp_seq.erase(std::remove(tmp_seq.begin(), tmp_seq.end(), '_'), tmp_seq.end());
    std::shared_ptr<SecondaryStructure> s = std::make_shared<SecondaryStructure>(header(), tmp_seq, first_pos());
    core::index4 i = 0;
    for (core::index4 j = 0; j < sequence.size(); ++i) {
      if ((sequence[j] != '-') && (sequence[j] != '-')) {
        s->fractions(i, fraction_H(j), fraction_E(j), fraction_C(j));
        ++i;
      }
    }
    return s;
  }
};

}
}
}

#endif
