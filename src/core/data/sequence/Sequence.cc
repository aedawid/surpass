#include <string>
#include <algorithm>

#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SequenceFilter.hh>

namespace core {
namespace data {
namespace sequence {

const char Sequence::gap_symbol = '-';

const core::chemical::Monomer & Sequence::get_monomer(core::index4 pos) const {
  if (sequence[pos] != gap_symbol || pos == 0) return core::chemical::Monomer::get(sequence[pos]);
  return (sequence[pos - 1] == gap_symbol) ? core::chemical::Monomer::GPE : core::chemical::Monomer::GAP;
}

std::string Sequence::fix_gaps(const char* & input_seq) {

  std::string s(input_seq);
  std::replace(s.begin(), s.end(), '_', gap_symbol);
  if (gap_symbol != '-') std::replace(s.begin(), s.end(), '-', gap_symbol);
  return s;
}

std::string Sequence::fix_gaps(const std::string & input_seq) {

  std::string s(input_seq);
  std::replace(s.begin(), s.end(), '_', gap_symbol);
  if (gap_symbol != '-') std::replace(s.begin(), s.end(), '-', gap_symbol);
  return s;
}

std::string Sequence::create_seq(const std::vector<core::chemical::Monomer> seq) {

    std::string out;
    for (const core::chemical::Monomer & m : seq) out += m.code1;
    return out;
}

bool Sequence::is_protein(const std::string &seq) {
  static const IsProteinSequence is_protein_sequence_;
  return is_protein_sequence_(seq, seq);
}

bool Sequence::is_nucleic(const std::string &seq) {
  static const IsNucleicSequence is_nucleic_sequence_;
  return is_nucleic_sequence_(seq, seq);
}


}
}
}
