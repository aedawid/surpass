#include <string>
#include <cstdlib>
#include <regex>

#include <core/index.hh>
#include <core/data/sequence/sequence_utils.hh>
#include <unordered_set>

namespace core {
namespace data {
namespace sequence {

core::index4 extract_gid(const std::string & sequence_header) {

  static const std::regex find_gi_regex("gi\\|([0-9]+)\\|");

  std::smatch base_match;
  if (std::regex_search(sequence_header, base_match, find_gi_regex)) {
    if (base_match.size() >= 2) {
      return (core::index4) atoi(base_match[1].str().c_str());
    }
  }

  return 0;
}

bool contains_aa_only(const std::string &sequence) {

  static const std::string accepted_letters = "ARNDCEQGHILKMFPSTWYVX";
  for (const char & c : sequence)
    if (accepted_letters.find(c) == std::string::npos) return false;
  return true;
}

core::index2 length_without_gaps(const std::string & sequence) {
  return std::count_if(sequence.begin(), sequence.end(), [](char c) {return (c!= '-')&&(c !='_');});
}

void remove_gaps(std::string & sequence_string) {

  sequence_string.erase(std::remove(sequence_string.begin(), sequence_string.end(), '-'), sequence_string.end());
  sequence_string.erase(std::remove(sequence_string.begin(), sequence_string.end(), '_'), sequence_string.end());
}

std::string remove_gaps(const std::string & sequence_string) {

  std::string out(sequence_string);
  out.erase(std::remove(out.begin(), out.end(), '-'), out.end());
  out.erase(std::remove(out.begin(), out.end(), '_'), out.end());
  return out;
}

core::index2 count_aa(const std::string &seq) {
  static const std::string aa = "ARNDCEQGHILKMFPSTWYVX";
  return count_if(seq.begin(), seq.end(), [](const char c) { return aa.find(c) != std::string::npos; });
}

core::index2 count_aa(const Sequence &seq) {
  core::index2 cnt = 0;
  for (core::index2 i = 0; i < seq.length(); ++i)
    cnt += seq.get_monomer(i).type == 'P';
  return cnt;
}

core::index2 count_nt(const std::string &seq) {
  static const std::string nt = "actguACTGU";
  return count_if(seq.begin(), seq.end(), [](const char c) { return nt.find(c) != std::string::npos; });
}

core::index2 count_nt(const Sequence &seq) {
  core::index2 cnt = 0;
  for (core::index2 i = 0; i < seq.length(); ++i) cnt += seq.get_monomer(i).type == 'N';
  return cnt;
}

std::vector<core::index4> & count_each_aa(const std::string &seq, std::vector<core::index4> & counts) {

  if(counts.size()<21) counts.resize(21);

  for(const char c : seq) {
    index2 id = core::chemical::Monomer::get(c).id;
    if(id<21) ++counts[id];
  }

  return counts;
}

std::vector<core::index4> & count_each_aa(const Sequence &seq, std::vector<core::index4> & counts) {

  return count_each_aa(seq.sequence,counts);
}

core::index1 count_aa_types(const std::string &seq) {

  std::unordered_set<char> s;
  for (const char c:seq) s.insert(c);

  return s.size();
}

core::index1 count_aa_types(const Sequence &seq) {

  return count_aa_types(seq.sequence);
}

}
}
}
