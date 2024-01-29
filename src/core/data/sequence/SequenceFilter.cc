#include <string>

#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SequenceFilter.hh>

namespace core {
namespace data {
namespace sequence {

utils::Logger FindInSequenceName::logger = utils::Logger("FindInSequenceName");

void FindInSequenceName::substrings_from_file(const std::string & file_name) {

  std::ifstream infile;
  infile.open(file_name);
  if (infile.is_open()) {
    substrings_from_stream(infile);
    infile.close();
  } else {
    logger << utils::LogLevel::SEVERE << "Can't open an input file: " << file_name << "\n";
  }
  logger << utils::LogLevel::FILE << wanted_substrings.size() << " strings for sequence filtering obtained from "
      << file_name << "\n";
}

void FindInSequenceName::substrings_from_stream(std::istream &in_stream) {

  std::string substr;
  std::string prev = "";
  wanted_substrings.clear();
  while (!in_stream.eof()) {
    in_stream >> substr;
    if (substr.length() > 3 && substr.compare(prev) != 0) {
      wanted_substrings.push_back(substr);
      prev = substr;
    }
  }
}

bool ByUnknownRatio::operator()(const Sequence & s) const {
  real ratio = real(s.length());
  for(core::index4 i=0;i<s.length();++i)
    ratio -= ((s.get_monomer(i) == core::chemical::Monomer::UNK) ||
              (s.get_monomer(i) == core::chemical::Monomer::UNG) ||
              (s.get_monomer(i) == core::chemical::Monomer::UNL));

  return ((ratio / real(s.length())) > known_fraction_);
}

bool ByUnknownRatio::operator()(const std::string & header, const std::string & sequence) const {

  real ratio = std::count(sequence.cbegin(),sequence.cend(),'X');
  ratio = sequence.length() - ratio;
  return ((ratio / real(sequence.length())) > known_fraction_);
}

}
}
}
