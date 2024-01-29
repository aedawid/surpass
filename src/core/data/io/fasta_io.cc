#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <fcntl.h>

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <stdexcept>

#include <core/data/io/fasta_io.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <utils/io_utils.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>
#include <core/alignment/on_alignment_computations.hh>

namespace core {
namespace data {
namespace io {

using core::data::sequence::Sequence;
using core::data::sequence::SecondaryStructure;

static utils::Logger fasta_io_logs("fasta_io");

std::vector<std::shared_ptr<Sequence>> & read_fasta_file(const std::string file_name,
    std::vector<std::shared_ptr<Sequence>> & sink) {

  size_t n = 0;

  fasta_io_logs << utils::LogLevel::FILE << "Reading FASTA from " << file_name << " ...\n";
  std::ifstream infile;
  utils::in_stream(file_name,infile);
  std::string line, current_sequence = "", current_id = "empty";

  core::data::sequence::Sequence_SP seq = nullptr;
  infile >> seq;
  while(seq!=nullptr) {
    sink.push_back(seq);
    infile >> seq;
    ++n;
  }

  fasta_io_logs << utils::LogLevel::FILE << "found " << n << " sequences\n";

  return sink;
}

std::istream & operator>>(std::istream & cin, core::data::sequence::Sequence_SP & seq) {

  seq = nullptr;
  if (cin.eof()) return cin;
  std::string header, seqstr, line;
  while (std::getline(cin, header)) {
    utils::trim(header);
    if (header[0] == '>') break;
  }

  if (header.length() == 0) return cin; // no sequence found because the header is empty! Returning nullptr.
  header = header.substr(1, header.size());
  while ((cin.peek() != '>') && (!cin.eof())) {
    std::getline(cin, line);
    utils::trim(line);
    if (line[0] != '#') seqstr += line;
  }

  std::replace_if(seqstr.begin(), seqstr.end(), [](const char c) { return (((c < 'A') || (c > 'Z')) && (c != '-') && (c != '_')); }, 'X');
  if (seqstr.size() > 0) seq = std::shared_ptr<Sequence>(new core::data::sequence::Sequence(header, seqstr, 0));

  return cin;
}

std::string create_fasta_string(const std::string & header, const std::string & sequence,
    const core::index2 line_width, const bool include_header) {

  std::string fasta = (include_header) ? ">" + header + "\n" : "";
  core::index2 len = sequence.length();
  core::index2 start = 0;
  do {
    core::index2 n = std::min(line_width, (core::index2) (len - start));
    fasta += sequence.substr(start, n) + "\n";
    start += n;
  } while (start < len);

  return fasta;
}

std::string create_fasta_string(const Sequence & seq, const core::index2 line_width, const bool include_header) {

  return create_fasta_string(seq.header(),seq.sequence,line_width);
}

std::string create_fasta_secondary_string(const SecondaryStructure & sec_str, const core::index2 line_width, const bool include_header) {

  std::string fasta = (include_header) ? ">" + sec_str.header() + " - secondary structure\n" : "";
  core::index2 len = sec_str.sequence.length();
  core::index2 start = 0;
  do {
    core::index2 n = std::min(line_width, (core::index2) (len - start));
    fasta += sec_str.str().substr(start, n) + "\n";
    start += n;
  } while (start < len);

  return fasta;
}

std::string create_fasta_string(const core::alignment::PairwiseAlignment &ali, const Sequence &query_sequence,
                                const Sequence &tmplt_sequence, const core::index2 line_width) {

  std::stringstream out;
  const std::string s1a = ali.get_aligned_query(query_sequence.sequence, '-');
  const std::string s2a = ali.get_aligned_template(tmplt_sequence.sequence, '-');
  core::index2 n_id = core::alignment::sum_identical(s1a, s2a);
  out << utils::string_format("# %9.2f %4d %6.2f", ali.alignment_score, n_id,
    core::real(n_id) * 2.0 / (ali.template_length() + ali.query_length())) << "\n";
  out << create_fasta_string(query_sequence.header(), s1a, line_width);
  out << create_fasta_string(tmplt_sequence.header(), s2a, line_width);

  return out.str();
}


}
}
}
