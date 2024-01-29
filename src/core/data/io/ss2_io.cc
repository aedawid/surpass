#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <utils/io_utils.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/io/ss2_io.hh>

namespace core {
namespace data {
namespace io {

static utils::Logger ss2_io_logs("ss2_io");

core::data::sequence::SecondaryStructure_SP read_ss2(std::istream & in_stream, const std::string & header) {

  std::string line;
  std::getline(in_stream, line); // Skip header

  int pos;
  int first_pos = -1000;
  char aa, ss;
  core::real h, e, c;
  std::vector<core::real> hv, ev, cv;
  std::string seq;
  while (std::getline(in_stream, line)) {
    if (line.length() < 15) continue;
    if (line[0] == '#') continue;
    std::stringstream ls(line);
    ls >> pos >> aa >> ss >> c >> h >> e;
    seq += aa;
    if (first_pos == -1000) first_pos = pos;
    hv.push_back(h);
    ev.push_back(e);
    cv.push_back(c);
  }
  sequence::SecondaryStructure_SP sec_str = std::make_shared<sequence::SecondaryStructure>(header, seq, first_pos);
  for (core::index2 i = 0; i < seq.size(); ++i)
    sec_str->fractions(i, hv[i], ev[i], cv[i]);

  ss2_io_logs << utils::LogLevel::FILE << "Secondary structure loaded: " << sec_str->str() << "\n";

  return sec_str;
}

core::data::sequence::SecondaryStructure_SP read_ss2(const std::string & in_fname, const std::string & header) {

  if(!utils::if_file_exists(in_fname)) {
    std::string msg = "Can't find SS2 file: " + in_fname + "!\n";
    ss2_io_logs<<utils::LogLevel::SEVERE<<msg;
    throw std::runtime_error(msg);
  }
  ss2_io_logs << utils::LogLevel::FILE << "Reading SS2 data from " << in_fname << "...\n";

  std::ifstream stream(in_fname);
  return read_ss2(stream, header);
}

std::vector<std::shared_ptr<core::data::sequence::Sequence>> & read_ss2(const std::string & in_fname,
    std::vector<std::shared_ptr<core::data::sequence::Sequence>> & sink, const std::string & header) {
  sink.push_back(read_ss2(in_fname, header));
  return sink;
}

void write_ss2(const core::data::sequence::SecondaryStructure & ss, std::ostream & out) {

  out << "# PSIPRED VFORMAT\n\n";
  for(core::index2 i=0;i<ss.length();++i)
    out
        << utils::string_format("%4d %c %c   %5.3f  %5.3f  %5.3f\n", ss.first_pos() + i, ss[i], ss.ss(i),
            ss.fraction_C(i), ss.fraction_H(i), ss.fraction_E(i));
}

}
}
}

