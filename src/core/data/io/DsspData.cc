#include <string>
#include <algorithm>
#include <iostream>
#include <vector>

#include <core/data/io/DsspData.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>
#include <utils/io_utils.hh>

namespace core {
namespace data {
namespace io {

const  std::vector<int> DsspDataLine::start_column = { 1, 7, 12, 14, 17, 27, 31, 35, 40, 47, 52, 58, 63, 69, 74, 80, 86,
    92, 98, 104, 110, 116, 123, 130, 11};
const  std::vector<int> DsspDataLine::end_column = { 6, 11, 13, 15, 18, 30, 34, 39, 46, 51, 57, 62, 68, 73, 79, 84, 92,
    98, 104, 110, 116, 123, 130, 137, 12};

const std::vector<int> DsspDataLine::start_column_classic = { 1, 7, 12, 14, 17, 27, 31, 35, 40, 45, 49, 54, 58, 63, 67,
    72, 76, 84, 90, 96, 102, 108, 115, 122, 11};

const  std::vector<int> DsspDataLine::end_column_classic = { 6, 11, 13, 15, 18, 30, 34, 39, 44, 49, 53, 58, 62, 67, 71,
    76, 84, 90, 96, 102, 108, 115, 120, 129, 12};
const int DsspDataLine::offset = 1;

DsspDataLine::DsspDataLine(const std::string & line,const bool ifClassicFormat) {

  const std::vector<int> & start = (ifClassicFormat) ? start_column_classic : start_column;
  const std::vector<int> & stop = (ifClassicFormat) ? end_column_classic : end_column;

  int i_column = 0;
  position = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;

  seq_position = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  chain_letter = line[start[i_column] - offset];
  i_column++;
  residue_letter = line[start[i_column] - offset];
  i_column++;
//    if (Character.isLowerCase(residueLetter))
//      residueLetter = 'C';
  structure = line[start[i_column] - offset];
  i_column++;
  bp1 = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  bp2 = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  acc = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  nho_partner1 = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  nho_energy1 = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  ohn_partner1 = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  ohn_energy1 = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  nho_partner2 = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  nho_energy2 = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  ohn_partner2 = utils::from_string<int>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  ohn_energy2 = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;

  tco = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  kappa = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  alpha = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  phi = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  psi = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  xCa = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  yCa = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  zCa = utils::from_string<double>(line, start[i_column] - offset, stop[i_column] - offset, -1);
  i_column++;
  icode =  line[start[i_column] - offset];
  i_column++;
}


bool DsspDataLine::is_my_residue(const core::data::structural::Residue & r) const {

  if((seq_position==r.id())&&(icode==r.icode())&&(chain_letter==r.owner()->id())) return true;
  return false;
}


DsspData::DsspData(const std::string & file_name, const bool classic_format) : if_classic_format(classic_format), logger("DsspData") {

  if (!utils::if_file_exists(file_name)) {
    std::string msg = "Can't find a file: " + file_name + "!\n";
    logger << utils::LogLevel::SEVERE << msg;
    throw std::runtime_error(msg);
  }

  std::ifstream in(file_name);
  if (!in) std::runtime_error("Can't read from the file: " + file_name);
  file_name_ = file_name;
  load(in, classic_format);
  in.close();
}


void DsspData::load(std::ifstream & in, const bool classic_format) {

  if (attempts < 2) {
    if (classic_format)
      logger << utils::LogLevel::FILE << "Attempting CLASSIC file format\n";
    else
      logger << utils::LogLevel::FILE << "Attempting NEW file format\n";
  } else {
    logger << utils::LogLevel::SEVERE << "No more file format variants to try - is it really a DSSP file?\n";
  }

  attempts++;
  std::string line = "";
  try {
    do {
      std::getline(in, line);
      if(line.find("HEADER")==0) {
        std::vector<std::string> tokens = utils::split(utils::trim(line),' ');
        code_ = tokens.back();
        if ((code_.length() < 2) && (tokens.size() > 1)) code_ = tokens[tokens.size() - 2];
      }
    } while (line[2] != '#');
    while (std::getline(in, line)) {
      if (line.length() < 120) continue;
      if (line[13] == '!') continue;
      data.emplace_back(line, classic_format);
      DsspDataLine & d = data.back();
      if (d.position > highest_residue_index) highest_residue_index = d.position;
    }
    return;
  } catch (std::ifstream::failure & e) {
    if ((!in.fail()) || (!in.eof()) || (in.bad()))
      logger << utils::LogLevel::CRITICAL << "Exception captured: " << e.what() << "\n";
  } catch (std::exception & e) {

  }
//  load(filename,!ifClassicFormat);
}

/** \brief Converts DSSP class into 3-state classification
 * @param dsspChar - one of the following: E, B, H, I, G, T, S
 * @return H, E or C
 */
char DsspData::dssp_to_hec(const char dsspChar) {

  switch (dsspChar) { // Convert sec. struct into HEC code
  case 'E':
    return 'E';
  case 'H':
  case 'I':
  case 'G':
    return 'H';
  case 'B':
  case 'T':
  case 'S':
  case ' ':
    return 'C';
  }
  return 'C';
}

std::vector<core::data::sequence::SecondaryStructure> DsspData::create_sequences() const {

  std::vector<core::data::sequence::SecondaryStructure> output;
  std::map<char,std::string> seq,ss2;
  std::map<char,int> first_pos;
  for(const auto & d : data) {
    char chain = d.chain_letter;
    if(seq.find(chain)==seq.end()) {
      seq[chain] = "";
      ss2[chain] = "";
      first_pos[chain] = d.seq_position;
    }
    seq[chain] += d.residue_letter;
    ss2[chain] += dssp_to_hec(d.structure);
  }
  for(const auto & p : seq) {
    char chain = p.first;
    if (chain != ' ') output.emplace_back(utils::string_format("%s chain %c", code_.c_str(), chain), seq[chain],
        first_pos[chain], ss2[chain]);
    else output.emplace_back(code_, seq[chain], first_pos[chain], ss2[chain]);
  }
  return output;
}

}
}
}
