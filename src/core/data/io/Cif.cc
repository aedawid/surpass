#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>

#include <core/data/io/Cif.hh>
#include <utils/string_utils.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace io {

utils::Logger Cif::logger = utils::Logger("Cif");

std::ostream &operator<<(std::ostream &out, const DataBlock &b) {

  out << "data_" << b.block_name << "\n";
  size_t len = b.longest_key_length();
  std::string fmt = "%-" + utils::string_format("%ds %%s\n", len);
  std::string s;
  for (const std::pair<std::string, std::string> &t : b.tokens) {
    s += utils::string_format(fmt, t.first.c_str(), t.second.c_str());
  }
  out << s;

  return out;
}

std::ostream &operator<<(std::ostream &out, const DataBlock *b) {

  out << "data_" << b->block_name << "\n";
  size_t len = b->longest_key_length();
  std::string fmt = "%-" + utils::string_format("%ds %%s\n", len);
  std::string s;
  for (const std::pair<std::string, std::string> &t : b->tokens) {
    s += utils::string_format(fmt, t.first.c_str(), t.second.c_str());
  }
  out << s;

  return out;
}

std::shared_ptr<DataBlock> Cif::read_block() {

  std::shared_ptr<DataBlock> recent_block = 0;
  bool loop_open = false;
  std::vector<std::string> tmp;

  if (last_line.compare(0, 5, "data_") == 0) {
    recent_block = std::shared_ptr<DataBlock>(new DataBlock(last_line.substr(5)));
    last_line = "";
  }
  std::string line;
  while (std::getline(*source, line)) {
    if ((line[0] == '#') || (line.length() < 5)) {
      loop_open = false;
      continue;
    }
    if (line.compare(0, 5, "data_") == 0) {
      if (recent_block != 0) {
        last_line = line;
        break;
      }
      recent_block = std::shared_ptr<DataBlock>(new DataBlock(line.substr(5)));
      continue;
    }
    if (line.compare(0, 5, "loop_") == 0) {
      loop_open = true;
      continue;
    }
    if (!loop_open) {
      tmp.clear();
      utils::from_string<std::string>(line, 2, tmp);
      if (tmp.size() == 2) recent_block->tokens.insert(std::pair<std::string, std::string>(tmp[0], tmp[1]));
    } else {
      // \todo_code implement loop-type block reading here !!!
    }
  }
  if (recent_block != 0) {
    logger << utils::LogLevel::FINE << "Finished data frame >"
           << recent_block->block_name << "< with "
           << recent_block->tokens.size() << " key:value pairs\n";
    return recent_block;
  }
  logger << utils::LogLevel::FINE << "Finished reading all the data\n";
  return 0;
}

}
}
}

