#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>

#include <core/algorithms/basic_algorithms.hh>
#include <utils/string_utils.hh>
#include <core/data/io/DataTable.hh>
#include <utils/Logger.hh>
#include <utils/io_utils.hh>

namespace core {
namespace data {
namespace io {

std::ostream &operator<<(std::ostream &out, const TableRow &row) {

  for (std::size_t i = 0; i < row.size() - 1; ++i) {out << row[i] << " ";}
  out << row.back();

  return out;
}

void DataTable::load(std::istream & source) {

  std::string line;
  while (std::getline(source, line)) {
    if (line.length() < 1) continue;
    utils::trim(line);
    if (line[0] == '#') continue;
    TableRow row;
    push_back(row);
    utils::split(line, back());
  }
  logger << utils::LogLevel::FILE << "found " << size() << " rows\n";
}

void DataTable::load(const std::string & fname) {

  if(!utils::if_file_exists(fname)) {
    std::string msg = "Can't find a file: " + fname + "!\n";
    logger<<utils::LogLevel::SEVERE<<msg;
    throw std::runtime_error(msg);
  }

  logger << utils::LogLevel::FILE << "loading tabular data from " << fname << "\n";
  std::ifstream in(fname);
  if (!in) std::runtime_error("Can't read from the file: "+fname);
  load(in);
  in.close();
}

}
}
}
