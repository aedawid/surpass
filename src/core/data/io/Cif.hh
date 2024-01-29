#ifndef CORE_DATA_IO_CIF_H
#define CORE_DATA_IO_CIF_H

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>

#include <utils/string_utils.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace io {

/** @brief Represents a <code>data</code> block of a CIF file
 *
 * @include ex_Cif.cc
 */
class DataBlock {
public:
  const std::string block_name; ///< Name of the block

  /// Prints a given CIF <code>data</code> block  nicely in CIF format
  friend std::ostream &operator<<(std::ostream &out, const DataBlock &e);

  /// Prints a given CIF <code>data</code> block  nicely in CIF format
  friend std::ostream &operator<<(std::ostream &out, const DataBlock *e);

  /** @brief Stores key-value pairs found in a <code>data</code> block
   */
  std::map<std::string, std::string> tokens;

  /// Creates an empty block of data
  DataBlock(const std::string &name) : block_name(name) {}

private:
  inline size_t longest_key_length() const {
    size_t m = 0;
    for (const std::pair<std::string, std::string> &k : tokens)
      m = std::max(m, k.first.length());
    return m;
  }
};

/** @brief For convenience, define a shortcut for a shared pointer to the DataBlock type
 */
typedef std::shared_ptr<DataBlock>  DataBlock_SP;

/** @brief Simple CIF file reader.
 *
 * @include ex_Cif.cc
 */
class Cif {
public:
  /** @brief Holds all data blocks found in a file
   *
   */
  std::vector<std::shared_ptr<DataBlock> > blocks;

  /** @brief Constructor opens CIF file for reading.
   *
   * No data is loaded by this constructor. Use read_block() to read a single CIF block of data
   * @param src - input stream
   */
  Cif(const std::string &file_name) {
    source = std::shared_ptr<std::istream>(new std::ifstream(file_name));
    logger << utils::LogLevel::FILE << "Reading CIF data from " << file_name
           << "\n";
  }

  /** @brief Constructor prepares CIF data stream for reading.
   *
   * No data is loaded by this constructor. Use read_block() to read a single CIF block of data
   * @param src - input stream
   */
  Cif(std::shared_ptr<std::istream> src) {
    source = src;
    logger << utils::LogLevel::FILE << "Reading CIF data from a stream\n";
  }

  /** @brief Reads the first block of data found in the file (or stream) opened by this object
   *
   * @return a data block object
   */
  std::shared_ptr<DataBlock> read_block();

private:
  std::string last_line;
  static utils::Logger logger;
  std::shared_ptr<std::istream> source;
};

}
}
}
/**
 * @example  ex_Cif.cc
 */
#endif
