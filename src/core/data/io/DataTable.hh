#ifndef CORE_DATA_IO_DataTable_H
#define CORE_DATA_IO_DataTable_H

#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>

#include <utils/string_utils.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace io {

/** @brief Represents a single row of data
 */
class TableRow : public std::vector<std::string> {
public:

  /// Returns a single entity from this row of data
  template<typename T>
  inline T get(const index2 index) { return utils::from_string<T>(operator [](index)); }

  /// Returns a single entity from this row of data
  template<typename T>
  inline T get(const index2 index) const { return utils::from_string<T>(operator [](index)); }
};

/** @brief Operator to print a data table row to a stream
 *
 * @param out - output stream
 * @param row - a table row instance
 * @return a reference to the given output stream
 */
std::ostream &operator<<(std::ostream &out, const TableRow &row);

/** Reads and stores data represented as a flat table.
 *
 * The data is read and stored row-wise. Each row from a file is stored
 * as a TableRow instance on a std::vector. This class provides also means to access data by-column
 *
 * todo_example Provide an example that read a column and e.g. computes basic statistics from it (avg, sdev, min, max, etc)
 */
class DataTable : public std::vector<TableRow> {
public:

  /// Create an empty data container
  DataTable() : logger("DataTable") {}

  /** Create a container and read data from a file into it.
   *
   * @param fname - name of the file with input data
   */
  DataTable(const std::string &fname) : logger("DataTable") { load(fname); }

  /** Loads new data into this container
   *
   * @param fname - name of the file with input data
   */
  void load(const std::string & fname);

  /** Loads new data into this container
   *
   * @param source - input stream
   */
  void load(std::istream & source);

  /** @brief Copies data from a single column into a given vector.
   *
   * All data in the column must be of the same type <code>T</code>
   * @param index - zero-referenced index of a column in this DataTable.
   * @param destination - vector where the data will be stored
   * @tparam T - the type of data stored in the requested column
   */
  template<typename T>
  std::vector<T> & column(const index2 index, std::vector<T> & destination) {
    std::transform(begin(), end(), std::back_inserter(destination),
        [index](const TableRow& r) {return utils::from_string<T>(r[index]);});
    return destination;
  }

private:
  utils::Logger logger;
};

}
}
}

#endif
