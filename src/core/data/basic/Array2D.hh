#ifndef CORE_DATA_BASIC_Array2D_H
#define CORE_DATA_BASIC_Array2D_H

#include <utility>
#include <tuple>
#include <vector>
#include <stdexcept>
#include <algorithm>		// for std::min_element() and std::make_tuple()

#include <core/data/basic/Vec3.hh>
#include <core/data/io/DataTable.hh>

#include <utils/string_utils.hh>

namespace core {
namespace data {
namespace basic {

/** @brief 2D matrix implemented with data stored as std::vector<T>
 *
 * @input ex_basic_algebra.cc
 * @tparam T - type of the data stored in the matrix
 */
template<typename T>
class Array2D {
public:
  /** @brief Creates a new matrix.
   *
   * @param n_rows - the number of rows
   * @param n_columns  - the number of columns
   */
  Array2D(const index2 n_rows, const index2 n_columns) : n_rows(n_rows), n_columns(n_columns) {

    if ((n_columns == 0) || (n_rows == 0)) throw std::invalid_argument("Invalid array size");
    the_data.resize(n_rows * n_columns);
    size = n_rows * n_columns;
  }

  /** @brief Creates a new matrix.
   *
   * @param n_rows - the number of rows
   * @param n_columns  - the number of columns
   * @param data - data used to fill the matrix
   */
  Array2D(const index2 n_rows, const index2 n_columns, const std::vector<T> & data) : n_rows(n_rows), n_columns(n_columns) {

    if ((n_columns == 0) || (n_rows == 0)) throw std::invalid_argument("Invalid array size");
    size = n_rows * n_columns;
    if (size != data.size()) throw std::invalid_argument("Array size inconsistent with the input content data");
    the_data.resize(size);
    std::copy(data.begin(), data.end(), the_data.begin());
  }

  /** @brief Creates a new \f$ 3 \times 3 \f$ matrix from given row vectors.
   *
   *  @param r1 - the first row of the new matrix
   *  @param r2 - the second row of the new matrix
   *  @param r3 - the third row of the new matrix
   *  @return a new matrix
   *
   * The resulting square matrix A is defined as
   * \f[
A = \left[
  \begin{array}{c}
r_1\\
r_2\\
r_3\\
  \end{array}
\right]
   * \f]
   */
  static Array2D from_rows(const Vec3 & r1,const Vec3 & r2,const Vec3 & r3) {

    Array2D a(3,3);
    a[0] = r1.x;
    a[1] = r1.y;
    a[2] = r1.z;
    a[3] = r2.x;
    a[4] = r2.y;
    a[5] = r2.z;
    a[6] = r3.x;
    a[7] = r3.y;
    a[8] = r3.z;

    return a;
  }

  /** @brief Creates a new \f$ 3 \times 3 \f$ matrix from given row vectors.
   *
   * @param r1 - the first row of the new matrix
   * @param r2 - the second row of the new matrix
   * @param r3 - the third row of the new matrix
   * @return a new matrix
   *
   * The resulting square matrix A is defined as
   * \f[
A = \left[
c_1, c_2, c_3
\right]
   * \f]
   *
   */
  static Array2D from_columns(const Vec3 & c1,const Vec3 & c2,const Vec3 & c3) {

    Array2D a(3,3);
    a[0] = c1.x;
    a[1] = c2.x;
    a[2] = c3.x;
    a[3] = c1.y;
    a[4] = c2.y;
    a[5] = c3.y;
    a[6] = c1.z;
    a[7] = c2.z;
    a[8] = c3.z;

    return a;
  }

  /** @brief Reads Array2D data from a file.
   *
   * The file must have a flat table format. All data entries from the file will be converted to the type T
   *
   * @param file_name - name of the input file
   * @return  array containing the data
   */
  static Array2D from_file(const std::string & file_name) {

    core::data::io::DataTable dt;
    dt.load(file_name);
    core::index2 nrow = dt.size();
    core::index2 ncol = dt[0].size();

    Array2D a(nrow, ncol);
    for (core::index2 i = 0; i < nrow; ++i) {
      const core::data::io::TableRow &tr = dt[i];
      for (core::index2 j = 0; j < ncol; ++j) a(i, j) = tr.get<T>(j);
    }

    return a;
  }

  /// Returns the number of rows of this matrix
  size_t count_rows() const { return n_rows; }

  /// Returns the number of columns of this matrix
  size_t count_columns() const { return n_columns; }

  /// Returns 1D style index of a particular (row,column) entry
  size_t to1D(const size_t row, const size_t column) const { return row * n_columns + column; }

  /// Returns 2D style (row,column) index of a particular  entry in the internal storage
  std::pair<size_t, size_t> to2D(const size_t index) const { return std::make_pair(index / n_columns, index % n_columns); }

  /// Returns a value stored in the 1D internal vector at a given index
  T const &operator[](const size_t index) const { return the_data[index]; }

  /// Returns a value stored in the 1D internal vector at a given index
  T &operator[](const size_t index) { return the_data[index]; }

  /// Returns a value stored in the certain element of the matrix
  T const &operator()(const size_t row, const size_t column) const { return the_data[to1D(row, column)]; }

  /// Returns a value stored in the certain element of the matrix
  T &operator()(const size_t row, const size_t column) { return the_data[to1D(row, column)]; }

  /// Returns a value stored in the certain element of the matrix
  T &get(const size_t row, const size_t column) { return the_data[to1D(row, column)]; }

  /// Returns a const-value stored in the certain element of the matrix
  const T & get(const size_t row,const size_t column) const { return the_data[to1D(row, column)]; }

  /// Sets a value for the certain element of the matrix
  void set(const size_t row,const  size_t column, T t) { the_data[to1D(row, column)] = t; }

  /** @brief Iterator pointing on the first element of this array.
   * The iterator provides all the data in this array row-by-row, so the total number of iterations will be
   * <code>n_columns * n_rows</code>
   * @return iterator to start loop over the data
   */
  typename std::vector<T>::iterator begin() { return the_data.begin(); }

  /** @brief Iterator pointing behind the last element of this array.
   * The iterator provides all the data in this array row-by-row, so the total number of iterations will be
   * <code>n_columns * n_rows</code>
   * @return iterator to end loop over the data
   */
  typename std::vector<T>::iterator end() { return the_data.end(); }

  /** Adds a given value to an element of this array.
   *
   * This method is provided for efficient accumulation of properties related to a 2D index
   *
   * @param row - row index
   * @param column - column index
   * @param inc_val - the value to be added
   */
  void add(const size_t row,const  size_t column, T inc_val) { the_data[to1D(row, column)] += inc_val; }

  /** @brief Copies the content from the given data vector into this matrix
   *
   * @param data - input data
   * @tparam D - a 1D data container type, e.g. <code>std::vector<double> </code>
   */
  template<typename D>
  void set(const D &data) {
    for(index4 i=0;i<the_data.size();++i) the_data[i] = data[i];
  }

  void set(const T e) {
    for(index4 i=0;i<the_data.size();++i) the_data[i] = e;
  }

  /** @brief Clears the matrix with the given zero-value.
   *
   * @param data - value used to fill the whole matrix
   */
  void clear(const T data) { for (index4 i=0;i<size;++i) the_data[i] = data; }

  /// Exposes the internal vector of data
  const std::vector<T> &expose_vector() const { return the_data; }

  /** @brief Returns the minimum value stroded in this matrix
   *
   * @return a tuple of three elements: row index, column index, and the minimum value found
   */
  std::tuple<int, int, T> min() {

    typename std::vector<T>::const_iterator min = std::min_element(the_data.begin(), the_data.end());
    int index = min - the_data.begin();
    return std::make_tuple((int) (index / n_columns), (int) (index % n_columns), *min);
  }

  /** @brief prints the matrix nicely
   *
   * @param format - format used to print a single matrix value
   * @param out_stream - stream to send the data
   */
  void print(const std::string &format, std::ostream &out_stream) {

    int k = -1;
	for (index2 i = 0; i < n_rows; i++) {
	  for (index2 j = 0; j < n_columns; j++) out_stream << utils::string_format(format, the_data[++k]);
	  out_stream << std::endl;
	}
  }

protected:
  std::vector<T> the_data;
  const index2 n_rows;
  const index2 n_columns;
  index4 size;

  size_t row_starts(const size_t row) const { return row * n_columns; }

  size_t row_ends(const size_t row) const { return row * (n_columns) + n_columns - 1; }
};

}
}
}
/**
 * @include  ex_basic_algebra.cc
 */
#endif
