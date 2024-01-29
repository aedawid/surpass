/** @file SimilarityMatrix.hh 
 * @brief Holds substitution / similarity matrix in a generic container.
 * Template speciacions provided for primitive data types
 */
#ifndef CORE_ALIGNMENT_SCORING_SimilarityMatrix_H
#define CORE_ALIGNMENT_SCORING_SimilarityMatrix_H

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <unordered_map>

#include <core/index.hh>
#include <core/SURPASSenvironment.hh>
#include <core/data/io/DataTable.hh>
#include <core/data/basic/Array2D.hh>

namespace core {
namespace alignment {
namespace scoring {

/** \brief Represents a substitution (aka substitution) matrix.
 *
 * The matrix data may be loaded into an object with <code>from_ncbi_file(const std::string)</code> method.
 * For convenience, this object provides also some basic substitution matrices
 */
template<typename T>
class SimilarityMatrix : public core::data::basic::Array2D<T> {
public:

  /// Creates an empty matrix with a predefined amino acid order (i.e. amino acid symbols assigned to columns)
  SimilarityMatrix(const std::vector<char> & aa);

  /// Create a new similarity matrix with a predefined amino acid order and fill it with data
  SimilarityMatrix(const core::data::basic::Array2D<T> & similarity_matrix, const std::vector<char> & character_mapping);

  /** \brief Prints this matrix into a given stream
   *
   * @param format - a formatting string used to print each matrix entry
   * @param total_field_width - width of each matrix entry string
   * @param out_stream - where to write the output
   */
  void print(const std::string &format,const int total_field_width, std::ostream &out_stream) const;

  /** \brief Converts a matrix index to a character scored
   *
   * @param index - matrix column (or row since it must be symmetric) index
   */
  char index_to_symbol(const index1 index) const { return  index_to_symbol_[index]; }

  /** \brief Converts a character (e.g. amino acid in 1-letter code) to a matrix column (row) index.
   *
   * If the given character is not defined in a substitution matrix, symbol that corresponds to residue UNK is returned.
   * @param symbol - single character denoting the scored residue
   */
  index1 symbol_to_index(const char symbol) const { return  (symbol_to_index_.find(symbol)==symbol_to_index_.end()) ?  symbol_to_index_.at('X') : symbol_to_index_.at(symbol); }

  /** @brief Returns amino acid symbols in the order they are assigned to columns  of this substitution matrix.
   */
  const std::vector<char> & aa_order() const { return index_to_symbol_; }

  /** \brief Loads the matrix from a file in NCBI format
   *
   * @param file_name - file name
   */
  static std::shared_ptr<SimilarityMatrix<T>> from_ncbi_file(const std::string & file_name);

  /** \brief Loads the matrix from data stream in NCBI format
   *
   * @param input - input stream with the data to be parsed
   */
  static std::shared_ptr<SimilarityMatrix<T>>  from_ncbi_stream(std::istream & input);

private:
  std::vector<char> index_to_symbol_;
  std::unordered_map<char, index2> symbol_to_index_;
};

/** @brief Defines NcbiSimilarityMatrix which is a <code>SimilarityMatrix</code> with entries of type <code>short int</code>
 * This is the type of BLOSUM62 or PAM250 matrix. When you need your scores to be floating point values, use
 * <code>SimilarityMatrixReal</code> type
 */
typedef SimilarityMatrix<short int> NcbiSimilarityMatrix;

/// Defines a shared pointer to <code>NcbiSimilarityMatrix</code>
typedef std::shared_ptr<NcbiSimilarityMatrix> NcbiSimilarityMatrix_SP;

/// Defines a <code>SimilarityMatrix</code> with floating point elements
typedef SimilarityMatrix<core::real> SimilarityMatrixReal;

/// Defines a shared pointer to <code>SimilarityMatrixReal</code>
typedef std::shared_ptr<SimilarityMatrixReal> SimilarityMatrixReal_SP;

/** @brief Converts NcbiSimilarityMatrix object to SimilarityMatrixReal.
 *
 * A new matrix that holds values of a type <code>core:;real</code> will be allocated and filled with data from the source
 * If the name does not refer to any matrix known to this factory, a null pointer is returned.
 * @param m - source matrix
 * @return newly allocated matrix of the new type
 */
SimilarityMatrixReal_SP create_real_similarity_matrix(const NcbiSimilarityMatrix & m);

}
}
}

#endif

