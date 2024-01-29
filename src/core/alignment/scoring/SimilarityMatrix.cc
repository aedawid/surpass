#include <vector>
#include <memory>

#include <core/index.hh>
#include <core/data/basic/Array2D.hh>
#include <core/alignment/scoring/SimilarityMatrix.hh>
#include <core/SURPASSenvironment.hh>
#include <core/data/io/DataTable.hh>

#include <utils/io_utils.hh>
#include <utils/Logger.hh>

namespace core {
namespace alignment {
namespace scoring {

SimilarityMatrixReal_SP create_real_similarity_matrix(const NcbiSimilarityMatrix & m) {

  core::data::basic::Array2D<core::real> data(m.aa_order().size(),m.aa_order().size());
  SimilarityMatrixReal_SP mm = std::make_shared<SimilarityMatrixReal>(data, m.aa_order());
  return mm;
}

template<typename T>
SimilarityMatrix<T>::SimilarityMatrix(const std::vector<char> & aa) :
  core::data::basic::Array2D<T>(aa.size(), aa.size()), index_to_symbol_(aa.size()) {
  for (index1 i = 0; i < aa.size(); ++i) {
    index_to_symbol_[i] = aa[i];
    symbol_to_index_[aa[i]] = i;
  }
}

template<typename T>
SimilarityMatrix<T>::SimilarityMatrix(const core::data::basic::Array2D<T> & similarity_matrix, const std::vector<char> & character_mapping) :
  SimilarityMatrix(character_mapping) {

  for(core::index2 i=0;i<character_mapping.size();++i)
    core::data::basic::Array2D<T>::operator[](i) = similarity_matrix[i];
}

template<typename T>
std::shared_ptr<SimilarityMatrix<T>> SimilarityMatrix<T>::from_ncbi_file(const std::string & file_name) {

  utils::Logger logs("from_ncbi_file");
  std::string fname = core::SURPASSenvironment::from_file_or_db(file_name, "alignments");
  logs << utils::LogLevel::FILE << "Loading a substitution matrix from an external file: " << fname << "\n";
  std::ifstream in(fname);
  return from_ncbi_stream(in);
}

template<typename T>
std::shared_ptr<SimilarityMatrix<T>> SimilarityMatrix<T>::from_ncbi_stream(std::istream & input) {

  utils::Logger logs("from_ncbi_stream");
  core::data::io::DataTable dt;
  dt.load(input);

  if ((dt.size() != 26) || (dt[0].size() != 25)) {
    logs << utils::LogLevel::CRITICAL << "Incorrect substitution matrix format\n";
    return nullptr;
  }

  index1 m_size = 25;

  std::vector<char> aa_order;
  for (index1 i = 0; i < m_size; ++i) aa_order.push_back(dt[0][i][0]);
  std::shared_ptr<SimilarityMatrix<T>> m = std::make_shared<SimilarityMatrix<T>>(aa_order);

  for (index1 i = 0; i < m_size; ++i)
    for (index1 j = 0; j < m_size; ++j)
      m->set(i, j, dt[i + 1].get<T>(j + 1));

  return m;
}

template<typename T>
void SimilarityMatrix<T>::print(const std::string &format, const int total_field_width, std::ostream &out_stream)  const {

  int k = -1;
  out_stream << ' ';
  for (size_t j = 0; j < core::data::basic::Array2D<T>::n_columns; j++)
    out_stream << std::setw(total_field_width) << index_to_symbol_[j];
  out_stream << "\n";

  for (size_t i = 0; i < core::data::basic::Array2D<T>::n_rows; i++) {
    out_stream << index_to_symbol_[i];
    for (size_t j = 0; j < core::data::basic::Array2D<T>::n_columns; j++)
      out_stream << utils::string_format(format, core::data::basic::Array2D<T>::the_data[++k]);
    out_stream << std::endl;
  }
}

template class SimilarityMatrix<short>;
template class SimilarityMatrix<int>;
template class SimilarityMatrix<float>;
template class SimilarityMatrix<double>;

}
}
}

