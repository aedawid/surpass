#include <cmath>
#include <string>

#include <core/index.hh>
#include <core/alignment/on_alignment_computations.hh>

namespace core {
namespace alignment {

//using core::data::sequence::Sequence;

static utils::Logger on_alignment_computations_logger("on_alignment_computations");

core::index2 sum_identical(const std::string &s1, const std::string &s2) {

  core::index2 match = 0;
  std::string::const_iterator it1 = s1.begin();
  std::string::const_iterator it2 = s2.begin();
  while (it1 != s1.end()) {
    if ((*it1 == *it2) && (*it1 != '-') && (*it1 != '_')) ++match;
    ++it1;
    ++it2;
  }
  return match;
}

core::index2 sum_identical(const Sequence &s1, const Sequence &s2) { return sum_identical(s1.sequence, s2.sequence); }

core::index2 sum_identical(const PairwiseAlignment &alignment, const std::string &query_sequence,
                           const std::string &tmplate_sequence) {

  const std::string tq = alignment.get_aligned_query(query_sequence, '-');
  const std::string tt = alignment.get_aligned_template(tmplate_sequence, '-');

  return sum_identical(tq,tt);
}

core::index2 sum_identical(const PairwiseSequenceAlignment &alignment) {

  return sum_identical(*alignment.alignment, alignment.query_sequence->sequence, alignment.template_sequence->sequence);
}

core::index2 sum_aligned(const std::string &s1, const std::string &s2) {

  core::index2 match = 0;
  std::string::const_iterator it1 = s1.begin();
  std::string::const_iterator it2 = s2.begin();
  while (it1 != s1.end()) {
    if ((*it1 != '-') && (*it1 != '_') && (*it2 != '-') && (*it2 != '_')) ++match;
    ++it1;
    ++it2;
  }

  return match;
}

core::index2 sum_aligned(const Sequence &s1, const Sequence &s2) { return sum_aligned(s1.sequence, s2.sequence); }

core::index2 sum_aligned(const PairwiseAlignment &alignment) {

  std::string p = alignment.to_path();
  return std::count(p.begin(), p.end(), '*');
}

core::index2 compute_aln_n(const PairwiseAlignment &reference_alignment, const core::index2 k,
                           const PairwiseAlignment &another_alignment) {

  if (another_alignment.query_length() != reference_alignment.query_length()) {
    on_alignment_computations_logger << utils::LogLevel::WARNING
                                     << "Cannot compare alignment to the reference because the query sequences differ in length\n";
    return 0;
  }

  size_t sum = 0;
  for (size_t i = 0; i < reference_alignment.query_length(); i++) {
    const int i_ref = reference_alignment.which_template_for_query(i);
    const int i_other = another_alignment.which_template_for_query(i);
    if ((i_ref < 0) || (i_other < 0))
      continue;
    if (abs(i_ref - i_other) <= k)
      ++sum;
  }

  return sum;
}


short int calculate_score(const core::data::sequence::Sequence &s1, const core::data::sequence::Sequence &s2,
                             const scoring::NcbiSimilarityMatrix & scoring, short int gap_open, short int gap_extend) {
  return calculate_score(s1.sequence, s2.sequence, scoring, gap_open, gap_extend);
}

short int calculate_score(const std::string &s1, const std::string &s2,
                             const scoring::NcbiSimilarityMatrix & scoring, short int gap_open, short int gap_extend) {

  short int score = 0;
  bool is_gap_open = false;
  for (core::index4 ipos = 0; ipos < s1.length(); ++ipos) {
    if ((s1[ipos] == '-') || (s1[ipos] == '_') || (s2[ipos] == '-') || (s2[ipos] == '_')) {
      if (is_gap_open) score += gap_extend;
      else {
        score += gap_open;
        is_gap_open = true;
      }
    } else {
      is_gap_open = false;
      score += scoring(scoring.symbol_to_index(s1[ipos]), scoring.symbol_to_index(s2[ipos]));
    }
  }

  return score;
}


}
}
