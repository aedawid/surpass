#include <vector>
#include <string>
#include <tuple>
#include <cctype>

#include <core/index.hh>
#include <core/alignment/AlignmentRow.hh>
#include <core/alignment/AlignmentBlock.hh>
#include <core/alignment/PairwiseAlignment.hh>

namespace core {
namespace alignment {

PairwiseAlignment::PairwiseAlignment(const core::index2 first_query_residue,
    const core::index2 first_template_residue,const core::real score,const std::string & alignment_path) :
    alignment_score(score), the_query(first_query_residue), the_template(first_template_residue), logs("PairwiseAlignment") {

  for(const char c : alignment_path) {
    if(c=='*') {
      append_match();
      continue;
    }
    if(c=='|') {
      append_query_gap();
      continue;
    }
    if(c=='-') {
      append_template_gap();
      continue;
    }
  }
}

PairwiseAlignment::PairwiseAlignment(const std::string &  aligned_query,
		const core::index2 first_query_residue, const std::string &  aligned_template,
		const core::index2 first_template_residue,const core::real score, bool expand_unaligned) :
		alignment_score(score), the_query(first_query_residue), the_template(first_template_residue), logs("PairwiseAlignment") {

	core::index2 n_unaligned = 0;
	for (core::index2 i = 0; i < aligned_query.size(); ++i) {
		if ((aligned_query[i] == '-') || (aligned_query[i] == '_')) { // there is a gap in the query sequence
			if ((aligned_template[i] == '-') || (aligned_template[i] == '_')) // both are gapped
				continue;
			clear_unaligned(n_unaligned, expand_unaligned);
			append_query_gap();
		} else { // no gap in the query
			if ((aligned_template[i] == '-') || (aligned_template[i] == '_')) { // but the template is gapped
				clear_unaligned(n_unaligned, expand_unaligned);
				append_template_gap();
			} else {	// a match
				if (islower(aligned_template[i]) && islower(aligned_query[i]))
					n_unaligned++;
				else {
					clear_unaligned(n_unaligned, expand_unaligned);
					append_match();
				}
			}
		}
	}
}

PairwiseAlignment::PairwiseAlignment(const core::data::sequence::Sequence &aligned_query,
									 const core::data::sequence::Sequence &aligned_template, const core::real score) :
	alignment_score(score), the_query(aligned_query.first_pos(), aligned_query.sequence),
	the_template(aligned_template.first_pos(), aligned_template.sequence), logs("PairwiseAlignment") {

	if (aligned_query.length() != aligned_template.length())
		throw std::length_error("Query row of the alignment differs in size from the template!");

}

PairwiseAlignment::PairwiseAlignment(const real score, const AlignmentRow &query, const AlignmentRow &tmplt) :
  alignment_score(score), the_query(first_index(query.first_pos,query,tmplt)),
  the_template(first_index(tmplt.first_pos,tmplt,query)), logs("PairwiseAlignment") {

  n_aligned_ = 0;
  if (the_query.alignment_to_sequence.size() != the_template.alignment_to_sequence.size())
    throw std::length_error("Query row of the alignment differs in size from the template!");

  for (core::index2 i = 0; i < query.length(); ++i) {
    if ((query.position_in_sequence(i) >= 0) && (tmplt.position_in_sequence(i) >= 0)) {
      the_query.append_aligned();
      the_template.append_aligned();
      ++n_aligned_;
      continue;
    }
    if ((query.position_in_sequence(i) >= 0) && (tmplt.position_in_sequence(i) < 0)) {
      the_query.append_aligned();
      the_template.append_gapped();
      continue;
    }
    if ((tmplt.position_in_sequence(i) >= 0) && (query.position_in_sequence(i) < 0)) {
      the_template.append_aligned();
      the_query.append_gapped();
      continue;
    }
  }
}

void PairwiseAlignment::clear_unaligned(core::index2 & n, const bool flag) {

	if (!flag)
		return;
	for (core::index2 i = 0; i < n; i++) {
		the_query.append_gapped();
		the_template.append_aligned();
	}
	for (core::index2 i = 0; i < n; i++) {
		the_template.append_gapped();
		the_query.append_aligned();
	}
	n = 0;
}

std::string PairwiseAlignment::to_path() const {

	std::string out = "";
	for (core::index2 i = 0; i < length(); i++) {
		if (the_query.alignment_to_sequence[i] < 0) {
			out += '|';
			continue;
		}
		if (the_template.alignment_to_sequence[i] < 0) {
			out += '-';
			continue;
		}
		out += '*';
	}

	return out;
}

std::vector<AlignmentBlock> PairwiseAlignment::aligned_blocks() const {

  std::vector<AlignmentBlock> out;
  std::vector<core::index2> tmp;
  if (!is_gapped(0)) {
    tmp.push_back(the_query.alignment_to_sequence[0]);
    tmp.push_back(the_template.alignment_to_sequence[0]);
  }
  for (core::index2 i = 1; i < length(); i++) {
    if ((!is_gapped(i - 1)) && (is_gapped(i))) {
      tmp.push_back(the_query.alignment_to_sequence[i - 1]);
      tmp.push_back(the_template.alignment_to_sequence[i - 1]);
    }
    if ((is_gapped(i - 1)) && (!is_gapped(i))) {
      tmp.push_back(the_query.alignment_to_sequence[i]);
      tmp.push_back(the_template.alignment_to_sequence[i]);
    }
  }
  if (tmp.size() % 4 != 0) {
    tmp.push_back(the_query.alignment_to_sequence[length() - 1]);
    tmp.push_back(the_template.alignment_to_sequence[length() - 1]);
  }
  for (core::index2 i = 0; i < tmp.size(); i += 4)
    out.push_back(AlignmentBlock(tmp[i], tmp[i + 2], tmp[i + 1], tmp[i + 3]));

  return out;
}

int PairwiseAlignment::first_index(int index, const AlignmentRow &the_row, const AlignmentRow &the_other_row) {

  core::index2 i = 0;
  while ((the_row.position_in_sequence(i) == -1) && (the_row.position_in_sequence(i) == -1)) {
    ++i;
    ++index;
  }
  return index;
}

} // ~alignment
} // ~core
