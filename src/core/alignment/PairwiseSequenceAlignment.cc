#include <string>

#include <core/real.hh>

#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/sequence_utils.hh>
#include <core/alignment/PairwiseAlignment.hh>
#include <core/alignment/PairwiseSequenceAlignment.hh>

namespace core {
namespace alignment {

PairwiseSequenceAlignment::PairwiseSequenceAlignment(const std::string &  query_name, const std::string & aligned_query,
    const core::index2 first_query_residue, const std::string &  template_name, const std::string & aligned_template, const core::index2 first_template_residue,
    const core::real score, bool expand_unaligned) {

  alignment = std::make_shared<PairwiseAlignment>(aligned_query,first_query_residue,aligned_template,first_template_residue,score,expand_unaligned);
  std::string tmp_q = aligned_query;
  core::data::sequence::remove_gaps(tmp_q);
  query_sequence = std::make_shared<core::data::sequence::Sequence>(query_name, tmp_q,first_query_residue);
  std::string tmp_t = aligned_template;
  core::data::sequence::remove_gaps(tmp_t);
  template_sequence = std::make_shared<core::data::sequence::Sequence>(template_name, tmp_t,first_template_residue);
}

PairwiseSequenceAlignment::PairwiseSequenceAlignment(const core::data::sequence::Sequence_SP aligned_query, const core::data::sequence::Sequence_SP aligned_template, const core::real score) {

  alignment = std::make_shared<PairwiseAlignment>(*aligned_query,*aligned_template,score);
  query_sequence = aligned_query->create_ungapped_sequence();
  template_sequence = aligned_template->create_ungapped_sequence();
}

std::ostream & operator<<(std::ostream & out, const PairwiseSequenceAlignment & seq_ali) {

  out << "# score: " << seq_ali.alignment_score() << " length: " << seq_ali.alignment->length() << "\n"
      << ">"<<seq_ali.query_sequence->header()<<"\n"<<seq_ali.get_aligned_query()<<"\n"
      << ">"<<seq_ali.template_sequence->header()<<"\n"<<seq_ali.get_aligned_template()<<"\n";

  return out;
}


} // ~alignment
} // ~core
