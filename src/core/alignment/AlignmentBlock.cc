#include <core/index.hh>
#include <utils/string_utils.hh>
#include <core/alignment/AlignmentBlock.hh>

namespace core {
namespace alignment {

std::ostream & operator<<(std::ostream & out, const AlignmentBlock & block) {

  out << utils::string_format("%4d %4d %4d %4d\n", block.first_query_pos, block.last_query_pos, block.first_tmplt_pos, block.last_tmplt_pos);

  return out;
}


std::ostream & operator<<(std::ostream & out, const AlignmentBlock_SP  block) {

  out << utils::string_format("%4d %4d %4d %4d\n", block->first_query_pos, block->last_query_pos, block->first_tmplt_pos, block->last_tmplt_pos);

  return out;
}

}
}
