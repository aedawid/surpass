#include <algorithm>
#include <set>
#include <iomanip>

#include <core/data/structural/Structure.hh>

#include <simulations/systems/AtomTypingBase.hh>

namespace simulations {
namespace systems {

utils::Logger AtomTypingBase::l("AtomTypingBase");

AtomTypingBase::AtomTypingBase(const std::vector<AtomNamePair> & atom_types) {

  int i = -1;
  for (auto & e : atom_types) {
    index_to_internal_name.push_back(e.internal_name);
    internal_name_to_index[e.internal_name] = ++i;
  }
}

AtomTypingBase::AtomTypingBase(const std::vector<std::string> & atom_names) {

  int i = -1;
  for (const std::string & e : atom_names) {
    index_to_internal_name.push_back(e);
    internal_name_to_index[e] = ++i;
  }
}

void AtomTypingBase::show_known_atom_types(std::ostream &output) const {

  std::vector<std::string> tmp;
  for (auto k : internal_name_to_index) tmp.push_back(k.first);
  std::sort(tmp.begin(), tmp.end());

  for (auto & k : tmp) {
    utils::replace_substring(k," ","_");
    output << std::setw(10)  << k << "\n";
  }
}

void AtomTypingBase::compile_types(const core::data::structural::Structure & structure) {

  std::set<AtomNamePair> atom_names;

  for(const auto cit : structure)
    for(const auto rit : *cit)
      for(const auto ait : *rit) {
        atom_names.insert(AtomNamePair(ait->atom_name(),core::chemical::AtomicElement::periodic_table[ait->element_index()].symbol));
      }
  int i = -1;
  for (auto & e : atom_names) {
    index_to_internal_name.push_back(e.internal_name);
    internal_name_to_index[e.internal_name] = ++i;
  }
}

}
}
