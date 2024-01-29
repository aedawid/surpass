#include <iostream>
#include <string>
#include <memory>

#include <core/data/structural/PdbAtom.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/structure_selectors.hh>

#include <utils/string_utils.hh> // for string_format()

namespace core {
namespace data {
namespace structural {

bool operator<(const Residue & ri, const Residue & rj) {

  if (ri.id() != rj.id()) return ri.id() < rj.id();
  if (ri.id() == rj.id()) return ri.icode() < rj.icode();
  return false;
}

bool operator<(const Residue_SP ri, const Residue_SP rj) {

  if (ri->id() != rj->id()) return ri->id() < rj->id();
  if (ri->id() == rj->id()) return ri->icode() < rj->icode();
  return false;
}

Residue_SP Residue::clone(const AtomSelector &which_atoms) const {

  Residue_SP ret = std::make_shared<Residue>(id_,residue_type_);
  ret->insertion_code_ = insertion_code_;
  ret->ss_type_ = ss_type_;

  for (const PdbAtom_SP &ai : *(this))
    if (which_atoms(*ai))
      ret->push_back(ai->clone());

  return ret;
}

void Residue::push_back(PdbAtom_SP a) { a->owner_ = shared_from_this(); std::vector<PdbAtom_SP>::push_back(a); }

std::ostream& operator<<(std::ostream &out, const Residue & r) {
	out << utils::string_format("%c%s %4d", r.insertion_code_, r.residue_type_.code3.c_str(), r.id_);
	return out;
}

utils::Logger &operator <<(utils::Logger &logger, const Residue & r) {

  logger << utils::string_format("%c%s %4d", r.insertion_code_, r.residue_type_.code3.c_str(), r.id_);
  return logger;
}


core::index2 Residue::count_heavy_atoms() const {

  core::index2 nnH = 0;
  for (size_t i = 0; i < size(); ++i)
    if (operator[](i)->element_index() != 1) ++nnH;
  return nnH;
}


void Residue::sort() {
  std::sort(begin(),end(),[](const PdbAtom_SP & ai,const PdbAtom_SP & aj){ return ai->id() < aj->id();});
}

const core::real Residue::min_distance(const Residue & another_residue) const {

  core::real min_val = std::numeric_limits<core::real>::max();

  for (core::index2 i = 0; i < size(); ++i) {
    const core::data::basic::Vec3 & vi = *std::static_pointer_cast<core::data::basic::Vec3>(operator [](i));
    for (core::index2 j = 0; j < another_residue.size(); ++j) {
      const core::data::basic::Vec3 & vj = *std::static_pointer_cast<core::data::basic::Vec3>(another_residue[j]);
      core::real d = vi.x - vj.x;
      core::real d2 = d * d;
      if (d2 > min_val) continue;
      d = vi.y - vj.y;
      d2 += d * d;
      if (d2 > min_val) continue;
      d = vi.z - vj.z;
      d2 += d * d;
      min_val = d2;
    }
  }

  return sqrt(min_val);
}
PdbAtom_SP Residue::find_atom(const std::string & atom_name)  {

  for(auto it=begin();it!=end();++it)
    if((*it)->atom_name().compare(atom_name)==0)
      return *it;

  return nullptr;
}

const PdbAtom_SP Residue::find_atom(const std::string & atom_name)  const {

  for(auto it=begin();it!=end();++it)
    if((*it)->atom_name().compare(atom_name)==0)
      return *it;

  return nullptr;
}

}
}
}
