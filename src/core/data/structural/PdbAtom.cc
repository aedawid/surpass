#include <string>

#include <utils/string_utils.hh> // for string_format()

#include <core/data/structural/Residue.hh>
#include <core/data/structural/Chain.hh>
#include <core/data/structural/PdbAtom.hh>

#include <core/chemical/Monomer.hh>

namespace core {
namespace data {
namespace structural {

std::string pdb_atom_name(const std::string &atom_name) {

  if (atom_name.size() == 1) return " " + atom_name + "  ";
  if (atom_name.size() == 2) return " " + atom_name + " ";
  if (atom_name.size() == 3) return " " + atom_name;
  return atom_name;
}

PdbAtom::PdbAtom(const core::index4 id, const std::string &atom_name, const core::index2 element_index) :
  id_(id), atom_name_(atom_name), element_index_(element_index) {
  alt_locator_ = ' ';
  occupancy_ = real(1.0);
  b_factor_ = real(99.99);
  is_heteroatom_ = false;
}

PdbAtom::PdbAtom(const core::index4 id, const std::string &atom_name,
                 const core::real cx, const core::real cy, const core::real cz, const core::real occupancy,
                 const core::real b_factor, const core::index2 element_index) :
  Vec3(cx, cy, cz), id_(id), atom_name_(atom_name), element_index_(element_index), alt_locator_(' '),
  occupancy_(occupancy), b_factor_(b_factor) {
  is_heteroatom_ = false;
}

PdbAtom_SP PdbAtom::clone() const {
  PdbAtom_SP ret = std::make_shared<PdbAtom>(id_,atom_name_,x,y,z,occupancy_,b_factor_,element_index_);
  ret->is_heteroatom_ = is_heteroatom_;
  return ret;
}

std::string PdbAtom::to_pdb_line() const {

  if (Residue_SP rr = owner_.lock()) {
    const Residue & r = *rr;
    return utils::string_format(core::data::io::Atom::atom_format_uncharged, id_, atom_name_.c_str(), ' ',
        r.residue_type().code3.c_str(), r.owner()->id(), r.id(), r.icode(), x, y, z, occupancy_, b_factor_,
        core::chemical::AtomicElement::periodic_table[element_index_].symbol.c_str());
  } else {
    return utils::string_format(core::data::io::Atom::atom_format_uncharged, id_, atom_name_.c_str(), ' ',
        "ALA", 'A', 0, ' ', x, y, z, occupancy_, b_factor_,
        core::chemical::AtomicElement::periodic_table[element_index_].symbol.c_str());
  }
}

}
}
}
