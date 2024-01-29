#include <memory>
#include <algorithm>

#include <core/index.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/structural/Chain.fwd.hh>
#include <core/data/structural/Chain.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/PdbAtom.hh>
#include <core/data/structural/structure_selectors.hh>
#include <core/data/sequence/SecondaryStructure.hh>

namespace core {
namespace data {
namespace structural {

Chain_SP Chain::clone(const ResidueSelector &which_residues)  const {

  Chain_SP ret = std::make_shared<Chain>(id_);

  for (const Residue_SP &ri : (*this))
    if (which_residues(*ri))
      ret->push_back(ri->clone(which_residues));

  return ret;
}

void Chain::push_back(Residue_SP r) {

  if (std::vector<Residue_SP>::size() != 0) {
    Residue_SP b = std::vector<Residue_SP>::back();
    if ((b->residue_type().type == 'P') && (r->residue_type().type == 'P')) {
      b->next(r);
      r->previous(b);
    }
    if ((b->residue_type().type == 'N') && (r->residue_type().type == 'N')) {
      b->next(r);
      r->previous(b);
    }
  }
  r->owner_ = shared_from_this();
  std::vector<Residue_SP>::push_back(r);
}

core::index2 Chain::count_aa_residues() const {

  return std::count_if(begin(),end(),[](const Residue_SP r){ return r->residue_type().type=='P'; });
}

core::index2 Chain::count_na_residues() const {

  return std::count_if(begin(),end(),[](const Residue_SP r){ return r->residue_type().type=='N'; });
}

Chain_SP Chain::create_ca_chain(const std::string & sequence, const char chain_id) {

  real dx = 3.67;
  real dy = 0.5;
  real sig = -1.0;
  Chain_SP c = std::make_shared<Chain>(chain_id);
  for (core::index2 i = 0; i < sequence.size(); ++i) {
    PdbAtom_SP a = std::make_shared<PdbAtom>(i + 1, std::string(" CA "));
    a->x = i * dx;
    a->y = (i + 1) * dy * sig;
    sig *= -1;
    a->z = 0.0;
    Residue_SP r = std::make_shared<Residue>(i + 1, core::chemical::Monomer::get(sequence[i]));
    r->push_back(a);
    c->push_back(r);
  }

  return c;
}

size_t chain_to_coordinates(const Chain_SP structure, core::data::basic::Coordinates & coordinates) {

  if (coordinates.size() != structure->size()) coordinates.resize(structure->size());
  int i = -1;
  for(auto a_it=structure->first_atom();a_it!=structure->last_atom();++a_it)
    coordinates[++i].set(*(*a_it));

  return (++i);
}

size_t chain_to_coordinates(const Chain_SP chain, core::data::basic::Coordinates & coordinates,
    AtomSelector_SP op) {

  if (op == nullptr) chain_to_coordinates(chain, coordinates);

  int i = -1;
  for (auto it = chain->first_atom(); it != chain->last_atom(); ++it)
    if ((*op)(**it)) coordinates[++i].set(**it);

  return (++i);
}

Residue_SP Chain::get_residue(const core::index2 residue_index) {

  for(auto & ri : *this)
    if(ri->id() == residue_index) return ri;
  return nullptr;
}

core::data::sequence::SecondaryStructure_SP Chain::create_sequence() const {

  std::string s;
  std::string ss;
  s.reserve(count_residues());
  ss.reserve(count_residues());
  for (auto r_it = begin(); r_it != end(); ++r_it) {
    if ((*r_it)->residue_type().type == 'P' || (*r_it)->residue_type().type == 'N') {
      s += (*r_it)->residue_type().code1;
      ss += (*r_it)->ss();
    }
  }
  return std::make_shared<core::data::sequence::SecondaryStructure>(owner()->code()+" "+id_,s,(*begin())->id(),ss);
}

void Chain::sort() {

  for(Residue_SP & r : *this)
    r->sort();
  std::sort(this->begin(),this->end());
}

std::ostream& operator<<(std::ostream &out, const Chain &c) {

  out << c.id();
  return  out;
}


}
}
}

