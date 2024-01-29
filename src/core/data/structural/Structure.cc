#include <memory>
#include <stdexcept>
#include <algorithm>

#include <core/index.hh>

#include <core/data/io/Pdb.hh>

#include <core/data/basic/Vec3Cubic.hh>

#include <core/data/structural/Chain.hh>
#include <core/data/structural/Structure.hh>
#include <core/data/structural/structure_selectors.hh>

namespace core {
namespace data {
namespace structural {

Structure_SP Structure::clone(const ChainSelector &which_chains) const {

  Structure_SP ret = std::make_shared<Structure>(code_);

  for (const Chain_SP &ci : (*this))
    if (which_chains(*ci))
      ret->push_back(ci->clone(which_chains));

  return ret;
}

void Structure::push_back(std::shared_ptr<Chain> c) { c->owner_ = shared_from_this(); std::vector<Chain_SP>::push_back(c);}

void Structure::sort() {

  std::sort(begin(),end(),[](const Chain_SP & ai,const Chain_SP & aj){ return ai->id() < aj->id();});
  for(Chain_SP c : *this) c->sort();
}

bool Structure::has_chain(const char code) const {
  for (const std::shared_ptr<Chain> & c : *this)
    if (c->id() == code) return true;

  return false;
}

const std::shared_ptr<Chain> Structure::get_chain(const char code) const {

  for(const std::shared_ptr<Chain> & c :*this) {
    if(c->id()==code) return c;
  }
  std::string codes = "";
  for(const Chain_SP & c:*this) codes +=c->id();
  std::string msg = std::string("Invalid chain code: '") + code + std::string("'; Registered chain codes: " + codes) + std::string("\n");
  throw std::invalid_argument(msg);
}

std::shared_ptr<Chain> Structure::get_chain(const char code) {

  for(std::shared_ptr<Chain> & c :*this) {
    if(c->id()==code) return c;
  }
  std::string codes = "";
  for(const Chain_SP & c:*this)
    codes +=c->id();
  std::string msg = std::string("Invalid chain code: '") + code + std::string("'; Registered chain codes: " + codes)
      + std::string("\n");
  logger << utils::LogLevel::CRITICAL << msg;
  throw std::invalid_argument(msg);
}

std::shared_ptr<Residue> Structure::get_residue(const char chain_code, const core::index2 residue_id, const char icode) {

  for(auto & r : *get_chain(chain_code))
    if((r->id()==residue_id)&&(r->icode()==icode)) return r;
  throw std::invalid_argument("Invalid residue request");
}

const std::shared_ptr<Residue> Structure::get_residue(const char chain_code, const core::index2 residue_id, const char icode) const {

  for(const auto & r : *get_chain(chain_code))
    if((r->id()==residue_id)&&(r->icode()==icode)) return r;
  throw std::invalid_argument("Invalid residue request");
}

std::shared_ptr<Residue> Structure::get_residue(const core::index2 residue_index) {

  core::index2 total_res = 0;
  short ichain=-1;
  do {
    total_res = operator[](++ichain)->count_residues();
  } while(total_res<residue_index);

  return operator[](ichain)->operator[](residue_index - total_res + operator[](ichain)->count_residues());
}

const std::shared_ptr<Residue> Structure::get_residue(const core::index2 residue_index) const {

  core::index2 total_res = 0;
  short ichain=-1;
  do {
    total_res = operator[](++ichain)->count_residues();
  } while(total_res<residue_index);

  return operator[](ichain)->operator[](residue_index - total_res + operator[](ichain)->count_residues());
}

std::vector<core::chemical::Monomer> & Structure::original_sequence(const char chain_code) {

  static const std::string key("SEQRES");
  auto it = pdb_header.find(key);
  if (it != pdb_header.end()) {
    const std::shared_ptr<io::Seqres> & seqr_sp = std::dynamic_pointer_cast<io::Seqres>((it->second));
    return seqr_sp->sequences[chain_code];
  } else throw std::runtime_error("Input PDB file has no SEQRES data");
}

size_t coordinates_to_structure(const core::data::basic::Coordinates & coordinates, Structure & structure) {

  int i = -1;
  std::for_each(structure.first_atom(), structure.last_atom(),
    [&](core::data::structural::PdbAtom_SP a) {(*a).set(coordinates[++i]);});

  return (++i);
}

size_t structure_to_coordinates(const Structure_SP structure, core::data::basic::Coordinates & coordinates) {

  if (coordinates.size() != structure->count_atoms()) coordinates.resize(structure->count_atoms());
  int i = -1;
  std::for_each(structure->first_atom(), structure->last_atom(),
    [&](core::data::structural::PdbAtom_SP a) {coordinates[++i].set(*a);});

  return (++i);
}

template <typename C>
size_t structure_to_coordinates(const Structure_SP structure, std::unique_ptr<C[]>  & coordinates) {

  int i = -1;
  std::for_each(structure->first_atom(), structure->last_atom(),
    [&](core::data::structural::PdbAtom_SP a) {coordinates[++i].set(*a);});

  return (++i);
}

template
size_t structure_to_coordinates(const Structure_SP structure, std::unique_ptr<basic::Vec3[]>  & coordinates);
template
size_t structure_to_coordinates(const Structure_SP structure, std::unique_ptr<basic::Vec3Cubic[]>  & coordinates);


template <typename C>
size_t coordinates_to_structure(const std::unique_ptr<C[]> &coordinates, const Structure_SP structure) {

  int i = -1;
  std::for_each(structure->first_atom(), structure->last_atom(), [&](core::data::structural::PdbAtom_SP a) {a->set(coordinates[++i]);});

  return (++i);
}

template
size_t coordinates_to_structure(const std::unique_ptr<basic::Vec3[]> &coordinates, Structure_SP structure);

template
size_t coordinates_to_structure(const std::unique_ptr<basic::Vec3Cubic[]> &coordinates, Structure_SP structure);


size_t structure_to_coordinates(const Structure_SP structure, core::data::basic::Coordinates & coordinates,
    const AtomSelector_SP op) {

  if (op == nullptr) return structure_to_coordinates(structure, coordinates);

  return structure_to_coordinates(structure,coordinates,*op);
}

size_t structure_to_coordinates(const Structure_SP structure, core::data::basic::Coordinates & coordinates,
    const AtomSelector & op) {

  utils::Logger l("structure_to_coordinates");

  size_t cnt = 0;
  for (auto it = structure->first_atom(); it != structure->last_atom(); ++it)
    if ((op)(**it)) cnt++;

  if (cnt != coordinates.size()) {
    l << utils::LogLevel::INFO << "Atomic coordinates for structure " << structure->code() << " " << cnt
        << " atoms found, vector resized from " << coordinates.size() << "\n";
    coordinates.resize(cnt);
  } else l << utils::LogLevel::FINE << cnt << " atomic coordinates for structure " << structure->code() << "\n";

  int i = -1;
  for (auto it = structure->first_atom(); it != structure->last_atom(); ++it) {
    if ((op)(**it))
      coordinates[++i].set(**it);
  }

  return (++i);
}

}
}
}
