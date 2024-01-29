#include <string>
#include <limits>
#include <iostream>
#include <cctype>
#include <memory>

#include <core/index.hh>
#include <core/chemical/Monomer.hh>
#include <core/data/structural/PdbAtom.hh>
#include <core/data/structural/structure_selectors.hh>

#include <utils/string_utils.hh>

namespace core {
namespace data {
namespace structural {

utils::Logger logger("structure_selectors");

//const std::string AtomSelector::selection_string_ = "";
const std::string SelectAllAtoms::selection_string_ = "*";
const std::string IsCA::selection_string_ = " CA ";
const std::string IsCB::selection_string_ = " CB ";
//const std::string IsBB::selection_string_ = "_N__+_CA_+_O__+_C__+_H__+_HA_+_HA1+_HA2+_HA3";
//const std::string IsBBCB::selection_string_ = "_N__+_CA_+_O__+_C__+_CB_+_H__+_HA_+_HA1+_HA2+_HA3";

bool IsNamedAtom::operator()(const PdbAtom & a) const {

  logger << utils::LogLevel::FINEST << "Selecting atom " << a.atom_name() << " : "
      << ((atom_name_[0] == '*') || (a.atom_name().compare(atom_name_) == 0)) << "\n";
  return ((atom_name_[0] == '*') || (a.atom_name().compare(atom_name_) == 0));
}

void IsNamedAtom::set(const std::string & atom_name) {

  atom_name_ = atom_name;
  std::replace(atom_name_.begin(), atom_name_.end(), '_', ' ');
  padded_atom_name_ = atom_name;
  logger << utils::LogLevel::FINE << "IsNamedAtom selector set to " << atom_name_ << "\n";
}

bool IsElement::operator()(const PdbAtom & a) const {

  return (a.element_index()==element_index_);
}

void IsElement::set(const std::string & element_name) {

  element_index_ = core::chemical::AtomicElement::by_symbol(element_name).z;
  logger << utils::LogLevel::FINE << "IsElement selector set to " << element_name << "\n";
}

bool ProperlyConnectedCA::operator()(const Residue &r) const {

  PdbAtom_SP my_ca = get_ca(r);
  if (my_ca == nullptr) return false;

  if (r.previous() != nullptr) {
    PdbAtom_SP prev_ca = get_ca(*r.previous());
    if (prev_ca == nullptr) return false;
    if (prev_ca->distance_to(*my_ca) > cutoff_) return false;
  }
  if (r.next() != nullptr) {
    PdbAtom_SP next_ca = get_ca(*r.next());
    if (next_ca == nullptr) return false;
    if (next_ca->distance_to(*my_ca) > cutoff_) return false;
  }

  return true;
}

PdbAtom_SP ProperlyConnectedCA::get_ca(const Residue & r) const {

  for(PdbAtom_SP ai : r)
    if (ai->atom_name().compare(" CA ") == 0) return ai;

  return nullptr;

}

bool ResidueIsTerminal::operator()(const Residue &r) const {

  if (r.next() == nullptr) return true;
  if (r.next()->residue_type().type != r.residue_type().type) return true;
  if (r.previous() == nullptr) return true;
  if (r.previous()->residue_type().type != r.residue_type().type) return true;

  return false;
}

bool IsFirstResidue::operator()(const Residue &r) const {

  if (r.previous() == nullptr) return true;
  if (r.previous()->residue_type().type != r.residue_type().type) return true;

  return false;
}

bool IsLastResidue::operator()(const Residue &r) const {

  if (r.next() == nullptr) return true;
  if (r.next()->residue_type().type != r.residue_type().type) return true;

  return false;
}

bool ResidueHasCA::operator ()(const Residue & r) const {

  for(PdbAtom_SP ai : r)
    if (ai->atom_name().compare(" CA ") == 0) return true;

  return false;
}

bool ResidueHasBB::operator ()(const Residue & r) const {

  if(r.size()==0) return false;
  index1 n = 0;
  for(PdbAtom_SP ai : r) {
    if (ai->atom_name().compare(" N  ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" CA ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" C  ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" O  ") == 0) { ++n; continue; }
  }

  if(n==4) return true;
  return false;
}

bool ResidueHasBBCB::operator ()(const Residue & r) const {

  if(r.size()==0) return false;
  index1 n = 0;
  for(PdbAtom_SP ai : r) {
    if (ai->atom_name().compare(" N  ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" CA ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" C  ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" O  ") == 0) { ++n; continue; }
    if (ai->atom_name().compare(" CB ") == 0) { ++n; continue; }
  }

  if ((n == 4) && (r.residue_type() == core::chemical::Monomer::GLY)) return true;
  if (n == 5) return true;
  return false;
}

core::index1 ResidueHasAllHeavyAtoms::expected_heavy_atoms(const Residue &r) const {

  if (r.residue_type().type == 'P')
    return (r.next() == nullptr) ? r.residue_type().n_heavy_atoms : r.residue_type().n_heavy_atoms - 1;

  if (r.residue_type().type == 'N')
    return (r.previous() == nullptr) ? r.residue_type().n_heavy_atoms : r.residue_type().n_heavy_atoms - 1;

  return r.residue_type().n_heavy_atoms;
}

bool ResidueHasAllHeavyAtoms::operator ()(const Residue & r) const {

  return (r.count_heavy_atoms() == expected_heavy_atoms(r));
}

bool IsAA::operator ()(const Residue & r) const {

  if(r.size()==0) return false;
  return r.residue_type().type == 'P';
}

bool IsAromaticAA::operator()(const Residue &r) const {

  using core::chemical::Monomer;
  if (IsAA::operator()(r))
    return (r.residue_type() == Monomer::TRP) || (r.residue_type() == Monomer::PHE) ||
           (r.residue_type() == Monomer::TYR) || (r.residue_type() == Monomer::HIS);
  return false;
}

bool IsNT::operator ()(const Residue & r) const {

  if(r.size()==0) return false;
  return r.residue_type().type == 'N';
}

bool ResidueSelector::operator()(const PdbAtom & a) const {

  if(a.owner()== nullptr) return false;
  return operator()(a.owner());
}

SelectResidueByName::SelectResidueByName(const std::string & selected_code3) {
  matching_code3.push_back(selected_code3);
  logger << utils::LogLevel::FINE << "selecting "<<selected_code3<<" residues\n";
}

SelectResidueByName::SelectResidueByName(const std::vector<std::string> & selected_code3) {

  std::string s;
  for(const std::string & n : selected_code3) {
    matching_code3.push_back(n);
    s += n + " ";
  }
  logger << utils::LogLevel::FINE << "selecting residues: "<<s<<"\n";
}

bool SelectResidueByName::operator()(const Residue &r) const {

  return (std::find(matching_code3.cbegin(), matching_code3.cend(), r.residue_type().code3) != matching_code3.cend());
}

ResidueSelector_SP residue_selector_from_string(const std::string & selection) {

  if (selection == "*") return std::make_shared<ResidueSelector>(); // --- '*' selects all residues
  if ((selection == "aa") || (selection == "AA")) return std::make_shared<IsAA>(); // --- 'aa' selects amino acids
  if ((selection == "nt") || (selection == "NT")) return std::make_shared<IsNT>(); // --- 'nt' selects nucleic acids

  return std::make_shared<SelectResidueRange>(selection);
}

bool ChainSelector::operator()(const Chain & c) const {

  logger << utils::LogLevel::FINER << "Selecting chain " << c.id() << " : "
      << ((chain_id_ == '*') || (c.id() == chain_id_)) << "\n";
  return ((chain_id_ == '*') || (c.id() == chain_id_));
}

bool SelectResidueRange::operator()(const Residue & r) const {

  if (first_residue > last_residue) return true; // --- This trick is used for selecting everything, i.e. to use '*' as a selector
  if ((r.id() > first_residue) && (r.id() < last_residue)) return true;
  if ((r.id() == first_residue) && (r.icode() >= first_residue_icode)) return true;
  if ((r.id() == last_residue) && (r.icode() <= last_residue_icode)) return true;

  return false;
}

void SelectResidueRange::set(const std::string & selector) {

  if ((selector.size() == 0)||(selector[0] == '*')) {
    first_residue = 0;
    last_residue = -1;
    return;
  }
  std::string s(selector);
  utils::trim(s);
  bool neg = false;
  if (s[0] == '-') {
    neg = true;
    s.erase(0, 1);
  }
  std::vector<std::string> tokens = utils::split(s, '-');

  if (isalpha(((std::string&) tokens[0]).back())) {
    first_residue_icode = ((std::string&) tokens[0]).back();
    ((std::string&) tokens[0]).pop_back();
  } else first_residue_icode = ' ';
  first_residue = utils::from_string<int>((std::string&) tokens[0]);
  if (neg) first_residue = -first_residue;
  if (tokens.size() == 1) {
    last_residue = first_residue;
    last_residue_icode = first_residue_icode;
  } else {
    if (isalpha(((std::string&) tokens[1]).back())) {
      last_residue_icode = ((std::string&) tokens[1]).back();
      ((std::string&) tokens[1]).pop_back();
    } else last_residue_icode = ' ';
    last_residue = utils::from_string<int>((std::string&) tokens[1]);
  }
  update_selector_string();
}

void SelectChainResidues::set(const std::string & selector) {

  if (selector.size() > 0) {
    std::vector<std::string> tokens = utils::split(selector, ':');
    logger << utils::LogLevel::FINE << "SelectChainResidueAtom selector set to >" << selector << "< with "
        << tokens.size() << " part(s)\n";
    if (chain_selector == nullptr) chain_selector = std::make_shared<ChainSelector>(tokens[0][0]);
    else chain_selector->set(tokens[0][0]);
    if (residue_selector == nullptr) {
      if (tokens.size() >= 2) residue_selector = std::make_shared<SelectResidueRange>(tokens[1]);
      else residue_selector = std::make_shared<SelectResidueRange>("*");
    }
  } else logger << utils::LogLevel::WARNING << "Empty string used to set SelectChainResidueAtom selector\n";

  selector_ = selector;
}



void ChainSelector::set(const std::string & chain_id) {

  if (chain_id.size() == 0) {
    this->chain_id_ = '*';
    selector = "*";
  } else {
    chain_id_ = chain_id[0];
    selector = chain_id;
  }
}

void SelectChainResidueAtom::set(const std::string & selector) {

  if (selector.size() > 0) {
    // ---------- Split the combined selector into parts
    std::vector<std::string> tokens = utils::split(selector, ':');
    logger << utils::LogLevel::FINE << "SelectChainResidueAtom selector set to >" << selector << "< with "
        << tokens.size() << " part(s)\n";
    // ---------- Create or set the chain selection part
    if (chain_selector == nullptr) chain_selector = std::make_shared<ChainSelector>(tokens[0][0]);
    else chain_selector->set(tokens[0][0]);
    // ---------- Create or set the residue selection part
    if (residue_selector == nullptr) {
      if (tokens.size() >= 2) residue_selector = residue_selector_from_string(tokens[1]);
      else residue_selector = std::make_shared<ResidueSelector>();
    }
    // ---------- Create or set the atom selection part
    if (atom_selector == nullptr) atom_selector = std::make_shared<IsNamedAtom>("*");
    if (tokens.size() == 3) {
      // ---------- Is it a simple named atom selector?
      if (tokens[2].find('+') == std::string::npos) atom_selector->set(tokens[2]);
      else {
        std::shared_ptr<CompositeSelector<IsNamedAtom>> cp = std::make_shared<CompositeSelector<IsNamedAtom>>('+');
        std::vector<std::string> subtokens = utils::split(tokens[2],'+');
        for(const std::string & st:subtokens)
          cp->add(st);
        atom_selector = cp;
      }
    }
    this->selector = chain_selector->selector_string() + ":" + residue_selector->selector_string() + ":"
        + atom_selector->selector_string();
    logger << utils::LogLevel::INFO << "SelectChainResidueAtom selector set to " << this->selector << "\n";
  } else logger << utils::LogLevel::WARNING << "Empty string used to set SelectChainResidueAtom selector\n";
}

std::ostream & operator<<(std::ostream &out, const ChainSelector &selector) {

  out << selector.chain_id_;
  return out;
}

std::ostream & operator<<(std::ostream &out, const AtomSelector &selector) {
  out << selector.selector_string();
  return out;
}

std::ostream & operator<<(std::ostream &out, const SelectResidueRange &selector) {

  if (selector.first_residue > selector.last_residue) {
    out << "*";
    return out;
  }

  out << selector.first_residue;
  if (selector.first_residue_icode != ' ') out << selector.first_residue_icode;
  out << "-" << selector.last_residue;
  if (selector.last_residue_icode != ' ') out << selector.last_residue_icode;
  return out;
}

std::ostream & operator<<(std::ostream &out, const SelectChainResidues &selector) {

  out << selector.chain_selector << ':' << selector.residue_selector;
  return out;
}

std::ostream & operator<<(std::ostream &out, const SelectChainResidueAtom &selector) {

  out << *(selector.chain_selector) << ':' << *(selector.residue_selector) << ':'
      << selector.atom_selector->selector_string();
  return out;
}

bool LogicalANDSelector::operator()(const PdbAtom &c) const {

  for (AtomSelector_SP sel : selectors) if (!(*sel)(c)) return false;
  return true;
}

}
}
}
