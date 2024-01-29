#include <utils/io_utils.hh>
#include <utils/string_utils.hh>
#include <core/data/io/Pdb.hh>
#include <core/data/structural/structure_selectors.hh>

#include <core/data/structural/Structure.hh>
//#include <core/data/structural/Chain.hh>
//#include <core/data/structural/Residue.hh>
//#include <core/data/structural/PdbAtom.hh>
#include <core/chemical/Monomer.hh>
#include <utils/Logger.hh>

#include <simulations/representations/surpass_utils.hh>

namespace simulations {
namespace representations {

utils::Logger surpass_utils_logger("surpass_utils");

core::data::structural::Structure &fix_surpass_ss_assignment(core::data::structural::Structure &surpass_strctr) {

  for (auto at_it = surpass_strctr.first_atom(); at_it != surpass_strctr.last_atom(); ++at_it) {
    if ((**at_it).atom_name() == " S  ") {
      (**at_it).owner()->ss('E');
      continue;
    }
    if ((**at_it).atom_name() == " H  ") {
      (**at_it).owner()->ss('H');
      continue;
    }
    (**at_it).owner()->ss('C');
  }
  return surpass_strctr;
}

core::data::structural::Structure_SP surpass_representation(const core::data::structural::Structure &strctr) {

  using namespace core::data::structural;
  // Parameters
  PdbAtom_SP next_CA_sp, next_SG;
  Residue_SP residue_sp, resid_sp, next_res_sp;
  Chain_SP ch_sp;
  Structure_SP structure = std::make_shared<Structure>(strctr.code());
  int resid_id, resid_number = 0;
  core::index2 n_res = 0;
  core::real BF = 0.0;
  char it_ss, ss_1, ss_2, ss_3, chain_id = 'o';
  bool is_OK = true;

// Iterate over all chains
  for (auto it_chain = strctr.begin(); it_chain != strctr.end(); ++it_chain) {
    Chain_SP chain_sp = *it_chain;
    if ((chain_sp->id() == chain_id) && (is_OK == false)) break;
    core::index2 chain_size = chain_sp->count_aa_residues();
    ch_sp = std::make_shared<core::data::structural::Chain>(chain_sp->id());
    structure->push_back(ch_sp);
    resid_number = (*((*it_chain)->begin()))->id() - 1;

// Iterate over all residues in the chain
    for (auto it_resid = (*it_chain)->begin(); it_resid != (*it_chain)->end(); ++it_resid) {
      PdbAtom_SP atom_CA_sp;
      ++n_res;
      if (n_res >= chain_size - 2) {
        chain_id = chain_sp->id();
        is_OK = false;
        break;
      }
      resid_sp = *it_resid;
      resid_id = resid_sp->id();
      if (resid_sp->residue_type().type == 'P') {
// Check the completeness of residues in the chain and atoms in the residue
        if (resid_number + 1 == resid_id) {
          if ((resid_number + 2 != (*chain_sp)[n_res]->id()) || (resid_number + 3 != (*chain_sp)[n_res + 1]->id()) ||
              (resid_number + 4 != (*chain_sp)[n_res + 2]->id())) {
            ++resid_number;
            continue;
          }
          ++resid_number;
          atom_CA_sp = resid_sp->find_atom(" CA ");
          if (atom_CA_sp == nullptr) {
            surpass_utils_logger << utils::LogLevel::WARNING << "Missing backbone atom CA from residue: " << resid_id
                                 << " " << resid_sp->residue_type().code3 << "\n";
            continue;
          }
        } else if (resid_number + 1 != resid_id) {
          surpass_utils_logger << utils::LogLevel::WARNING << "Missing all residue(s) from number " << resid_number + 1
                               << " to " << resid_id - 1 << "\n";
          resid_number = resid_id;
          atom_CA_sp = resid_sp->find_atom(" CA ");
          if ((resid_number + 1 != (*chain_sp)[n_res]->id()) || (resid_number + 2 != (*chain_sp)[n_res + 1]->id()) ||
              (resid_number + 3 != (*chain_sp)[n_res + 2]->id())) {
            ++resid_number;
            continue;
          }
          if (atom_CA_sp == nullptr) {
            surpass_utils_logger << utils::LogLevel::WARNING << "Missing backbone atom CA from residue: " << resid_id
                                 << " " << resid_sp->residue_type().code3 << "\n";
            continue;
          }
        }
// Make new representation (one ball for 4 subsequent residues; if sequence size is N, than you have N-3 balls)
        PdbAtom_SP atom_SG_sp = std::make_shared<core::data::structural::PdbAtom>(1, " SU ");
        atom_SG_sp->occupancy(1.00);
        atom_SG_sp->b_factor(0.00);
        atom_SG_sp->owner(resid_sp);
        *(atom_SG_sp) += *(atom_CA_sp);
        BF += atom_CA_sp->core::data::structural::PdbAtom::b_factor();
        it_ss = resid_sp->ss();
        ss_1 = (chain_sp->Chain::get_residue(resid_id + 1))->ss();
        ss_2 = (chain_sp->Chain::get_residue(resid_id + 2))->ss();
        ss_3 = (chain_sp->Chain::get_residue(resid_id + 3))->ss();
        if (!((it_ss == ss_1) && (it_ss == ss_2) && (it_ss == ss_3))) {
          if ((it_ss == 'C') && (ss_1 == ss_2) && (ss_1 == ss_3)) it_ss = ss_1;
          else if ((ss_3 == 'C') && (ss_1 == ss_2) && (ss_1 == it_ss)) it_ss = ss_1;
          else it_ss = 'C';
        }
        for (core::index2 i = 1; i < 4; ++i) {
          next_res_sp = chain_sp->get_residue(resid_id + i);
          next_CA_sp = next_res_sp->find_atom(" CA ");
          *(atom_SG_sp) += *(next_CA_sp);
          BF += next_CA_sp->core::data::structural::PdbAtom::b_factor();
        }
        *(atom_SG_sp) /= 4;
        if (it_ss == 'H') { atom_SG_sp->atom_name(" H  "); }
        else if (it_ss == 'E') { atom_SG_sp->atom_name(" S  "); }
        else if (it_ss == 'C') { atom_SG_sp->atom_name(" C  "); }
        atom_SG_sp->b_factor(BF / 4);
        residue_sp = std::make_shared<core::data::structural::Residue>(n_res, core::chemical::Monomer::GLY);
        atom_SG_sp->id(n_res);
        (*residue_sp).push_back(atom_SG_sp);
        residue_sp->ss(it_ss);
        (*ch_sp).push_back(residue_sp);
      }
    }
  }
  return structure;
}

core::data::sequence::SecondaryStructure_SP surpass_representation(
  const core::data::sequence::SecondaryStructure &secstr,
  const bool check_if_converted) {

  if (check_if_converted) {
    for (char c : secstr.sequence)
      if (c != 'G') {
        std::make_shared<core::data::sequence::SecondaryStructure>(secstr.header(), secstr.sequence, secstr.first_pos(),
          secstr.str());
      }
  }

  using namespace core::data::sequence;
  std::string seq(secstr.length() - 3, 'G');
  std::string sec(secstr.length() - 3, 'C');
  core::data::sequence::SecondaryStructure_SP out = std::make_shared<SecondaryStructure>(secstr.header(), seq,
    secstr.first_pos(), sec);
  std::string tmp = secstr.str();
  for (core::index2 i = 0; i < secstr.length() - 3; ++i) {
    if (tmp.compare(i, 4, "HHHH") == 0) {
      out->fractions(i, 1.0, 0, 0);
      continue;
    }
    if (tmp.compare(i, 4, "EEEE") == 0) {
      out->fractions(i, 0, 1.0, 0);
      continue;
    }
    if ((tmp.compare(i, 4, "CHHH") == 0) || (tmp.compare(i, 4, "HHHC") == 0)) {
      out->fractions(i, 0.75, 0, 0.25);
      continue;
    }
    if ((tmp.compare(i, 4, "EHHH") == 0) || (tmp.compare(i, 4, "HHHE") == 0)) {
      out->fractions(i, 0.75, 0.25, 0);
      continue;
    }
    if ((tmp.compare(i, 4, "CEEE") == 0) || (tmp.compare(i, 4, "EEEC") == 0)) {
      out->fractions(i, 0.0, 0.75, 0.25);
      continue;
    }
    if ((tmp.compare(i, 4, "HEEE") == 0) || (tmp.compare(i, 4, "EEEH") == 0)) {
      out->fractions(i, 0.25, 0.75, 0);
      continue;
    }
  }

  return out;
}

bool is_surpass_model(const core::data::structural::Structure &strctr) {

  using namespace core::data::structural;

  bool is_OK = true;
  for (auto it_chain = strctr.begin(); it_chain != strctr.end(); ++it_chain) {
    for (auto it_resid = (*it_chain)->begin(); it_resid != (*it_chain)->end(); ++it_resid) {
      if ((*it_resid)->residue_type().type == 'P') {
        if ((*it_resid)->count_atoms() != 1) {
          is_OK = false;
        }
        if (((**it_resid)[0]->atom_name() != " H  ") & ((**it_resid)[0]->atom_name() != " S  ") &
            ((**it_resid)[0]->atom_name() != " C  ")) {
          is_OK = false;
        }
      }
      if (is_OK == false) {
        surpass_utils_logger << utils::LogLevel::WARNING << "Chain  " << (*it_chain)->id()
                             << " is not in SURPASS representation!\n";
      }
      break;
    }

  }
  return is_OK;
}

}
}


