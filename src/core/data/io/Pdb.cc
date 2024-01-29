#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <fstream>
#include <vector>

#include <core/real.hh>
#include <core/data/io/Pdb.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/structural/PdbAtom.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/Chain.hh>
#include <core/data/structural/Structure.hh>
#include <core/chemical/AtomicElement.hh>

#include <utils/Logger.hh>
#include <utils/exit.hh>
#include <utils/string_utils.hh>
#include <utils/io_utils.hh>
#include <utils/exceptions/NoSuchFile.hh>

namespace core {
namespace data {
namespace io {

using core::real;
using namespace core::data::structural;

utils::Logger Pdb::logger = utils::Logger("Pdb");

const PdbLineFilter invert_filter(const PdbLineFilter f1) {

  auto f = [=](const std::string line) {

    return !f1(line);
  };
  return f;
}

const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2) {

  auto f = [=](const std::string line) {

    return f1(line)&&f2(line);
  };
  return f;
}

const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2, const PdbLineFilter f3) {

  auto f = [=](const std::string & line) {

    return f1(line)&&f2(line)&&f3(line);
  };
  return f;
}

const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2, const PdbLineFilter f3, const PdbLineFilter f4) {

  auto f = [=](const std::string line) {

    return f1(line)&&f2(line)&&f3(line)&&f4(line);
  };
  return f;
}

const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2, const PdbLineFilter f3, const PdbLineFilter f4, const PdbLineFilter f5) {

  auto f = [=](const std::string line) {

    return f1(line)&&f2(line)&&f3(line)&&f4(line)&&f5(line);
  };
  return f;
}


Pdb::Pdb(const std::string & fname, const PdbLineFilter & predicate,
         const bool if_parse_header, const bool first_model_only) {

  logger << utils::LogLevel::FINER << "Attempting " << fname << "\n";
  if (if_parse_header)
    logger << utils::LogLevel::FINER << "PDB file header will be parsed\n";

  if(!utils::if_file_exists(fname)) {
    std::string msg = "Can't find PDB file: " + fname + "!\n";
    logger<<utils::LogLevel::SEVERE<<msg;
    throw utils::exceptions::NoSuchFile(fname);
  }
  std::shared_ptr<std::istream> in;
  if (utils::has_suffix(fname, ".gz")) {
    std::string tmp;
    utils::load_binary_file(fname, tmp);
    if (tmp.size() == 0) EXIT("Zero bytes loaded from a gzip'ed file: " + fname + "\n");

    std::shared_ptr<std::stringstream> ssp = std::make_shared<std::stringstream>();
    utils::ungzip_string(tmp,*ssp);
    in = std::static_pointer_cast<std::istream>(ssp);
    logger << utils::LogLevel::FILE << "Reading PDB from a gzip'ed file: " << fname << "\n";
  } else {
    in = std::static_pointer_cast<std::istream>(std::make_shared<std::ifstream>(fname));
    if (!*in) logger << utils::LogLevel::SEVERE << "Can't open file: " << fname << "\n";
    logger << utils::LogLevel::FILE << "Reading PDB from: " << fname << "\n";
  }
  fname_ = fname;
  read_pdb(*in, predicate, if_parse_header, first_model_only);
}

void Pdb::read_pdb(std::istream & infile, const PdbLineFilter & predicate,
                   const bool if_parse_header, const bool first_model_only) {

  std::string line;
  std::vector<Atom> v;
  atoms.push_back(std::make_shared<std::vector<Atom>>());
  std::shared_ptr<std::vector<Atom>> current_model = atoms.back();
  auto start = std::chrono::high_resolution_clock::now();
  std::shared_ptr<Seqres> seqres_sp = std::make_shared<Seqres>();
  while (std::getline(infile, line)) {
    if (line.size() < 3) continue;   // Silently skip very short (or empty) lines

//    std::cerr<<line<<"\n";

    if (line.compare(0, 6, "ENDMDL") == 0) {
      atoms.emplace_back(std::make_shared<std::vector<Atom>>());
      current_model = atoms.back();
      if(first_model_only) break;
      continue;
    }
    if (if_parse_header) {
      if (line.compare(0, 6, "SEQRES") == 0) {
        seqres_sp->add_line(line);
      }

      if (line.compare(0, 6, "MODRES") == 0) {
        std::shared_ptr<PdbField> p = std::static_pointer_cast<PdbField>(std::make_shared<Modres>(line));
        header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("MODRES", p));
        continue;
      }

      if (line.compare(0, 6, "HEADER") == 0) {
        std::shared_ptr<PdbField> p = std::static_pointer_cast<PdbField>(std::make_shared<Header>(line));
        header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("HEADER", p));
        continue;
      }

      if (line.compare(0, 6, "KEYWDS") == 0) {
        std::shared_ptr<PdbField> p = std::static_pointer_cast<PdbField>(std::make_shared<Keywords>(line));
        header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("KEYWDS", p));
        continue;
      }

      if (line.compare(0, 6, "DBREF ") == 0) {
        std::shared_ptr<PdbField> p = std::static_pointer_cast<PdbField>(std::make_shared<DBRef>(line));
        header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("DBREF", p));
        continue;
      }

      if (line.compare(0, 6, "HELIX ") == 0) {
        std::shared_ptr<PdbField> p = std::static_pointer_cast<PdbField>(std::make_shared<HelixField>(line));
        header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("HELIX", p));
        continue;
      }

      if (line.compare(0, 6, "SHEET ") == 0) {
        std::shared_ptr<PdbField> p = std::static_pointer_cast<PdbField>(std::make_shared<SheetField>(line));
        header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("SHEET", p));
        continue;
      }
    }
    if (!predicate(line)) {
      logger << utils::LogLevel::FINE << "Skipping line due to predicate: " << line << "\n";
      continue;
    }
    if (line.compare(0, 4, "ATOM") == 0) {
      if (line.size() < 53) {
        logger << utils::LogLevel::WARNING << "Skipping incomplete ATOM line: " << line << "\n";
        continue;
      }
      current_model->emplace_back(line);
      continue;
    }
    if (line.compare(0, 6, "HETATM") == 0) {
      if (line.size() < 53) {
        logger << utils::LogLevel::WARNING << "Skipping incomplete HETATM line: " << line << "\n";
        continue;
      }
      current_model->emplace_back(line);
      continue;
    }
  }

  if(seqres_sp!=nullptr)
    header.insert(std::pair<std::string, std::shared_ptr<PdbField>>("SEQRES", seqres_sp));

  if (current_model->size() == 0) {
    if (atoms.size() == 1)
      logger << utils::LogLevel::WARNING << "The input file stream contains no valid ATOM lines !\n";
    atoms.pop_back();
  }
  //---------- show timer and stats
  auto end = std::chrono::high_resolution_clock::now();
  logger << utils::LogLevel::INFO << "Found " << atoms.size() << " models, " << atoms[0]->size() << " atoms each after "
      << (double) (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) << " [ms]\n";
}

Structure_SP Pdb::create_structure(const core::index2 which_model) {

  Structure_SP structure = std::make_shared<Structure>(pdb_code());
  if(structure->code().size()==0) structure->code(pdb_code_from_file_name(fname_));
  char last_chain_code = 0;
  Chain_SP last_chain;
  int last_resid_id = -65535;
  char last_icode = ' ';
  Residue_SP last_residue;
  for (Atom aline : *atoms[which_model]) {
// ---------- do we have a new chain?
    if (aline.chain != last_chain_code) {
      last_chain_code = aline.chain;
      if (!structure->has_chain(aline.chain)) {
        last_chain = std::make_shared<Chain>(last_chain_code);
        structure->push_back(last_chain);
        logger << utils::LogLevel::FINE << "Creating a new chain: " << aline.chain << "\n";
      } else last_chain = structure->get_chain(aline.chain);
    }
    // ---------- do we have a new residue?
    if ((aline.residue_id != last_resid_id) || (aline.i_code != last_icode)) {
      last_icode = aline.i_code;
      last_resid_id = aline.residue_id;
      last_residue = std::make_shared<Residue>(aline.residue_id, aline.res_name);
      last_residue->icode(last_icode);
      last_chain->push_back(last_residue );
      logger << utils::LogLevel::FINER << "Creating a new residue: " << aline.residue_id <<" "<< aline.res_name << "\n";
    }

    PdbAtom_SP a =
        ((aline.element[0] != ' ') || (aline.element[1] != ' ')) ?
            std::make_shared<PdbAtom>(aline.serial, aline.name, aline.x, aline.y, aline.z, aline.occupancy,
                aline.temp_factor, core::chemical::AtomicElement::by_symbol(utils::trim(aline.element)).z) :
            std::make_shared<PdbAtom>(aline.serial, aline.name, aline.x, aline.y, aline.z, aline.occupancy,
                aline.temp_factor, core::chemical::AtomicElement::DUMMY.z);
    a->alt_locator(aline.alt_loc);
    a->is_heteroatom(aline.is_heteroatom);
    last_residue->push_back(a);
    logger << utils::LogLevel::FINEST << "Creating a new atom: " << *a << "\n";
  }

  // ---------- copy header items ----------
  structure->pdb_header.insert(header.begin(), header.end());

  // ---------- annotate residues with secondary structure ----------
  auto range = header.equal_range("HELIX");
  for (Pdb::HeaderIterator it = range.first; it != range.second; ++it) {
    std::shared_ptr<HelixField> h = std::static_pointer_cast<HelixField>((*it).second);
    for (auto rit = structure->first_residue(); rit != structure->last_residue(); ++rit)
      if (h->is_my_residue(**rit)) (*rit)->ss('H');
  }

  range = header.equal_range("SHEET");
  for (Pdb::HeaderIterator it = range.first; it != range.second; ++it) {
    std::shared_ptr<SheetField> h = std::static_pointer_cast<SheetField>((*it).second);
    for (auto rit = structure->first_residue(); rit != structure->last_residue(); ++rit)
      if (h->is_my_residue(**rit)) (*rit)->ss('E');
  }

  return structure;
}

void Pdb::create_structures(std::vector<core::data::structural::Structure_SP> & structures) {

  for(core::index4 i=0;i<count_models();++i) structures.push_back( create_structure(i) );
}

std::string Pdb::pdb_code() const {

  auto it = header.find("HEADER");
  if (it == header.end()) return "";
  else {
    std::shared_ptr<Header> h = std::static_pointer_cast<Header>(it->second);
    return h->pdb_id;
  }
}

std::string Header::to_pdb_line() const {

  return "HEADER    " + classification + date + "   " + pdb_id;
}

Keywords::Keywords(const std::string &pdb_line) {

  utils::split(pdb_line.substr(10), keywords, ',');
}

std::string Keywords::to_pdb_line() const {

  std::string line = "KEYWDS   ";
  std::string full_line = "";
  if (keywords.size() < 1) return line;
  line += keywords[0];
  std::cout << "\n" << full_line << "\n";
  int cnt = 1;
  for (size_t i = 1; i < keywords.size(); i++) {
    const std::string &s = keywords[i];
    if (line.length() + s.length() >= 76) {
      cnt++;
      full_line += line;
      line = utils::string_format("\nKEYWORDS  %2d", cnt) + " " + s;
    } else line += ", " + s;
  }

  full_line += line;

  return full_line;
}

DBRef::DBRef(const std::string &pdb_line) {

  pdb_code = pdb_line.substr(7, 4);
  chain = pdb_line[12];
  residue_from = utils::from_string<int>(pdb_line, 14, 17, -1);
  insert_from = pdb_line[18];
  residue_to = utils::from_string<int>(pdb_line, 20, 23, -1);
  insert_to = pdb_line[24];
  db_name = pdb_line.substr(26, 6);
  db_accession = pdb_line.substr(33, 8);
  db_code = pdb_line.substr(42, 12);
  db_residue_from = utils::from_string<int>(pdb_line, 55, 59, -1);
  db_insert_from = pdb_line[60];
  db_residue_to = utils::from_string<int>(pdb_line, 62, 66, -1);
  db_insert_to = pdb_line[67];
}

std::string DBRef::to_pdb_line() const {

  return utils::string_format("DBREF   %s %c %4d%c %4d%c %6s %8s %12s %4d%c %4d%c", pdb_code.c_str(), chain,
      residue_from, insert_from, residue_to, insert_to, db_name.c_str(), db_accession.c_str(), db_code.c_str(),
      db_residue_from, db_insert_from, db_residue_to, db_insert_to);
}

TVect::TVect(const std::string &pdb_line) {

  serial = utils::from_string<index2>(pdb_line, 7, 10, 1);
  x = utils::from_string<double>(pdb_line, 10, 20, 0.0);
  y = utils::from_string<double>(pdb_line, 20, 30, 0.0);
  z = utils::from_string<double>(pdb_line, 30, 40, 0.0);
}

std::string TVect::to_pdb_line() const {

  return utils::string_format("TVECT %4d%10.3f%10.3f%10.3f", serial, x, y, z);
}

HelixField::HelixField(const std::string &pdb_line) {

  serial = utils::from_string<int>(pdb_line, 7, 3, -1);
  helix_id = pdb_line.substr(11, 3);
  residue_name_from = pdb_line.substr(15, 3);
  chain_from = pdb_line[19];
  residue_id_from = utils::from_string<int>(pdb_line, 21, 4, -1);
  insert_from = pdb_line[25];

  residue_name_to = pdb_line.substr(27, 3);
  residue_id_to = utils::from_string<int>(pdb_line,33, 3,-1);
  chain_to = pdb_line[31];
  insert_to = pdb_line[37];
  helix_class = utils::from_string<int>(pdb_line,39, 2,-1);
  length =  utils::from_string<int>(pdb_line,72, 5,-1);
  comment = pdb_line.substr(40, 30);
}

std::string HelixField::to_pdb_line() const {

  return utils::string_format("HELIX  %3d %s %s %c %4d%c %3s %c %4d%c%2d %s %4d", serial, helix_id.c_str(),
      residue_name_from.c_str(), chain_from, residue_id_from, insert_from, residue_name_to.c_str(), chain_to,
      residue_id_to, insert_to, helix_class, comment.c_str(), length);
}

bool HelixField::is_my_residue(const core::data::structural::Residue & r) {

  if ((r.owner()->id() != chain_from) && (r.owner()->id() != chain_to)) return false;
  if ((r.id() > residue_id_from) && (r.id() < residue_id_to)) return true;
  if ((r.id() == residue_id_from) && (r.icode() >= insert_from)) return true;
  if ((r.id() == residue_id_to) && (r.icode() <= insert_to)) return true;
  return false;
}

SheetField::SheetField(const std::string &pdb_line) {

  strand_id = utils::from_string<int>(pdb_line, 7, 3, -1);
  sheet_id = pdb_line.substr(11, 3);
  n_strands = utils::from_string<int>(pdb_line, 14, 2, -1);

  residue_name_from = pdb_line.substr(17, 3);
  chain_from = pdb_line[21];
  residue_id_from = utils::from_string<int>(pdb_line, 22, 4, -1);
  insert_from = pdb_line[26];

  residue_name_to = pdb_line.substr(28, 3);
  chain_to = pdb_line[32];
  residue_id_to = utils::from_string<int>(pdb_line, 33, 4, -1);
  insert_to = pdb_line[37];

  sense = utils::from_string<int>(pdb_line, 38, 2, -1);

  if (sense != 0) {
    register_my_atom = pdb_line.substr(41, 4);
    register_my_residue_name = pdb_line.substr(45, 3);
    register_my_chain = pdb_line[49];
    register_my_residue_id = utils::from_string<int>(pdb_line, 51, 4, -1);
    register_my_insert = pdb_line[54];

    register_previous_atom = pdb_line.substr(56, 4);
    register_previous_residue_name = pdb_line.substr(60, 3);
    register_previous_chain = pdb_line[64];
    register_previous_residue_id = utils::from_string<int>(pdb_line, 65, 4, -1);
    register_previous_insert = pdb_line[69];
  } else {
    register_my_atom = "    ";
    register_my_residue_name = "   ";
    register_my_chain = ' ';
    register_my_residue_id = -1;
    register_my_insert = ' ';

    register_previous_atom = "    ";
    register_previous_residue_name = "   ";
    register_previous_chain = ' ';
    register_previous_residue_id = -1;
    register_previous_insert = ' ';
  }
}

std::string SheetField::to_pdb_line() const {

  if (sense != 0) {
    return utils::string_format("SHEET  %3d %3s%2d %3s %c%4d%c %3s %c%4d%c%2d %4s%3s %c%4d%c %4s%3s %c%4d%c", strand_id,
        sheet_id.c_str(), n_strands, residue_name_from.c_str(), chain_from, residue_id_from, insert_from,
        residue_name_to.c_str(), chain_to, residue_id_to, insert_to, sense, register_my_atom.c_str(),
        register_my_residue_name.c_str(), register_my_chain, register_my_residue_id, register_my_insert, register_previous_atom.c_str(),
        register_previous_residue_name.c_str(), register_previous_chain, register_previous_residue_id,
        register_previous_insert);
  } else {
    return utils::string_format("SHEET  %3d %3s%2d %3s %c%4d%c %3s %c%4d%c%2d", strand_id, sheet_id.c_str(),
        n_strands, residue_name_from.c_str(), chain_from, residue_id_from, insert_from, residue_name_to.c_str(),
        chain_to, residue_id_to, insert_to, sense);
  }
}

bool SheetField::is_my_residue(const core::data::structural::Residue & r) {

  if ((r.owner()->id() != chain_from) && (r.owner()->id() != chain_to)) return false;
  if ((r.id() > residue_id_from) && (r.id() < residue_id_to)) return true;
  if ((r.id() == residue_id_from) && (r.icode() >= insert_from)) return true;
  if ((r.id() == residue_id_to) && (r.icode() <= insert_to)) return true;
  return false;
}

Modres::Modres(const std::string &pdb_line) {

  pdb_id = pdb_line.substr(7, 4);
  res_name = pdb_line.substr(12, 3);
  chain = pdb_line[16];
  residue_id = utils::from_string<int>(pdb_line, 18, 21, -1);
  i_code = pdb_line[23];
  std_res_name = pdb_line.substr(24, 3);
  comment = pdb_line.substr(29, 40);
}

void Seqres::add_line(const std::string & seqres_line) {

  char chain = seqres_line[11];
  core::index2 n = utils::from_string<core::index2>(seqres_line, 12, 17, 0);
  std::vector<std::string> aa = utils::split(seqres_line.substr(17),' ');
  if(sequences.find(chain)==sequences.end()) {
    std::vector<core::chemical::Monomer> v;
    sequences.insert(std::pair<char,std::vector<core::chemical::Monomer>>(chain,v));
  }
  std::vector<core::chemical::Monomer> & v = sequences.at(chain);
  for(const std::string & s:aa) {
    if(s.length() == 3) v.push_back(core::chemical::Monomer::get(s));
    else {
      if(s.length() == 2) v.push_back(core::chemical::Monomer::get(" " + s));
      else v.push_back(core::chemical::Monomer::get("  " + s));
    }
  }
  if(v.size()>n) {
    logger<<utils::LogLevel::SEVERE<<"Corrupted SEQRES data in a PDB file!\n";
  }
}

std::string Modres::to_pdb_line() const {

  return utils::string_format("MODRES %4s %3s %c %4d%c %3s %s", pdb_id.c_str(), res_name.c_str(), chain, residue_id,
      i_code, std_res_name.c_str(), comment.c_str());
}

const std::string Conect::conect_format_2 = std::string("CONECT %4d %4d\n");
const std::string Conect::conect_format_3 = std::string("CONECT %4d %4d %4d\n");
const std::string Conect::conect_format_4 = std::string("CONECT %4d %4d %4d %4d\n");
const std::string Conect::conect_format_5 = std::string("CONECT %4d %4d %4d %4d %4d\n");
const std::string Conect::conect_format_6 = std::string("CONECT %4d %4d %4d %4d %4d %4d\n");


Conect::Conect(const std::string &pdb_line) {
  atom_ids.push_back(utils::from_string<index4>(pdb_line, 8, 12, 0));
}

Conect::Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id) {
  atom_ids.push_back(this_atom_id);
  atom_ids.push_back(bonded_atom_id);
}

Conect::Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id,
               const core::index4 bonded_atom_id_2) {
  atom_ids.push_back(this_atom_id);
  atom_ids.push_back(bonded_atom_id);
  atom_ids.push_back(bonded_atom_id_2);
}

Conect::Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2,
               const core::index4 bonded_atom_id_3) {
  atom_ids.push_back(this_atom_id);
  atom_ids.push_back(bonded_atom_id);
  atom_ids.push_back(bonded_atom_id_2);
  atom_ids.push_back(bonded_atom_id_3);
}

Conect::Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2,
               const core::index4 bonded_atom_id_3, const core::index4 bonded_atom_id_4) {
  atom_ids.push_back(this_atom_id);
  atom_ids.push_back(bonded_atom_id);
  atom_ids.push_back(bonded_atom_id_2);
  atom_ids.push_back(bonded_atom_id_3);
  atom_ids.push_back(bonded_atom_id_4);
}

Conect::Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2,
               const core::index4 bonded_atom_id_3, const core::index4 bonded_atom_id_4,
               const core::index4 bonded_atom_id_5) {
  atom_ids.push_back(this_atom_id);
  atom_ids.push_back(bonded_atom_id);
  atom_ids.push_back(bonded_atom_id_2);
  atom_ids.push_back(bonded_atom_id_3);
  atom_ids.push_back(bonded_atom_id_4);
  atom_ids.push_back(bonded_atom_id_5);
}

std::string Conect::to_pdb_line() const {

  if (atom_ids.size() == 2) return utils::string_format(conect_format_2, atom_ids[0], atom_ids[1]);
  if (atom_ids.size() == 3) return utils::string_format(conect_format_3, atom_ids[0], atom_ids[1], atom_ids[2]);
  if (atom_ids.size() == 4)
    return utils::string_format(conect_format_4, atom_ids[0], atom_ids[1], atom_ids[2], atom_ids[3]);
  if (atom_ids.size() == 5)
    return utils::string_format(conect_format_5, atom_ids[0], atom_ids[1], atom_ids[2], atom_ids[3], atom_ids[4]);
  return utils::string_format(conect_format_6, atom_ids[0], atom_ids[1], atom_ids[2], atom_ids[3], atom_ids[4],
      atom_ids[5]);
}


Atom::Atom(const std::string &pdb_line) {

  serial = utils::from_string<int>(pdb_line, 6, 10, -1);
  is_heteroatom = (pdb_line[0]=='H');
  name = pdb_line.substr(12, 4);
  alt_loc = pdb_line[16];
  res_name = pdb_line.substr(17, 3);
  chain = pdb_line[21];
  residue_id = utils::from_string<int>(pdb_line, 22, 25, -1);
  i_code = pdb_line[26];
  x = utils::from_string<real>(pdb_line, 30, 37, 0.0);
  y = utils::from_string<real>(pdb_line, 38, 45, 0.0);
  z = utils::from_string<real>(pdb_line, 46, 53, 0.0);
  if (pdb_line.length() >= 59) {
    occupancy = utils::from_string<real>(pdb_line, 54, 59, 1.0);
    if (pdb_line.length() >= 65) {
      temp_factor = utils::from_string<real>(pdb_line, 60, 65, 0.0);
      if (pdb_line.length() >= 78) {
        element = pdb_line.substr(76, 2);
        charge = pdb_line.substr(78, 2);
      }
    }
  }
}

std::string Atom::to_pdb_line() const {

  return utils::string_format(atom_format, serial, name.c_str(), alt_loc, res_name.c_str(), chain, residue_id, i_code,
      x, y, z, occupancy, temp_factor, element.c_str(), charge.c_str());
}

const std::string Atom::atom_format = std::string(
    "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s");
const std::string Atom::atom_format_uncharged = std::string(
    "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ");

const std::string Atom::hetatm_format = std::string(
    "HETATM%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s");
const std::string Atom::hetatm_format_uncharged = std::string(
    "HETATM%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ");

std::string Hetatm::to_pdb_line() const {

  return utils::string_format("HETATM%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", serial,
      name.c_str(), alt_loc, res_name.c_str(), chain, residue_id, i_code, x, y, z, occupancy, temp_factor,
      element.c_str(), charge.c_str());
}

core::index4 Pdb::read_coordinates(const std::string & fname, std::vector<core::data::basic::Vec3> & destination,
    const bool if_push_back, const PdbLineFilter & predicate) {

  std::ifstream infile;
  infile.open(fname);
  if (!infile) logger << utils::LogLevel::SEVERE << "Can't open file: " << fname << "\n";
  logger << utils::LogLevel::FILE << "Reading PDB coordinates from: " << fname << "\n";

  return read_coordinates(infile,destination, if_push_back, predicate);
}

core::index4 Pdb::read_coordinates(std::istream &infile, std::vector<core::data::basic::Vec3> &destination,
                                   const bool if_push_back, const PdbLineFilter &predicate) {

  core::index4 nc = 0;
  std::string line;
  if (if_push_back) destination.clear();
  while (std::getline(infile, line)) {
    if ((line.compare(0, 4, "ATOM") != 0) && (line.compare(0, 6, "HETATM") != 0)) continue;  // --- not an atom
    if (!predicate(line)) continue;

    if (!if_push_back) {
      if (nc >= destination.size()) return nc; // --- receiver is full
      destination[nc].x = utils::to_double(line.substr(30, 8).c_str());
      destination[nc].y = utils::to_double(line.substr(38, 8).c_str());
      destination[nc].z = utils::to_double(line.substr(46, 8).c_str());
    } else {
      Vec3 v(utils::to_double(line.substr(30, 8).c_str()),
        utils::to_double(line.substr(38, 8).c_str()),
        utils::to_double(line.substr(46, 8).c_str()));
      destination.push_back(v);
    }
    nc++;
  }

  return nc;
}

std::string Seqres::to_pdb_line() const {

  std::string lines;
  for(const auto & p :sequences) {
    const char chain_id = p.first;
    const std::vector<core::chemical::Monomer> & seq = p.second;
    core::index2 row=1;
    lines += utils::string_format("SEQRES %3d %c%5d ",row,chain_id,seq.size());
    core::index2 n=0;
    for(const core::chemical::Monomer & m : seq) {
      lines += " " + m.code3;
      ++n;
      if(n==13) {
        ++row;
        lines += utils::string_format("\nSEQRES %3d %c%5d ",row,chain_id,seq.size());
        n = 0;
      }
    }
    lines +="\n";
  }
  return lines;
}

static std::string prefixes[] = { "pdb", "PDB", "pdb", "" };
static std::string  sufixes[] = { ".ent", ".ent.gz", ".gz", ".pdb", ".PDB", ".pdb.gz", "" };

std::string find_pdb(const std::string & pdbCode, const std::string & pdbPath) {

  std::string path(pdbPath);
  if (path.back() != utils::dir_separator) path += utils::dir_separator;
  std::string code_lo(pdbCode);
  utils::to_lower(code_lo);
  std::string code_up(pdbCode);
  utils::to_upper(code_up);

  const std::string subdir = code_lo.substr(1, 3) + utils::dir_separator;

  std::vector<std::string> tested;
  for (const std::string & p : prefixes) {
    for (const std::string & s : sufixes) {
      std::string name = path + p + code_lo + s;
      if (utils::if_file_exists(name)) return name;
      else tested.push_back(name);
      name = path + p + code_up + s;
      if (utils::if_file_exists(name)) return name;
      else tested.push_back(name);
      name = path + subdir + p + code_lo + s;
      if (utils::if_file_exists(name)) return name;
      else tested.push_back(name);
      name = path + subdir + p + code_up + s;
      if (utils::if_file_exists(name)) return name;
      else tested.push_back(name);
    }
  }

  static utils::Logger logger("find_pdb");
  logger << utils::LogLevel::WARNING
      << "Cannot find a pdb file for a given code: " + pdbCode + " at path: " + path
          + "\n\tThe following possibilities were tested:\n\t" + utils::to_string(tested," ");

  return "";
}

std::string pdb_code_from_file_name(const std::string & code) {

  std::string code_(utils::basename(code));
  for (const std::string & p : sufixes) utils::replace_substring(code_, p, "");
  for (const std::string & s : prefixes) utils::replace_substring(code_, s, "");

  return code_;
}

bool is_pdb_code(const std::string & code) {

  if (code.length() != 4) return false;
  if (!isdigit(code[0])) return false;
  return true;
}

void write_pdb(const core::data::structural::Structure_SP structure, std::ostream & out) {

  for (auto a_it = structure->first_atom(); a_it != structure->last_atom(); ++a_it)
    out << (*a_it)->to_pdb_line() << "\n";
}


} // ~ io
} // ~ data
} // ~ core

