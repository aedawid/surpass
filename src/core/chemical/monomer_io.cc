#include <memory>
#include <string>
#include <unordered_map>

#include <core/SURPASSenvironment.hh>
#include <core/chemical/monomer_io.hh>
#include <core/chemical/Monomer.hh>
#include <core/data/io/Cif.hh>

#include <utils/Logger.hh>
#include <utils/string_utils.hh>

namespace core {
namespace chemical {

const std::string comp_id = "_chem_comp.id";
const std::string comp_name = "_chem_comp.name";
const std::string comp_type = "_chem_comp.type";
const std::string comp_formula = "_chem_comp.formula";
const std::string comp_parent = "_chem_comp.mon_nstd_parent_comp_id";
const std::string comp_charge = "_chem_comp.pdbx_formal_charge";
const std::string comp_ambig = "_chem_comp.pdbx_ambiguous_flag";
const std::string comp_code1 = "_chem_comp.one_letter_code";
const std::string comp_code3 = "_chem_comp.three_letter_code";

using namespace core::data::io;
utils::Logger monomer_io_logs = utils::Logger("monomer_io");

void read_monomers_cif(const std::string &cif_filename) {

  core::data::io::Cif reader(cif_filename);
  std::shared_ptr<DataBlock> b = 0;
  core::index2 index = 10000;
  core::index2 parent_code = 10000;
  const int beg = std::distance(Monomer::cbegin(), Monomer::cend());
  do {
    b = reader.read_block();
    if (b != 0) {
      char code1 = b->tokens[comp_code1][0];
      std::string code3 = b->tokens[comp_code3];
      if (code3.length() == 2) code3 = ' ' + code3;
      if (code3.length() == 1) code3 = "  " + code3;
      std::string parent = b->tokens[comp_parent];
      std::string type = b->tokens[comp_type];
      utils::to_upper(type);
      char c_type = 'U';
      if (type.find("PEPTIDE") != std::string::npos) c_type = 'P';
      if ((type.find("DNA") != std::string::npos) || (type.find("RNA") != std::string::npos)) {
        c_type = 'N';
        code1 = tolower(code1);
      }
      if (type.find("SACHARIDE") != std::string::npos) c_type = 'S';

      if (Monomer::is_known_monomer(parent)) parent_code = Monomer::get(parent).id;
      else {
        parent_code = Monomer::UNL.id;
      }
      if (Monomer::is_known_monomer(code3)) index = Monomer::get(code3).id;
      else {
        index = parent_code + 1000;
      }
      bool ambig_flag = utils::from_string<bool>(b->tokens[comp_ambig]);
      core::real charge = utils::from_string<core::real>(b->tokens[comp_charge]);
//			blocks.push_back(b);
      Monomer m(index, code1, code3, c_type, 255, 255, ambig_flag, charge, parent_code);
      Monomer::register_monomer(m);
    }
  } while (b != 0);
  monomer_io_logs << utils::LogLevel::INFO << (int) (std::distance(Monomer::cbegin(), Monomer::cend()) - beg)
      << " new monomers loaded from: " << cif_filename << "\n";
}

void read_monomers_binary(const std::string & filename) {

  Monomer m(0, ' ', "   ", ' ', 0, 0, false, 0.0, 0);
  monomer_io_logs << utils::LogLevel::FILE << "Loading all monomer definitions from a binary file: " << filename
      << "\n";
  const int beg = std::distance(Monomer::cbegin(), Monomer::cend());
  std::ifstream file;
  file.open(filename, std::ios::app | std::ios::in);
  file.seekg(0, std::ios::beg);
  while (!file.eof()) {
    if (!file.read(reinterpret_cast<char*>(const_cast<Monomer*>(&m)), sizeof(Monomer))) break;
    if (!Monomer::is_known_monomer(m.code3)) Monomer::register_monomer(m);
  }
  file.close();
  monomer_io_logs << utils::LogLevel::INFO << (int) (std::distance(Monomer::cbegin(), Monomer::cend()) - beg)
      << " new monomers loaded from: " << filename << "\n";
}

void write_monomers_binary(const std::string & filename) {

  Monomer m(0, ' ', "   ", ' ', 0, 0, false, 0.0, 0);
  std::fstream file;
  file.open(filename, std::ios::app | std::ios::out | std::ios::in);
  for (auto it = Monomer::cbegin(); it != Monomer::cend(); ++it) {
    file.write(reinterpret_cast<char*>(const_cast<Monomer*>(&(*it))), sizeof(Monomer));
  }
  file.close();
}

void read_monomers_txt(const std::string &txt_filename) {

  std::ifstream file;
  file.open(txt_filename);
  std::string line;
  while (std::getline(file, line)) {
    Monomer m(line);
    Monomer::register_monomer(m);
  }
}

void write_monomers_txt(const std::string &txt_filename) {

  std::ofstream file;
  file.open(txt_filename);
  for (auto it = Monomer::cbegin(); it != Monomer::cend(); ++it) {
    file
        << utils::string_format("%5d %2d %c %c %3d %1d %4.1f %3s\n", it->id, it->parent_id, it->code1, it->type, it->n_atoms, it->is_ambiguous,
            it->charge, it->code3.c_str());
  }
  file.close();
}

void load_monomers_from_db() {
  read_monomers_txt(core::SURPASSenvironment::from_file_or_db("monomers.txt"));
}

}
}
