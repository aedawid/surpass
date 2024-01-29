#include <iostream>

#include <core/data/io/ss2_io.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/io/Pdb.hh>
#include <simulations/representations/surpass_utils.hh>

#include <core/chemical/Monomer.hh>
#include <core/chemical/monomer_io.hh>
#include <core/data/structural/Structure.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/PdbAtom.hh>

using namespace core::data::structural;
using namespace core::data::sequence;
using namespace core::data::io;


void print_usage(const char* program_name) {
    std::cerr 
              << "Reads an all-atom structure from a PDB file and produces a structure in SURPASS representation.\n\n"
              << "Options:\n"
              << "  -h, --help          Show this help message and exit\n"
              << "Arguments:\n"
              << "  <PDB file>          Input PDB file containing the structure.\n"
              << "  [SS2 file]          Optional SS2 file for secondary structure prediction.\n\n"
              << "Usage: " << program_name << " <PDB file> [SS2 file] > output.pdb\n";
}


/** @brief Reads an all-atom structure from a PDB file and produces a structure in SURPASS representation.
 */
int main(const int argc, const char* argv[]) {

  // Check for help flag or no input
  if (argc == 1 || (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
    print_usage(argv[0]);
    return 1;
  }

  // --- Read the input PDB and create a structure object
  core::data::io::Pdb reader(argv[1], is_not_alternative, true);
  core::data::structural::Structure_SP strctr = reader.create_structure(0);

  // --- Check whether loaded structure is in the SURPASS representation
  if (simulations::representations::is_surpass_model(*strctr)) {
    std::cerr<<"Loaded structure of "<<argv[1]<<" has SURPASS representation! Load fullatom model.\n";
  }
  else {
  // --- Convert the Structure into SURPASS representation and write the result in the PDB format
    Structure_SP structure_sp = simulations::representations::surpass_representation(*strctr);

    if (argc==2) {
      std::cerr<<"The secondary structure of "<<argv[1]<<" is based on header from .pdb file.\n";
      for (auto atom_sp = structure_sp->first_atom(); atom_sp != structure_sp->last_atom(); ++atom_sp)
        std::cout << (*atom_sp)->to_pdb_line() << "\n";
    } else if (argc==3) {
      std::cerr<<"The secondary structure of "<<argv[1]<<" is based on prediction from .ss2 file.\n";
      core::data::sequence::SecondaryStructure_SP ss2 = core::data::io::read_ss2(argv[2],"");
      core::index2 id=0;
      core::data::sequence::SecondaryStructure_SP ss2_surpass = simulations::representations::surpass_representation(*ss2);
      for (auto atom_sp = structure_sp->first_atom(); atom_sp != structure_sp->last_atom(); ++atom_sp) {
        if ((*ss2_surpass).ss(id)=='H') {
	  (*atom_sp)->atom_name(" H  ");
	  (*atom_sp)->owner()->ss('H');
        } else if ((*ss2_surpass).ss(id)=='E') {
	  (*atom_sp)->atom_name(" S  ");
	  (*atom_sp)->owner()->ss('E');
        } else if ((*ss2_surpass).ss(id)=='C') {
	  (*atom_sp)->atom_name(" C  ");
	  (*atom_sp)->owner()->ss('C');
        }
        ++id;
        std::cout << (*atom_sp)->to_pdb_line() << "\n";
      }
    }
    auto prev_atom_sp = structure_sp->first_atom();
    for (auto atom_sp = (++(structure_sp->first_atom())); atom_sp != structure_sp->last_atom(); ++atom_sp) {
      core::data::io::Conect cn((*prev_atom_sp)->id(),(*atom_sp)->id());
      std::cout << cn.to_pdb_line();
      ++prev_atom_sp;
    }
  }
}
