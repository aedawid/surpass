#include <iostream>
#include <iomanip>
#include <core/data/io/fasta_io.hh>
#include <core/data/io/Pdb.hh>
#include <utils/string_utils.hh>


void print_usage(const char* program_name) {
    std::cerr << "Reads an all-atom structure from a PDB file and produces FASTA sequence file.\n\n"
              << "Options:\n"
              << "  -h, --help          Show this help message and exit\n"
              << "Arguments:\n"
              << "  <PDB file>          Input PDB file containing the structure.\n\n"
              << "Usage: " << program_name << " <PDB file> > output.fasta\n";
}

int main(const int argc, const char* argv[]) {

  if (argc == 1 || (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
    print_usage(argv[0]);
    return 1;
  }

  using namespace core::data::io;
  using namespace core::data::sequence;
  using core::data::sequence::Sequence_SP;
  using namespace core::data::io;
      
  core::data::io::Pdb reader(argv[1],is_standard_atom,true);
  core::data::structural::Structure_SP strctr = reader.create_structure(0);

// Iterate over all residues in the structure
  int n = 0;
  int N = strctr->count_residues();
  std::cout << ">"<<reader.pdb_code()<<", length: "<<N<<"\n";
  for (auto it_resid = strctr->first_residue(); it_resid!=strctr->last_residue(); ++it_resid) {
    const core::chemical::Monomer & m = (*it_resid)->residue_type();
    if (m.type == 'P') {
      ++n;
      if (n==60) {
        std::cout<<utils::string_format("%c", m.code1)<<"\n";
        n = 0;
      } else std::cout<<utils::string_format("%c", m.code1);
    }
  }
}
