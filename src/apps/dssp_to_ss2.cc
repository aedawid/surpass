#include <iostream>
#include <core/data/io/ss2_io.hh>
#include <core/data/io/DsspData.hh>


void print_usage(const char* program_name) {
    std::cerr << "Reads the output from DSSP and produces the secondary structure assignment (in PsiPred format).\n\n"
              << "Options:\n"
              << "  -h, --help          Show this help message and exit\n"
              << "Arguments:\n"
              << "  <DSSP file>         Input DSSP file containing the secondary structure assignment.\n\n"
              << "Usage: " << program_name << " <DSSP file> > output.ss2\n";
}


/** @brief Reads a DSSP file and prints the secondary structure of each chain in SS2 format.
 *@see ex_DsspData.cc converts DSSP to FASTA format
 */
int main(const int argc, const char* argv[]) {

  if (argc == 1 || (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
    print_usage(argv[0]);
    return 1;
  }

  core::data::io::DsspData dssp(argv[1], true);
  for (const auto & ss2 : dssp.create_sequences())
    core::data::io::write_ss2(ss2,std::cout);
}
