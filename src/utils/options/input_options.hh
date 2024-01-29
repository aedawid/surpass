/** \file input_options.hh
 * @brief defines option objects (command line flags) that provide input data to programs
 */
#ifndef UTILS_OPTIONS_input_options_HH
#define UTILS_OPTIONS_input_options_HH

#include <utils/options/Option.hh>
#include <utils/options/OptionParser.hh>
#include <utils/Logger.hh>
#include <utils/LogManager.hh>
#include <core/SURPASSenvironment.hh>

namespace utils {
namespace options {

static Option input_file("-i", "-in:file", "provide an input file");
static Option data_column("-c", "-in:column", "which data column to use", "1");
static Option data_columns("-cc", "-in:columns", "which data columns to use", "1,2,3");
static Option movers_config("-m", "-in:movers", "provide a config file that defines movers");

/********** Load additional monomers from a file **********/
static Option input_cif_monomers("-in:monomers:cif", "-in::monomers::cif", "load definitions of additional monomers (standard CIF format)");
static Option input_bin_monomers("-in:monomers:bin", "-in::monomers::bin", "load definitions of additional monomers (internal binary format)");
static Option input_txt_monomers("-in:monomers:txt", "-in::monomers::txt", "load definitions of additional monomers (internal text format)");

/********** Load sequence(s) or a sequence profile**********/
static Option input_fasta("-if", "-in:fasta", "provide an input file in FASTA format");
static Option input_pir("-in::pir", "-in:pir", "provide an input file in PIR format");
static Option input_ss2("-in::ss2", "-in:ss2", "provide an input secondary structure in PsiPred's SS2 format");
static Option input_clustalw("-iw", "-in:clustalw", "provide an input file in ClustalW format");
static Option input_native_fasta("-ifn", "-in:fasta:native", "provide the native (or reference) sequence (or alignment) in FASTA format");
static Option input_chk("-ib", "-in:profile:chk", "provide an input sequence profile in the binary CHK format (legacy blastpgp output)");
static Option input_asn1("-ia", "-in:profile:asn1", "provide an input sequence profile in the ASN.1 format (blast+ output)");

static Option input_n_atoms("-n", "-in:n_atoms", "the number of atoms in the input structure(s)");

/********** Load file(s) in PDB format **********/
static Option input_pdb_path("-ippath", "-in:pdb:path", "search for pdb files in this directory");
static Option input_pdb("-ip", "-in:pdb", "provide an input protein structure(s) in PDB format");
static Option input_pdb_header("-in:pdb:header", "-in:pdb:header", "parse header from input PDB files (turned off by default)");
static Option input_models("-in:pdb:models", "-in::pdb::models", "provide an input protein structure(s) in PDB format");
static Option input_pdb_sort("-pdb_sort", "-in:pdb:sort", "sort atoms, residues and chains found in input PDB structures");
static Option input_pdb_hydrogens("-pdb_sort", "-in:pdb:hydrogens", "read in also hydrogen atoms");
static Option input_trax("-ix", "-in:trax", "input trajectory in TRAX format");

/********** Load abstract alignment **********/
static Option input_alipath("-alipath", "-in:alignment:as_path", "provide input alignment as a path through an alignment matrix");

/// Load PsiBlast output
static Option input_blast_output("-blastout", "-in:blast:outfile", "load results file produced by (psi)blast program");
static Option input_blast_iteration("-blast_iter", "-in:blast:which_iteration", "load results only from the requested iteration of psiblast");
static Option input_blast_last_iteration("-blast_last_iter", "-in:blast:last_iteration", "load results only from the last iteration of psiblast");

/** @name Input alignment files
 */
///@{
static Option input_alignment_fasta("-in:alignment:fasta", "-in::alignment::fasta", "input alignment in FASTA");
static Option input_alignment_edinburgh("-in:alignment:edinburgh", "-in::alignment::edinburgh", "input alignment in Edinburgh format");
static Option input_alignment_pir("-in:alignment:pir", "-in::alignment::pir", "input alignment in PIR format");
///@}

/* \brief Load a native (reference) structure from a PDB file.
 *
 * The snippet below shows the example usage of this option:
 * @code
 *   Structure_SP native = nullptr;
 *   if (utils::options::input_pdb_native.was_used()) {
 *     core::data::io::Pdb pdb(option_value<std::string>(input_pdb_native), core::data::io::is_ca);
 *     native = pdb.create_structure(0);
 *   }
 * @endcode
 *
 * The more elaborate code (which also uses PDB line filters) is available as: <code>utils::options::native_from_cmdline()</code>
 * declared in input_utils.hh
 */
static Option input_pdb_native("-ipn", "-in:pdb:native", "provide the native (or reference) protein structure in PDB format");

static Option input_pdb_list("-ipl", "-in:pdb:listfile",  "read all PDB files listed in the given listfile (as a single column with file names)");

static Option input_chainid_list("-icl", "-in:pdb:chainid::list",  "read PDB files and extract chains; input IDs should look like 2azaA or 2aza:A");

static Option input_pdb_modelslist("-modelslist", "-in:pdb:modelslist", "read all PDB model files listed in the given listfile (as a single column)");

static Option input_pdb_models("-models", "-in:pdb:models", "read all PDB structures as models (all models must be exactly the same molecule)");

/********** Filters for the PDB formatted input **********/
static Option all_models("-all_models", "-in:pdb:all_models", "read all models from the given PDB files (otherwise only the first model will be loaded)");
static Option ca_only("-ca_only", "-in:pdb:ca_only", "read only CA atoms from the given PDB files");
static Option include_hydrogens("-with_h", "-in:pdb:with_hydrogens", "read in all hydrogen atoms while reading the given PDB files; by default it does not use these atoms");
static Option bb_only("-bb_only", "-in:pdb:bb_only",  "read only backbone heavy atoms from the given PDB files");
static Option keep_alternate("-no_alt", "-in:pdb:no_alt",  "remove all alternate locations for atoms, leave just one variant per atom (by default apps keep only variant 'A' if there is more than one)");
static Option keep_water("-water", "-in:pdb:keep_water",  "do not remove water molecules when reading in a PDB file");

class SetDBPath : public Option {
public:
  SetDBPath() :
      Option("-d", "-in:database", "path to parameters directory", "", false) {
  }
  virtual void execute() {
    core::SURPASSenvironment::surpass_db_path(option_value<std::string>("-in:database"));
  }
}static db_path;

}
}

#endif
