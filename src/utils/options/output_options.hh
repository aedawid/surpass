/** \file output_options.hh
 * @brief defines option objects (command line flags) that provide output files for programs
 */
#ifndef UTILS_OPTIONS_output_options_HH
#define UTILS_OPTIONS_output_options_HH

#include <utils/options/Option.hh>
#include <utils/Logger.hh>
#include <utils/LogManager.hh>

namespace utils {
namespace options {

static Option output_file("-o", "-out:file", "provide an output file");
static Option output_name_prefix("-out::prefix", "-out::prefix", "a string attached in front of the name of any output file produced by a program");

// ------------- sequence stuff ----------
static Option output_generic("-o", "-out", "redirect the generic output of the program to a file, otherwise it would be printed on the screen");
static Option output_fasta("-of", "-out:fasta", "provide an output file to write sequences in FASTA format", "stdout");
static Option output_fasta_secondary("-fasta_ss", "-out:fasta:secondary", "write secondary structure in FASTA format when available");
static Option output_blocks("-out:sequence:blocks", "-out:sequence:blocks", "provide an output file to write gapless blocks found in a given alignment(s)");
static Option output_gaps("-out:sequence:gaps", "-out:sequence:gaps", "provide an output file to write the location of each gapped fragment found in a given alignment(s)");
static Option output_seq_width("-w", "-out:sequence:width", "break each sequence to lines of this width");
static Option out_profile_txt("-out:profile:txt", "-out:profile:txt", "writes a sequence profile in as a flat text table");

// ------------- structure stuff ----------
static Option output_pdb("-op", "-out:pdb", "provide an output file to write structure in PDB format");
static Option output_pdb_min("-out:pdb:min_en", "-out:pdb:min_en", "provide an output file to write low-energy structures in PDB format");
static Option output_pdb_min_value("-out:pdb:min_en::value", "-out:pdb:min_en::value", "the highest energy value for a structure to be recorded with -out:pdb:min_en option");
static Option output_pdb_min_fraction("-out:pdb:min_en::fraction", "-out:pdb:min_en::fraction", "say 0.15 to record structures worse by 15% of energy than the currently lowest ");
static Option output_trax("-ox", "-out:trax", "provide a file name to write output trajectory in TRAX format");
static Option output_pdb_header("-out:pdb:header", "-out:pdb:header", "write a header when writing a PDB file");
static Option out_sse("-sse","-out:sse",  "prints a list of secondary structure elements");

// -------------  Save alignment as an abstract path -------------
static Option output_alipath("-out:alignment:as_path", "-out:alignment:as_path", "write alignment as a path through an alignment matrix");

/********** Options that define representation for output structures **********/
static Option define_output_representation("-orep", "-out:representation:cfg_file", "define representation in which structures will be written: the representation will be defined by a cfg file");

//  -------------  HSP hits output
static Option output_blast_nhits("-n", "-out:blast:nhits", "print only certain number of hits, ordered by the selected sorting criteria (e-value is used by default)");

}
}

#endif
