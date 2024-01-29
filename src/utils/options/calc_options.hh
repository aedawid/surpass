/** \file calc_options.hh
 * @brief defines option objects (command line flags) that switch on calculations
 */
#ifndef UTILS_OPTIONS_calc_options_HH
#define UTILS_OPTIONS_calc_options_HH

#include <utils/options/Option.hh>
#include <utils/Logger.hh>
#include <utils/LogManager.hh>

namespace utils {
namespace options {

static Option sequence_identity("-calc:sequence:identity", "-calc::sequence::identity", "calculate sequence identity between all possible sequence pairs");
static Option calc_aln_k("-calc:alignment:aln_k", "-calc::alignment::aln_k", "calculate Aln_k measure between two pairwise alignment");

/** @name Crmsd and superimposition related options
 *  Control how superimpositions are performed and crmsd computed
 */
///@{
static Option calc_crmsd("-calc:crmsd", "-calc::crmsd", "calculate crmsd between structures");
static Option calc_crmsd_matching_atoms("-calc:crmsd:matching_atoms", "-calc::crmsd::matching_atoms", "define atoms for the superposition on which the crmsd will be computed");
static Option calc_crmsd_matching_atoms_file("-calc:crmsd:matching_atoms:file", "-calc::crmsd::matching_atoms::file", "load values for -calc:crmsd:matching_atoms option from a file");
static Option calc_crmsd_rotated_atoms("-calc:crmsd:rotated_atoms", "-calc::crmsd::rotated_atoms", "define atoms that will be transformed by the rototranslation transformation");
static Option calc_crmsd_rotated_atoms_file("-calc:crmsd:rotated_atoms:file", "-calc::crmsd::rotated_atoms::file", "load -calc::crmsd::rotated_atoms option values from a file");
///@}

/** @name Options to calculate tmscore
 */
///@{
static Option calc_tmscore("-calc:tmscore", "-calc::tmscore", "calculate tmscore between structures");
static Option tmscore_ref_length("-calc:tmscore:length", "-calc::tmscore::length", "reference length for tmscore calculations");
///@}

/// Use DSSP algorithm
static Option calc_dssp("-dssp", "-calc::dssp", "use DSSP algorithm to calculate protein secondary structure and architecture");

/// Option to re-orient a structure
static Option orient_structure("-calc:orient", "-calc::orient", "reorients each structure so it is placed in the origin");

/** @name Options for simple structural calculations
 *  Enforce calculations of structural properties: distances, angles, etc.
 */
///@{
static Option calc_a13("-a13","-calc:a13",  "calculate A13 - the planar angle between three subsequent CA");
static Option calc_r13("-r13","-calc:r13",  "calculate R13 - the distance between i-th and (i+2)nd CA");
static Option calc_t14("-t14","-calc:t14",  "calculate phi and psi dihedral angles");
static Option calc_r14("-r14","-calc:r14",  "calculate R14 - the distance between i-th and (i+3)rd CA");
static Option calc_r15("-r15","-calc:r15",  "calculate R15 - the distance between i-th and (i+4)th CA");
static Option calc_chi("-chi","-calc:chi",  "calculate all chi dihedral angles");
static Option calc_phi_psi("-phi_psi","-calc:phi_psi",  "calculate all Phi and Psi dihedral angles");
static Option calc_tau("-tau","-calc:tau",  "calculate all tau (N-CA-C) planar angles");
static Option calc_omega("-omega","-calc:omega",  "calculate all omega dihedral angles");
static Option calc_abego("-abego","-calc:abego",  "finds ABEGO classification for each residue");
static Option calc_lambda("-lambda","-calc:lambda",  "calculates lambda angles i.e. dihedrals between a peptide plate and a respective plane of 3 CA atoms");
static Option calc_local("-calc::local","-calc:local",  "calculates all the local structural properties of a protein chain");
static Option calc_hbonds_bb("-bb_hbonds","-calc:hbonds:bb",  "detects plausible hydrogen bonds within a protein backbone");
///@}

/** @name Options for computing various distance maps
 */
///@{
static Option calc_distmap_ca("-calc:distmap:ca", "-calc::distmap::ca", "calculate distance map based on CA atoms");
static Option calc_distmap_min("-calc:distmap:min", "-calc::distmap::min", "calculate map of minimal interatomic distances between residues");
static Option calc_distmap_allatom("-calc:distmap:allatom", "-calc::distmap::allatom", "calculate distance map based on all atoms in a given structure");
static Option calc_distmap_described("-calc:distmap:describe", "-calc::distmap::describe", "print detailed information for atoms whose distances are reported");
static Option longest_distance_shown("-calc:distmap:shorter_than", "-calc::distmap::shorter_than", "print only short range distances");
static Option distmap_seq_separation("-calc:distmap:seq_separation", "-calc::distmap::seq_separation", "print only distances from residues separated by at least this number of other residues along a sequence");
///@}
}
}

#endif
