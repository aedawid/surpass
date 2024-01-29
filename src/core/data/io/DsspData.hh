#ifndef CORE_DATA_IO_DsspData_H
#define CORE_DATA_IO_DsspData_H

#include <string>
#include <iostream>
#include <vector>

#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/structural/Residue.hh>
#include <core/data/structural/Chain.hh>
#include <utils/Logger.hh>

namespace core {
namespace data {
namespace io {

/** @brief Stores data extracted form a single DSSP line
 */
class DsspDataLine {
public:
  int position; ///< cumulative position in the structure
  int seq_position; ///< position in the chain
  char icode; ///< insertion code for the residue
  char chain_letter; ///< chain_id (a single character)
  char residue_letter; ///< residue type. (in a one-letter code)
  char structure; ///< structure type (H,E,T,G,B.. - DSSP code)
  int bp1; ///< first bridge partner resnum.
  int bp2; ///< second bridge partner resnum.
  int acc; ///< residue water exposed surface in Angstrom**2.

  int nho_partner1; ///< N-H-->O partner (1)
  double nho_energy1; ///< N-H-->O E(kcal/mol) (1)
  int ohn_partner1; ///< N-H-->O partner (1)
  double ohn_energy1; ///< N-H-->O E(kcal/mol) (1)
  int nho_partner2; ///< N-H-->O partner (2)
  double nho_energy2; ///< N-H-->O E(kcal/mol) (2)
  int ohn_partner2; ///< N-H-->O partner (2)
  double ohn_energy2; ///< N-H-->O E(kcal/mol) (2)
  double tco; ///< tco angle around this residue
  /** \brief kappa angle
   *
   * virtual bond angle (bend angle) defined by the three
   * C-alpha atoms of residues I-2,I,I+2. Used to define bend (structure code 'S').
   */
  double kappa;
  /** \brief alpha angle
   *
   * Virtual torsion angle (dihedral angle) defined by
   * the four C-alpha atoms of residues I-1,I,I+1,I+2.
   * Used to define chirality (structure code '+' or '-').
   */
  double alpha;
  double phi; ///< IUPAC Phi peptide backbone torsion angles
  double psi; ///< IUPAC Psi peptide backbone torsion angles
  double xCa; ///< X coordinate of CA atom
  double yCa; ///< Y coordinate of CA atom
  double zCa; ///< Z coordinate of CA atom

  /** \brief Constructor parses a single line from a DSSP file and initializes fields of this class
   * @param line - a single line in the DSSP format
   * @param ifClassicFormat - if true, the line will be parsed as in the classic (original) format. Otherwise
   * the EMBL (the wider) format will be assumed.
   */
  DsspDataLine(const std::string & line, bool ifClassicFormat);

  /**\brief Returns true if a given residue matches this DSSP data line.
   *
   * This method tests is residue ID, residue insertion code and chain ID are the same
   * @param r - a residue
   */
  bool is_my_residue(const core::data::structural::Residue & r) const;

  /// Creates a string in the DSSP classic format based on the data stored in this object
  std::string  toString() {

    return utils::string_format("%5d %4d %c %c %c %8.3f %8.3f %8.3f %8.3f %8.3f",position,seq_position,
        chain_letter,residue_letter,structure,tco,alpha,kappa,phi,psi);
  }

private:
  /// Start column for each field of the line - the new format. Starts at 1.
  static const std::vector<int> start_column;
  /// End column for each field of the line - the new format. Starts at 1.
  static const std::vector<int> end_column;
  /// Start column for each field of the line - the classic format
  static const std::vector<int> start_column_classic;
  /// End column for each field of the line - the classic format
  static const std::vector<int> end_column_classic;
  static const int offset;
};

/** @brief Reads and parses a DSSP file.
 *
 * The following example converts DSSP to fasta:
 * @include ex_DsspData.cc
 *
 * The second example converts DSSP to SS2:
 * @include ex_dssp_to_ss2.cc
 */
class DsspData {
public:
  /** \brief Whether the input DSSP file is in the new or in the old classic format.
   *
   * If the flag is true, the parser will attempt to parse the old 'classic' file format.
   * Otherwise, the CMBI (the wider) file format will be assumed
   */
  bool if_classic_format;

  /// Constructor reads and parses the data
  DsspData(const std::string & file_name, const bool if_classic_format = true);

  /// Creates a secondary structure object for every chain found in the DSSP file
  std::vector<core::data::sequence::SecondaryStructure> create_sequences() const;

  /// Returns PDB code of the source protein - if available
  const std::string & code() const { return code_; }

  /// Converts a single character from DSSP convention to the HEC code
  static char dssp_to_hec(const char dsspChar);

private:
  std::string code_;
  utils::Logger logger;
  unsigned char attempts = 0;
  std::string file_name_;
  std::vector<DsspDataLine> data;
  int highest_residue_index = -1;
  void load(std::ifstream & in, const bool if_classic_format);
};

}
}
}

#endif
/**
 * @example ex_DsspData.cc
 * @example ex_dssp_to_ss2.cc
 */

