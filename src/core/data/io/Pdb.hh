/** \file Pdb.hh
 * @brief provides classes and methods for reading and parsing PDB files
 */
#ifndef CORE_DATA_IO_PDB_H
#define CORE_DATA_IO_PDB_H

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <functional>
#include <iostream>

#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/structural/Structure.hh>
#include <core/data/structural/Residue.fwd.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/io/PdbField.hh>

#include <utils/Logger.hh>

namespace core {
namespace data {
namespace io {

using core::real;

/** @brief A data type for predicated that filter PDB lines (given as strings) while a PDB file is loaded in.
 *
 * This will remove unwanted atoms up-front prior parsing PDB data making this process faster
 */
typedef std::function<bool(const std::string&)> PdbLineFilter;

/** @name Predicates for filtering PDB lines.
 *
 * A few handy predicates were declared as a static const pointer to a function (aka  PdbLineFilter instances)
 * Each of these predicates takes a string as an argument and returns true if the string is a PDB line of a particular type.
 * The purpose for the predicates is to filter PDB data prior parsing and hence speed up the PDB reading process.
 *
 * <strong>Note:</strong> any predicate does not check if a given line is actually in the PDB format; this is taken for granted.
 *
 * Given an example string:
 * @code
std::string line = "ATOM    181  C   LYS A  24     -17.473  23.415  -7.068  1.00 14.12           C";
 * @endcode
 *
 * the following statement:
 * @code
 * core::data::io::is_bb(line);
 * @endcode
 * returns true. Moreover, also the following predicated:
 *     - core::data::io::is_standard_atom
 *     - core::data::io::is_not_hydrogen
 *     - core::data::io::is_not_water
 *     - core::data::io::is_not_alternative
 *     - core::data::io::keep_all
 *
 * will return <code>true</code>. But core::data::io::is_hetero_atom and core::data::io::is_ca will return <code>false</code>.
 *
 * The following simple example removes or water molecules and alternate atom locations while reading a PDB file:
 * \include ex_PdbLineFilter.cc
 */
///@{

/** @brief A method that creates a chain filter to read-in a particular chain
 *
 * @tparam C - the template code (single character)
 */
template<char C>
const PdbLineFilter create_chain_filter() {
  auto f = [&](const std::string line) {
    if(line[21]!=C) return false;
    return true;
  };
  return f;
}

/** @brief Predicate that negates another filter
 */
const PdbLineFilter invert_filter(const PdbLineFilter f1);

/** @brief Predicate that returns true if both <code>f1</code> and <code>f2</code> returned true.
 */
const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2);

/** @brief Predicate that returns true if all the predicates returned true.
 */
const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2, const PdbLineFilter f3);

/** @brief Predicate that returns true if all the predicates returned true.
 */
const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2, const PdbLineFilter f3, const PdbLineFilter f4) ;

/** @brief Predicate that returns true if all the predicates returned true.
 */
const PdbLineFilter all_true(const PdbLineFilter f1, const PdbLineFilter f2, const PdbLineFilter f3, const PdbLineFilter f4, const PdbLineFilter f5);

/** @brief Predicate that returns true if a given PDB line holds a standard atom.
 */
static const PdbLineFilter is_standard_atom = [](const std::string line) {
  if(line[0]!='A') return false;
  if(line[1]!='T') return false;
  return true;
};

/** @brief Predicate that returns true if a given PDB line holds a hetero-atom
 */
static const PdbLineFilter is_hetero_atom = [](const std::string line) {
  if(line[0]!='H') return false;
  if(line[3]!='A') return false;
  return true;
};

/** @brief Predicate that returns true if a given PDB line holds a non-hydrogen atoms
 */
static const PdbLineFilter is_not_hydrogen = [](const std::string line) {

  if(line[13]=='H') return false;
  if(line[12]=='H') return false;
  return true;
};

/** @brief Predicate that returns true if a given PDB line holds a non-hydrogen atoms
 */
static const PdbLineFilter is_not_water = [](const std::string line) {
  if((line[17]=='H')&&(line[18]=='O')&&(line[19]=='H')) return false;
  return true;
};

/** @brief Predicate that returns true if a given PDB line holds a non-hydrogen atoms
 */
static const PdbLineFilter is_water = [](const std::string line) {
  if((line[17]=='H')&&(line[18]=='O')&&(line[19]=='H')) return true;
  return false;
};

/** @brief Predicate that returns true if a given PDB line holds an atom that has no alternatives OR is the first of the set of alternative atom locations
 */
static const PdbLineFilter is_not_alternative = [](const std::string line) {
  if((line[16]=='A')||(line[16]==' ')) return true;
  return false;
};

/** @brief Predicate that returns true if a given PDB line holds C\f$\alpha\f$ atom; false otherwise
 */
static const PdbLineFilter is_ca = [](const std::string line) {
  if(line[13]!='C') return false;
  if(line[14]!='A') return false;
  return true;
};

/** @brief Predicate that returns true if a given PDB line holds one of the backbone heavy atoms; false otherwise
 */
static const PdbLineFilter is_bb = [](const std::string line) {
  if((line[17]=='H')&&(line[18]=='O')) return false; // --- Skip water molecules
  if((line[13]=='N')&&(line[14]==' ')) return true;
  if((line[13]=='C')&&(line[14]=='A')) return true;
  if((line[13]=='C')&&(line[14]==' ')) return true;
  if((line[13]=='O')&&(line[14]==' ')) return true;
  return false;
};

/** @brief Predicate that always returns true.
 *
 * When this predicate is used, all PDB lines are loaded
 */
static const PdbLineFilter keep_all = [](const std::string line) {
  return true;
};
///@}


class Atom : public PdbField {
public:
  int serial;
  int residue_id;
  std::string name;
  char alt_loc;
  char i_code;
  std::string res_name;
  real x, y, z;
  real occupancy, temp_factor;
  char chain;
  bool is_heteroatom;
  std::string element;
  std::string charge;

  Atom(const std::string &pdb_line);
  std::string to_pdb_line() const;
  const static std::string atom_format;
  const static std::string atom_format_uncharged;
  const static std::string hetatm_format;
  const static std::string hetatm_format_uncharged;
};

class Hetatm : public Atom {
public:
  Hetatm(const std::string &pdb_line) :
      Atom(pdb_line) {
  }
  std::string to_pdb_line() const;
};

class Header : public PdbField {
public:
  std::string classification;
  std::string date;
  std::string pdb_id;

  Header(const std::string & pdb_line) :
      classification(pdb_line.substr(10, 40)), date(pdb_line.substr(50, 9)), pdb_id(pdb_line.substr(62, 4)) {
  }
  std::string to_pdb_line() const;
};

class Modres : public PdbField {
public:
  int residue_id;
  char chain;
  char i_code;
  std::string pdb_id;
  std::string res_name;
  std::string std_res_name;
  std::string comment;

  Modres(const std::string &pdb_line);
  std::string to_pdb_line() const;
};

class Conect : public PdbField {
public:
  std::vector<index4> atom_ids;

  Conect(const std::string &pdb_line);

  Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id);

  Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2);

  Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2,
         const core::index4 bonded_atom_id_3);

  Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2,
         const core::index4 bonded_atom_id_3, const core::index4 bonded_atom_id_4);

  Conect(const core::index4 this_atom_id, const core::index4 bonded_atom_id, const core::index4 bonded_atom_id_2,
         const core::index4 bonded_atom_id_3, const core::index4 bonded_atom_id_4,
         const core::index4 bonded_atom_id_5);

    std::string to_pdb_line() const;

  const static std::string conect_format_2;
  const static std::string conect_format_3;
  const static std::string conect_format_4;
  const static std::string conect_format_5;
  const static std::string conect_format_6;
};

class Keywords : public PdbField {
public:
  std::vector<std::string> keywords;

  Keywords(const std::string &pdb_line);
  std::string to_pdb_line() const;
};


class TVect : public PdbField {
public:
  index2 serial;
  double x;
  double y;
  double z;

  TVect(const std::string &pdb_line);

  TVect(const index2 serial, const double x, const double y, const double z) : serial(serial), x(x), y(y), z(z) {}

  std::string to_pdb_line() const;
};

class DBRef : public PdbField {
public:
  std::string pdb_code;
  char chain;
  int residue_from;
  char insert_from;
  int residue_to;
  char insert_to;
  std::string db_name;
  std::string db_accession;
  std::string db_code;
  int db_residue_from;
  char db_insert_from;
  int db_residue_to;
  char db_insert_to;

  DBRef(const std::string &pdb_line);
  std::string to_pdb_line() const;
};

/** @brief Represents a single <code>HELIX</code> entry of a PDB file header.
 *
 */
class HelixField : public PdbField {
public:
  core::index2 serial;
  std::string helix_id;  ///< ID of this helix

  /** @name First residue for this helix
   *  The fields provide chain_d as well as id, name and icode of the residue
   */
  //@{
  char chain_from; ///< chain ID
  int residue_id_from; ///< residue ID
  std::string residue_name_from; ///< residue name
  char insert_from;///< insertion code for the first residue
  //@}
  /** @name Last residue for this helix
   *  The fields provide chain_d as well as id, name and icode of the residue
   */
  //@{
  char chain_to; ///< chain ID
  int residue_id_to; ///< residue ID
  std::string residue_name_to; ///< residue name
  char insert_to;///< insertion code for the last residue
  //@}
  char helix_class; ///< Class of this helix
  std::string comment; ///< A comment string
  core::index2 length; ///< The number of residues that belong to this helix

  /// Creates the data structure based on a HELIX line from a PDB file
  HelixField(const std::string &pdb_line);

  virtual std::string to_pdb_line() const;

  /// Returns true if a given residue belongs to this helix
  bool is_my_residue(const core::data::structural::Residue & r);
};

/** @brief Represents a single <code>SHEET</code> entry of a PDB file header.
 *
 */
class SheetField : public PdbField {
public:
  core::index2 strand_id;  ///< ID of this strand
  std::string sheet_id;  ///< ID of the sheet this strand belong to
  core::index2 n_strands;   ///< the number of strands in this sheet

  /** @name first residue for this strand
   *  The fields provide chain_d as well as id, name and icode of the residue
   */
  //@{
  char chain_from; ///< chain ID
  int residue_id_from; ///< residue ID
  std::string residue_name_from; ///< residue name
  char insert_from;///< insertion code for the last residue
  //@}

  /** @name Last residue for this strand
   *  The fields provide chain_d as well as id, name and icode of the residue
   */
  //@{
  char chain_to; ///< chain ID
  int residue_id_to; ///< residue ID
  std::string residue_name_to; ///< residue name
  char insert_to;///< insertion code for the last residue
  //@}

  int sense;

  std::string register_my_atom;
  char register_my_insert;
  char register_my_chain;
  std::string register_my_residue_name;
  int register_my_residue_id;
  std::string register_previous_atom;
  char register_previous_insert;
  char register_previous_chain;
  std::string register_previous_residue_name;
  int register_previous_residue_id;

  /// Creates the data structure based on a SHEET line from a PDB file
  SheetField(const std::string &pdb_line);
  virtual std::string to_pdb_line() const;

  /// Returns true if a given residue belongs to this strand
  bool is_my_residue(const core::data::structural::Residue & r);
};

/** @brief Represents <code>SEQRES</code> field from a PDB file header.
 *
 * \include ex_Seqres.cc
 */
class Seqres : public PdbField {
public:

  /// Map that holds a list of residues for each chain, identified by a single character
  std::unordered_map<char,std::vector<core::chemical::Monomer>> sequences;

  /// Bare constructor
  Seqres() : logger("Seqres") {}

  /// Bare virtual destructor
  inline virtual ~Seqres() {};

  /// Accumulates next <code>SEQRES</code> line from a file header
  void add_line(const std::string & seqres_line);

  /// Writes this <code>SEQRES</code> field in the PDB format
  virtual std::string to_pdb_line() const;

  /// Returns the number of chains found a PDB header
  core::index2 count_chains() { return sequences.size(); }

  /** @brief Creates a Sequence object for a chain.
   *
   * Note, that a corresponding Chain object will not be created.
   * @param chain_id - ID of a requested chain
   * @param sequence_header - a string used to describe this sequence
   */
  core::data::sequence::Sequence_SP create_sequence(const char chain_id, const std::string & sequence_header) {
    return std::make_shared<core::data::sequence::Sequence>(sequence_header.c_str(), sequences[chain_id]);
  }

private:
  utils::Logger logger;
};

/** \brief PDB class reads and writes biomolecular structures in the PDB file format.
 *
 * The following example reads a PDB file (backbone atoms only) and prints some basic statistics:
 * \include ex_Pdb.cc
 */
class Pdb {
public:
  /** @brief Contains objects created by parsing PDB file header.
   *
   * Because there might be several entries of the same type in a header (e.g. "SHEET"), the header data container
   * must be implemented as <code>std::multimap</code> rather than <code>std::map</code>. Here it is how one can
   * iterate over one type of fields:
   * <code><pre>
  std::pair<Pdb::HeaderIterator, Pdb::HeaderIterator> range = reader.header.equal_range("HELIX");
  for (Pdb::HeaderIterator it = range.first; it != range.second; ++it)
    std::cerr << (*it).second->to_pdb_line() << "\n";
</pre></code>
   */
  std::multimap<std::string, std::shared_ptr<PdbField>> header;

  /// For convenience, lets define the iterator type for the header data container.
  typedef typename std::multimap<std::string,std::shared_ptr<PdbField>>::iterator HeaderIterator;

  std::vector<std::shared_ptr<std::vector<Atom>>> atoms;  ///< All the atoms found in the PDB file; not yet assigned to any residue

  /** \brief Constructor reads and parses PDB file but does not create any core::data::structural::Structure object
   *
   * For example, to read only protein backbone and parse the header of a PDB file, use:
   * @code
   * Pdb reader("infile.pdb",is_bb,true); // 'true' says that the header should be parsed
   * @endcode
   * Now we can ask the reader to create a structure (i.e. an object that represents a protein) from the first (indexed by 0) model stored in the PDB file.
   * @code
   * core::data::structural::Structure_SP backbone = reader.create_structure(0);
   * @endcode
   *
   * @param fname - input file name
   * @param predicate - object that says whether a line should be parsed or not; by default all lines are parsed
   * @param if_parse_header - header is not parsed by default; say <code>true</code> to parse it
   * @param first_model_only - load only first model from a file and stop reading
   */
  Pdb(const std::string & fname, const PdbLineFilter & predicate = keep_all,
      const bool if_parse_header = false, const bool first_model_only = false);

  /** \brief Constructor reads and parses PDB data from a stream; it does not create any core::data::structural::Structure object
   *
   * @param infile - input stream with PDB-formatted data
   * @param predicate - object that says whether a line should be parsed or not; by default all lines are parsed
   * @param if_parse_header - header is not parsed by default; say <code>true</code> to parse it
   * @see Pdb(const std::string, const PdbLineFilter &, const bool)
   */
  Pdb(std::istream & infile, const PdbLineFilter & predicate = keep_all) {
    read_pdb(infile, predicate);
  }

  /** @brief Returns the PDB code of a deposit : four-character string.
   *
   * @returns four-character string, e.g. 2AZA, 2GB1, 1HLB etc
   */
  std::string pdb_code() const;

  /// Returns the number of models in the PDB deposit that has been loaded
  core::index2 count_models() const {
    return atoms.size();
  }

  /** @brief Creates a Structure object from a model stored in the PDB data structure
   *
   * @param which_model - model index; starts from 0
   * @returns pointer to the newly created Structure object
   */
  core::data::structural::Structure_SP create_structure(const core::index2 which_model);

  /** @brief Creates all Structure objects that are available in this PDB stream
   *
   * @param structures - vector to hold created structures
   */
  void create_structures(std::vector<core::data::structural::Structure_SP> & structures);

  /** @brief Creates Vec3 object for each valid <code>ATOM</code> line in the input PDB data.
   *
   * Resulting Vec3 objects will be inserted into a given std::vector container.
   *
   * @param infile - stream providing input data
   * @param destination - vector where the data should go; should be large nough to fit the coordinates,
   *	otherwise say if_push_back = true
   * @param if_push_back - when false (which is the default), this method assumes the vector of coordinates is 
   *	long enough to accomodate all atoms, If <code>true</code>, then the given destination vector is emptied
   *	and all the atoms are pushed back to it
   * @param predicate - coordinates will be parsed if and only if a given PDB line satisfies a predicate. By default PdbLineFilter::is_ca predicate is used
   * so only C\f$\alpha\f$ atom are parsed and inserted into the destination vector.
   * @returns the number of atoms created
   */
  static core::index4 read_coordinates(std::istream & infile, std::vector<core::data::basic::Vec3> & destination,
    const bool if_push_back = false, const PdbLineFilter & predicate = is_ca);

  /** @brief Creates Vec3 object for each valid <code>ATOM</code> line in the input PDB data.
   *
   * Resulting Vec3 objects will be inserted into a given std::vector container.
   *
   * @param fname - name of the input PDB file
   * @param destination - vector where the data should go
   * @param if_push_back - when false (which is the default), this method assumes the vector of coordinates is 
   *	long enough to accomodate all atoms, If <code>true</code>, then the given destination vector is emptied
   *	and all the atoms are pushed back to it
   * @param predicate - coordinates will be parsed if and only if a given PDB line satisfies a predicate. By default PdbLineFilter::is_ca predicate is used
   * so only C\f$\alpha\f$ atom are parsed and inserted into the destination vector.
   * @returns the number of atoms created
   */
  static core::index4 read_coordinates(const std::string & fname, std::vector<core::data::basic::Vec3> & destination,
    const bool if_push_back = false, const PdbLineFilter & predicate = is_ca);

private:
  static utils::Logger logger;
  std::string fname_;
  void read_pdb(std::istream & infile, const PdbLineFilter & predicate = keep_all,
                const bool if_parse_header = true, const bool first_model_only = false);
};

/** \brief Tries to find a PDB file in a directory
 *
 * <p>Looks in the specified directory for a file with a given PDB data, identified by
 * a given PDB code. For a given 4-character ID (digit + 3 letters) the method checks
 * the following possibilities:
 *     - given_path/1abc
 *     - given_path/1ABC
 *     - given_path/1abc.pdb
 *     - given_path/1ABC.pdb
 *     - given_path/1ABC.PDB
 *     - given_path/pdb1abc
 *     - given_path/PDB1ABC
 *     - given_path/pdb1abc.ent
 *     - given_path/PDB1ABC.ent
 *     - given_path/pdb1abc.ent.gz
 *     - given_path/PDB1ABC.ent.gz
 *     - given_path/ab/pdb1abc.ent
 *     - given_path/ab/pdb1abc.ent.gz
 * where 1abc and 1ABC denote a lower-case and an upper-case PDB ID, respectively.</p>
 *
 * @param pdbCode four-character PDB ID
 * @param pdbPath directory to look into
 *
 * @return a name of PDB file that was found or an empty string
 */
std::string find_pdb(const std::string & pdbCode, const std::string & pdbPath);

/** @brief Writes a given structure in the PDB format.
 *
 * This method simply calls <code>PdbAtom::to_pdb_line()</code> for each atom fonud in the given structure.
 * @param structure - points to the structure to be written
 * @param out - output stream
 */
void write_pdb(const core::data::structural::Structure_SP structure, std::ostream & out);

/// Returns true if the given string is a valid  PDB code
bool is_pdb_code(const std::string & code);

/** @brief Attempts to extract a PDB code from a file name.
 * This method only strips off commonly used prefixes such as "pdb" and extensions (e.g., .pdb, .pdb.gz, .ent and so on)
 * Actually, it uses the same set of prefixes and suffixes as <code>find_pdb()</code> method
 * It is not guaranteed this method returns a string that would be accepted by <code>is_pdb_code()</code> method
 */
std::string pdb_code_from_file_name(const std::string & code);

} // ~ io
} // ~ data
} // ~ core

#endif

/**
 *  \example ex_Pdb.cc
 *  \example ex_Seqres.cc
 *  \example ex_PdbLineFilter.cc
 */
