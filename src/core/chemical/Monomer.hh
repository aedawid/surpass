#ifndef CORE_CHEMICAL_Monomer_HH
#define CORE_CHEMICAL_Monomer_HH

#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <stdexcept>

#include <core/chemical/monomer_io.hh>
#include <core/index.hh>
#include <core/real.hh>
#include <utils/Logger.hh>

namespace core {
namespace chemical {

/** @brief Represents a monomer type that may be found in a PDB file, e.g. an amino acid residue type or a ligand molecule
 *
 * It is possible to iterate over all monomers, that are defined in the library. One can also
 * access monomer by its three-letter code, one-letter code or by an internal integer index. Finally,
 * the standard monomers may be accessed directly as variables because they are declared as static objects.
 * The features were shown in the example below:
 *
 * \include ex_Monomer.cc
 */
class Monomer {
public:
  core::index2 id;        ///< Unique integer identifier for this monomer type
  core::index2 parent_id; ///< For derived monomers : ID of the parent monomer type; otherwise equal to the <code>id</code> field
  char code1;             ///< One-letter code to represent this monomer, e.g. 'A' for alanine

  /** \brief Define the type of this monomer
   *
   *     - P - protein residue (amino acid)
   *     - S - sugar
   *     - N - nucleic residue
   *     - U - unknown (other than the two above)
   */
  char type;
  unsigned char n_atoms;  ///< the number of atoms in this monomer
  unsigned char n_heavy_atoms;  ///< the number of non-hydrogen atoms in this monomer
	bool is_ambiguous;    ///< ambiguity flag; when true, this particular three-letter code represents more than one chemical entity
	core::real charge;  ///< formal charge of the monomer
	std::string code3;  ///< Three letter code of this monomer, e.g. "ALA" for alanine

  /** @brief Constructor initializes this object with provided data
   *
   * @param id - index of this monomer; counting starts from 0
   * @param code1 - one letter code for this monomer; e.g. 'P' for proline
   * @param code2 - three-letter code for this monomer; e.g. 'PRO' for proline
   * @param type - either 'P', 'N' or 'U', see above
   * @param n_atoms - the number of non-leaving atoms in this monomer
   * @param n_heavy_atoms - the number of non-hydrogen atoms in this monomer
   * @param ambig_flag - whether the monomer code3 is ambiguous or not
   * @param charge - formal charge of the monomer
   * @param parent - id of the parent monomer
   */
  Monomer(const core::index2 id, const char code1, const std::string code3, const char type,
		  const unsigned char n_atoms, const unsigned char n_heavy_atoms,
		  const bool ambig_flag, const core::real charge, const core::index2 parent);

  /** @brief Constructor initializes this object with provided data
   *
   * @param id - index of this monomer; counting starts from 0
   * @param code1 - one letter code for this monomer; e.g. 'P' for proline
   * @param code2 - three-letter code for this monomer; e.g. 'PRO' for proline
   * @param type - either 'P', 'N' or 'U', see above
   * @param n_atoms - the number of non-leaving atoms in this monomer
   * @param n_heavy_atoms - the number of non-hydrogen atoms in this monomer
   * @param ambig_flag - whether the monomer code3 is ambiguous or not
   * @param charge - formal charge of the monomer
   * @param parent - id of the parent monomer
   */
  Monomer(const core::index2 id,const char code1,const char* code3,const char type,const unsigned char n_atoms,const unsigned char n_heavy_atoms,
      const bool ambig_flag,const core::real charge,const core::index2 parent);

  /// Copying constructor
  Monomer(const Monomer &m);

  /// Constructor creates an object from a string in an internal format
  Monomer(const std::string & line) ;

  /// Returns true if this is a standard monomer
  bool is_standard() const { return id<1000; }

  /** @brief Static method to register a new monomer in the library.
   *
   * @param m - the new monomer to be inserted a the monomers' set
   * @return true if the monomer was actually inserted; false if it was already present in the set
   */
  static bool register_monomer(Monomer & m) {

    if(is_known_monomer(m.code3)) return false;
    by_code3.insert(std::pair<std::string,Monomer>(m.code3,m));
    by_code1.insert(std::pair<char,Monomer>(m.code1,m));
    by_id().push_back(m);
    logs << utils::LogLevel::FINEST<<"Registering a new monomer with code3 >"<<m.code3<<"< and code1: "<<char(m.code1)<<"\n";
    return true;
  }

  /** @brief Returns a monomer for a given single-character code.
   *
   * Note, that only a standard monomer (i.e. one of the 20 amino acids or one of the 5 nucleotides) may be accessed by this method
   * @param code1 - monomer's code, e.g. 'A' for alanine
   * @return the monomer requested
   */
	static const Monomer & get(const char code1) { return by_code1.at(code1); }

  /** @brief Returns a monomer for a given three-character code.
   *
   * @param code3 - monomer's code, e.g. 'ALA' for alanine
   * @return the monomer requested
   */
	static const Monomer & get(const std::string & code3) {
#ifdef DEBUG
    if(by_code3.find(code3)==by_code3.end())
      logs <<utils::LogLevel::WARNING<<"Can't find a monomer type for code: "<<code3<<", loading the DB file\n";
#endif
    try {
		return by_code3.at(code3);
    } catch(std::out_of_range & e) {
      logs << utils::LogLevel::FINE<<"extended monomer requested: "<<code3<<", loading the database file\n";
      load_monomers_from_db();
#ifdef DEBUG
      if(by_code3.find(code3)==by_code3.end()) {
        logs <<utils::LogLevel::SEVERE<<"Can't find a monomer type for code: "<<code3<<", even though the full database was loaded\n"
        <<"\tUpdate the database file or include the monomer information in the PDB header!\n";
      }
#endif
      return by_code3.at(code3);
    }
	}

  /** @brief Returns a monomer for a given integer id
   *
   * @param id - monomer's index, according to the internal indexing
   * @return the monomer requested
   */
	static const Monomer & get(const core::index2 id) {
		return by_id()[id];
	}

  /** @brief Check whether there is a monomer registered under a given three-letter code
   *
   * @param code3 - monomer's code, e.g. 'ALA' for alanine
   * @return true, if the monomer has been added to the library
   */
	static bool is_known_monomer(const std::string &code3) {
		if (by_code3.find(code3) == by_code3.end()) return false;
		return true;
	}

  /** @brief begin constant iterator for monomers.
   *
   * @return constant iterator pointing to the very first monomer
   */
  static std::vector<Monomer>::const_iterator cbegin() { return by_id().cbegin(); }

  /** @brief end constant iterator for monomers
   *
   * @return constant iterator pointing to the very first monomer
   */
  static std::vector<Monomer>::const_iterator cend() { return by_id().cend(); }

  /** @brief begin constant iterator for standard amino acids
   *
   * @return constant iterator pointing to the very first standard amino acid monomer
   */
  static std::vector<Monomer>::const_iterator aa_cbegin() { return by_id().cbegin(); }

  /** @brief end constant iterator for standard amino acids
   *
   * @return constant iterator pointing to the very last standard amino acid monomer
   */
  static std::vector<Monomer>::const_iterator aa_cend() { return by_id().cbegin() + 20; }

  /** @name Standard monomers : amino acids and nucleotides
   *  Standard monomers are declared as static object so they are readily available. All other
   *  monomers are loaded from a file when needed.
   */
  const static Monomer ALA;
	const static Monomer ARG;
	const static Monomer ASN;
	const static Monomer ASP;
	const static Monomer CYS;
	const static Monomer GLN;
	const static Monomer GLU;
	const static Monomer GLY;
	const static Monomer HIS;
	const static Monomer ILE;
	const static Monomer LEU;
	const static Monomer LYS;
	const static Monomer MET;
	const static Monomer PHE;
	const static Monomer PRO;
	const static Monomer SER;
	const static Monomer THR;
	const static Monomer TRP;
	const static Monomer TYR;
	const static Monomer VAL;
	const static Monomer UNK;
	const static Monomer a;
	const static Monomer c;
	const static Monomer g;
	const static Monomer t;
	const static Monomer u;
	const static Monomer GAP;
    const static Monomer GPE;
	const static Monomer UNL;
  	const static Monomer UNG;
  //@}

  friend std::ostream &operator<<(std::ostream &out, const Monomer &m);

  friend std::ostream &operator<<(std::ostream &out, const Monomer *m);

  friend void read_monomers_cif(const std::string &cif_filename);

  bool operator==(const Monomer &m) const { return id == m.id; }

  /** @brief Returns <code>true</code> if <code>id</code> of this monomer is lower than <code>m.id</code>
   * @param m - another monomer
   */
  bool operator<(const Monomer &m) const { return id < m.id; }

private:
  static std::unordered_map<std::string, Monomer> by_code3;
  static std::unordered_map<char, Monomer> by_code1;
  static std::vector<Monomer> & by_id();
  static utils::Logger logs;

  static std::unordered_map<std::string, Monomer> create_map3();
  static std::unordered_map<char, Monomer> create_map1();
};

std::ostream& operator<<(std::ostream &out, const Monomer &m);
std::ostream& operator<<(std::ostream &out, const Monomer *m);
}
}

#endif

/**
 *\example ex_Monomer.cc
 */
