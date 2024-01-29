/** @brief Provides methods to read and write monomer definitions from files in different formats
 *
 */
#include <string>

namespace core {
namespace chemical {

/** @brief Read monomers from a CIF file
 * @param cif_filename - name of the input file
 */
void read_monomers_cif(const std::string &cif_filename);

/** @brief Read monomers from a file in a binary form (internal format)
 * @param bin_filename - name of the input file
 */
void read_monomers_binary(const std::string &bin_filename);

/** @brief Stores monomers in in a binary file
 * @param bin_filename - name of the output file
 */
void write_monomers_binary(const std::string &bin_filename);

/** @brief Read monomers from a flat-text file (internal format)
 * @param txt_filename - name of the input file
 */
void read_monomers_txt(const std::string &txt_filename);

/** @brief Stores monomers in in a flat-text file
 * @param txt_filename - name of the output file
 */
void write_monomers_txt(const std::string &txt_filename);

/** @brief Loads monomers from the database.
 *
 * The standard monomers (the 20 amino acids and the five nucleotides) are defined as static objects.
 * All the others are loaded from a file when needed. Call this method to trigger monomers loading.
 */
void load_monomers_from_db();

}
}
