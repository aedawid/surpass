/** \file sequence_utils.hh
 * @brief Some utility methods that operate on protein and nucleic sequences.
 *
 */
#ifndef CORE_DATA_SEQUENCE_sequence_utils_H
#define CORE_DATA_SEQUENCE_sequence_utils_H

#include <string>
#include <cstdlib>
#include <regex>

#include <core/index.hh>
#include <core/data/sequence/Sequence.hh>

namespace core {
namespace data {
namespace sequence {

/// Extracts gi-number form a sequence header string.
/**
 * @param sequence_header - a sequence header which contains a gi number, formatted as gi|123456|
 * @return gi number if it was found in the given sequence header; 0 otherwise
 */
core::index4 extract_gid(const std::string & sequence_header);

/** @brief Returns true if a given string contains <strong>only</strong> characters corresponding to amino acids (one-letter code)
 * @param sequence - a putative sequence
 * @return true if this is an amino acid sequence
 */
bool contains_aa_only(const std::string &sequence);

/** @brief Returns the number of residues in a sequence that are not gaps.
 * @param sequence - an input sequence
 * @return the number of non-gap symbols in that sequence
 */
core::index2 length_without_gaps(const std::string & sequence);

/** @brief Removes all gap characters ('_' and '-') from a given sequence string.
 * @param sequence_string - an input sequence
 */
void remove_gaps(std::string & sequence_string);

/** @brief Removes all gap characters ('_' and '-') from a given sequence string.
 * @param sequence_string - an input sequence
 */
std::string remove_gaps(const std::string & sequence_string);

/** @brief Counts amino acid characters in the given sequence
 * @param seq - an input sequence
 * @return the number of characters corresponding to amino acids
 */
core::index2 count_aa(const std::string &seq);

/** @brief Counts amino acid residues in the given sequence
 * @param seq - an input sequence
 * @return the number of residues recognized as amino acid monomers
 */
core::index2 count_aa(const Sequence &seq);

/** @brief Counts nucleic acid residues in the given sequence
 * @param seq - an input sequence
 * @return the number of characters corresponding to nucleic acids (A, C, T, G, U, a, c, t, g or u)
 */
core::index2 count_nt(const std::string &seq);

/** @brief Counts amino acid residues in the given sequence
 * @param seq - an input sequence
 * @return the number of residues recognized as nucleic acid monomers
 */
core::index2 count_nt(const Sequence &seq);

/** @brief Counts how many times each amino acid type appears in a given sequence
 * @param seq - an input sequence
 * @param counts - vector of size 21 to store the statistics. This method does not clear the vector, but if its shorter than 21,
 *      it will be resized
 * @return the vector with (updated) statistics
 */
std::vector<core::index4> & count_each_aa(const Sequence &seq, std::vector<core::index4> & counts);

/** @brief Counts how many times each amino acid type appears in a given sequence
 * @param seq - an input sequence
 * @param counts - vector of size 21 to store the statistics. This method does not clear the vector, but if its shorter than 21,
 *      it will be resized
 * @return the vector with (updated) statistics
 */
std::vector<core::index4> & count_each_aa(const std::string &seq, std::vector<core::index4> & counts);

/** @brief Counts how many different amino acid types show up in a given amino acid sequence.
 * @param seq - an input sequence
 * @return the number of different amino acid types, e.g. for the sequence "GPGPGPGPGP" this method will return 2
 */
core::index1 count_aa_types(const std::string &seq);

/** @brief Counts how many different amino acid types show up in a given amino acid sequence.
 * @param seq - an input sequence
 * @return the number of different amino acid types, e.g. for the sequence "GPGPGPGPGP" this method will return 2
 */
core::index1 count_aa_types(const Sequence &seq) ;

}
}
}

#endif
