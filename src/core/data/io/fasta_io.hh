/** \file fasta_io.hh
 * @brief Provides various I/O methods for FASTA file format.
 *
 * These methods can be used to write sequences aligned with each other but do not write alignment scores etc.
 */
#ifndef CORE_DATA_IO_fasta_io_H
#define CORE_DATA_IO_fasta_io_H

#include <string>
#include <iostream>
#include <vector>

#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <core/alignment/PairwiseAlignment.hh>

namespace core {
namespace data {
namespace io {

using core::data::sequence::Sequence;
using core::data::sequence::SecondaryStructure;

/** @brief Reads a file in the FASTA format.
 *
 * The method creates Sequence objects from the data and places them on a given vector
 *
 * @include ex_split_fasta.cc
 * @include ex_find_in_fasta.cc
 *
 * @param file_name - name of the input file
 * @param sink - where to store the newly created sequences
 * @return a reference to the <code>sink</code> vector
 */
std::vector<std::shared_ptr<Sequence>> & read_fasta_file(const std::string file_name,
    std::vector<std::shared_ptr<Sequence>> & sink);

/** @brief Converts a given header-sequence pair of strings into the FASTA format.
 *
 * @param header - sequence name, e.g.
 * @param sequence - the sequence itself
 * @param line_width - the length of each sequence line (the sequence string will be broken to match the width)
 * @param include_header - if true (the default behavior), FASTA header is included as the first line of the returned string
 * @return a string in FASTA format
 */
std::string create_fasta_string(const std::string & header, const std::string & sequence,
    const core::index2 line_width = 65535, const bool include_header = true);

/** @brief Converts a given Sequence object into the FASTA-formatted string.
 *
 * @param seq - the sequence to be converted
 * @param line_width - the length of each sequence line (the sequence string will be broken to match the width)
 * @param include_header - if true (the default behavior), FASTA header is included as the first line of the returned string
 * @return a string in FASTA format
 */
std::string create_fasta_string(const Sequence & seq, const core::index2 line_width = 65535, const bool include_header = true);

/** @brief Converts a given secondary structure stored in a SecondaryStructure object into a FASTA-formatted string.
 *
 * @param sec_str - the secondary structure to be converted
 * @param line_width - the length of each sequence line (the string will be broken to match the width)
 * @param include_header - if true (the default behavior), FASTA header is included as the first line of the returned string
 * @return a string in FASTA format
 */
std::string create_fasta_secondary_string(const SecondaryStructure & sec_str, const core::index2 line_width = 65535, const bool include_header = true);

/** @brief Converts a given sequence alignment into a FASTA-formatted string.
 *
 * @param sec_str - the secondary structure to be converted
 * @param line_width - the length of each sequence line (the string will be broken to match the width)
 * @param include_header - if true (the default behavior), FASTA header is included as the first line of the returned string
 * @return a string in FASTA format
 */
std::string create_fasta_string(const core::alignment::PairwiseAlignment &ali, const Sequence &query_sequence,
                                const Sequence &tmplt_sequence, const core::index2 line_width = 65535);

std::istream & operator>>(std::istream & cin, core::data::sequence::Sequence_SP & seq);

}
}
}

/**
 * @example ex_split_fasta.cc
 * @example ex_find_in_fasta.cc
 */
#endif
