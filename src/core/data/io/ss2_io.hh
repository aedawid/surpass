/** \file ss2_io.hh
 * @brief Provides various I/O methods for reading in secondary structure profiles in the SS2 (i.e. PsiPred) file format
 */
#ifndef CORE_DATA_IO_ss2_io_H
#define CORE_DATA_IO_ss2_io_H

#include <string>
#include <memory>
#include <iostream>
#include <vector>
#include <utils/io_utils.hh>
#include <core/data/sequence/Sequence.hh>
#include <core/data/sequence/SecondaryStructure.hh>

namespace core {
namespace data {
namespace io {

/** @brief Reads in secondary structure probability table in PsiPred's SS2 format.
 *
 * @param in_stream - input stream
 * @param header - header string of the sequence that will be created. The header will be used e.g. when the sequence
 * is printed in the FASTA format
 * @return newly created SecondaryStructure instance, containing the data obtained from the given stream
 */
core::data::sequence::SecondaryStructure_SP read_ss2(std::istream & in_stream, const std::string & header = "");

/** @brief Reads in secondary structure probability table in PsiPred's SS2 format.
 *
 * @param in_fname - input file name
 * @param header - header string of the sequence that will be created. The header will be used e.g. when the sequence
 * is printed in the FASTA format
 * @return newly created SecondaryStructure instance, containing the data obtained from the given file
 */
core::data::sequence::SecondaryStructure_SP read_ss2(const std::string & in_fname, const std::string & header = "");

/** @brief Reads in secondary structure probability table in PsiPred's SS2 format.
 *
 * The created object, which is derived from Sequence class, will be inserted to the given vector.
 * @param in_fname - input file name
 * @param sink - where to store the newly created sequence (it will be appended at the end of the vector)
 * @param header - header string of the sequence that will be created. The header will be used e.g. when the sequence
 * is printed in the FASTA format
 * @return a reference to the <code>sink</code> vector
 */
std::vector<std::shared_ptr<core::data::sequence::Sequence>> & read_ss2(const std::string & in_fname,
    std::vector<std::shared_ptr<core::data::sequence::Sequence>> & sink, const std::string & header = "");

/** @brief Writes a SecondaryStructure object into the SS2 (PsiPred) data format
 * @param ss - object to be written
 * @param out - output stream
 *
 * The following example reads a DSSP file and writes the data in SS@ format:
 * @include ex_dssp_to_ss2.cc
 */
void write_ss2(const core::data::sequence::SecondaryStructure & ss, std::ostream & out);

}
}
}

#endif
