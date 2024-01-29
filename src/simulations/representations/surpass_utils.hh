/** \file surpass_utils.hh
 * @brief Provide utility methods for SURPASS model
 */
#ifndef SIMULATIONS_REPRESENTATIONS_surpass_utils_HH
#define SIMULATIONS_REPRESENTATIONS_surpass_utils_HH

#include <core/data/structural/Structure.hh>

namespace simulations {
namespace representations {

/** @brief Assigns secondary structure for all residues in the given structure according to Surpass atom names
 * @param strctr - structure to be updated
 */
core::data::structural::Structure & fix_surpass_ss_assignment(core::data::structural::Structure &surpass_strctr);

/** @brief Converts an all-atom structure to SURPASS representation.
 *
 * The simple program shows how to do the conversion:
 * \include ex_surpass_representation.cc
 * @param strctr - a full atom structure (at least all the backbone atoms must be present)
 * @return a structure in the SURPASS representation
 */
core::data::structural::Structure_SP surpass_representation(const core::data::structural::Structure & strctr);

/** @brief Converts a SecondaryStructure object to SURPASS representation.
 *
 * The features of the resulting object:
 *   - the amino acid sequence will be three residues shorter
 *   - the sequence is changed, e.g. to all - GLY
 *   - secondary structure is modified according to SURPASS conversion rules
 * @param secstr input secondary structure for a regular protein sequence
 * @param check_if_converted - if set to true (the default behavior), the method will check if the sequence is all-GLY.
 * If so, it assumes the input object has already been converted to SURPASS representation and no conversion is required.
 * @return a secondary structure in the SURPASS representation
 */
core::data::sequence::SecondaryStructure_SP surpass_representation(const core::data::sequence::SecondaryStructure &secstr,
  const bool check_if_converted = true);

/** @brief Returns true if a given structure is a correct model in SURPASS representation
 * @param strctr - an input structure
 * @return true if the argument in a SURPASS model
 */
bool is_surpass_model(const core::data::structural::Structure & strctr);

}
}

#endif

/**
 * \example ex_surpass_representation.cc
 */
