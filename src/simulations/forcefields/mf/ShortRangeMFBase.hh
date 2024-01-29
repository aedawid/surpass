#ifndef SIMULATIONS_CARTESIAN_FF_MF_ShortRangeMFBase_HH
#define SIMULATIONS_CARTESIAN_FF_MF_ShortRangeMFBase_HH

#include <cmath>

#include <iostream>
#include <memory>
#include <vector>
#include <map>

#include <core/index.hh>
#include <core/calc/numeric/interpolators.hh>
#include <core/calc/numeric/Interpolate1D.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/sequence/SecondaryStructureAnnotation.hh>
#include <utils/string_utils.hh>

#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/systems/ResidueChain.hh>
#include <simulations/forcefields/ShortRangeEnergyBase.hh>
#include <simulations/forcefields/mf/MeanFieldDistributions.hh>

namespace simulations {
namespace forcefields {
namespace mf {

using namespace core::calc::numeric;
using core::data::sequence::SecondaryStructure;
using core::data::sequence::SecondaryStructure_SP;

/** @brief A base class for secondary structure dependent mean-field short range energy types.
 *
 * This is CABS-like local (along amino acid sequence) energy function that depends both
 * on an amino acid sequence and a secondary structure.
 */
template<typename C>
class ShortRangeMFBase : public ShortRangeEnergyBase<C> {
public:

  const signed char first_aa_pos;
  const signed char second_aa_pos;
  const core::index2 n_residues;

  /** @brief Constructor loads energy terms (spline function parameters) and repacks them into vectors indexes by residue index.
   *
   * The re-packing of these parameters takes care both of the sequence and secondary structure of the scored protein chain.
   * @param system - the system whose energy will be evaluated
   * @param ff_file - config file used by MeanFieldDistribution1D base class to load the distributions
   * @param scored_secondary - (a pointer to) the secondary structure of the scored chain
   * @param first_aa_pos - relative index of the first residue the scoring function depends on
   * @param secnd_aa_pos - relative index of the first residue the scoring function depends on
   * @param property_span - the number of residues involved in a single structural property measurement, e.g. three in the case of planar angle or 5 for \f$R_{15}\f$ distance
   * @param pseudocounts - pseudocounts value used by MeanFieldDistribution1D base class to convert probabilities into energy
   *
   * @see MeanFieldDistributions
   */
  ShortRangeMFBase(const systems::ResidueChain<C> &system, const std::string &ff_file,
                   const SecondaryStructure_SP &scored_secondary, const signed char first_aa_pos,
                   const signed char secnd_aa_pos, const unsigned char property_span, const double pseudocounts = -1) :
     ShortRangeEnergyBase<C>(system, property_span), first_aa_pos(first_aa_pos), second_aa_pos(secnd_aa_pos),
     n_residues(system.count_residues()), scored_secondary(scored_secondary), ff_for_sequence_(9) {

    prepare_ff(ff_file,pseudocounts);
  }

  /// Bare virtual destructor to obey the rules
  virtual ~ShortRangeMFBase() {}

  /// Returns the secondary structure object used for the energy calculations
  const SecondaryStructure_SP sequence() const { return scored_secondary; }

  core::real score_property(core::index2 which_residue,core::real value) const {

    const core::data::sequence::HecFractions ss_first = sequence()->fractions(which_residue + first_aa_pos);
    const core::data::sequence::HecFractions ss_secnd = sequence()->fractions(which_residue + second_aa_pos);

    core::real en = 0;
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        en += (*ff_for_sequence_[ss_to_index(j, k)][which_residue])(value) * ss_first[j] * ss_secnd[k];
      }
    }

    return en;
  }

protected:
  const core::data::sequence::SecondaryStructure_SP scored_secondary;
  std::vector<std::vector<EnergyComponent_SP>> ff_for_sequence_; // 9 elements for each HEC combination

  /** \brief Prepares local energy functions for each position in the sequence.
   *
   * To each residue in the simulated chain, exactly three functions: for H, C and E secondary structure type
   * are assigned to score the local property such as a distance or an angle by a derived class. The property
   * scoring function assigned to the index <code>i</code> will depend on amino acid types at positions
   * <code>i+first_aa_pos</code> and <code>i+secnd_aa_pos</code>.
   */
  void prepare_ff(const std::string & ff_file, const core::real pseudocounts) {

    std::shared_ptr<MeanFieldDistributions> mf = load_1D_distributions(ff_file, pseudocounts);

    for (int i = 0; i < 9; ++i) ff_for_sequence_[i].resize(ShortRangeEnergyBase<C>::last_positions_scored + 1);
    const char codes[] = {'H', 'E', 'C'};

    std::string key("__.__");
    for (int i = 0; i <= ShortRangeEnergyBase<C>::last_positions_scored; ++i) {
      key[0] = char(scored_secondary->sequence[i + first_aa_pos]);
      key[1] = char(scored_secondary->sequence[i + second_aa_pos]);
      for(int j=0;j<3;++j) {
        key[3] = codes[j];
        for(int k=0;k<3;++k) {
          key[4] = codes[k];
          if (!mf->contains_distribution(key))
            throw std::runtime_error(missing_error_msg(key,mf));
          ff_for_sequence_[ss_to_index(j,k)][i] = mf->at(key);
        }
      }
    }
  }

private:

  inline core::index1 ss_to_index(const core::index1 first_ss,const core::index1 second_ss) const { return first_ss * 3 + second_ss; }

  std::string missing_error_msg(const std::string &key, const std::shared_ptr<MeanFieldDistributions> mf) const {

    const std::vector<std::string> valid_keys = mf->known_distributions();
    std::stringstream o;
    o << "Can't find the distribution: " << key << ". Distributions known for this ShortRangeEnergy instance :\n";
    std::copy(valid_keys.cbegin(), valid_keys.cend(), std::ostream_iterator<std::string>(o," "));
    return o.str();
  }
};

}
}
}

#endif

