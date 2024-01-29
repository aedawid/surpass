#ifndef SIMULATIONS_FORCEFIELDS_surpass_SurpassA13_HH
#define SIMULATIONS_FORCEFIELDS_surpass_SurpassA13_HH

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

#include <core/calc/numeric/interpolators.hh>
#include <core/calc/numeric/Interpolate1D.hh>
#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/io/ss2_io.hh>
#include <utils/string_utils.hh>
#include <core/calc/structural/angles.hh>

#include <simulations/forcefields/mf/ShortRangeMFBase.hh>
#include <simulations/systems/surpass/SurpassModel.hh>
#include <simulations/representations/surpass_utils.hh>
#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/systems/ResidueChain.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

using namespace core::calc::numeric;
using namespace core::data::sequence;

/** @brief Knowledge based potential that evaluates energy of a \f$A_{13}\f$ planar angle.
 */
template<typename C>
class SurpassA13 : public mf::ShortRangeMFBase<C> {
public:

  /** \brief Creates SurpassR13 energy function based on string parameters.
   *
   * @param system - the system whose energy will be evaluated
   * @param scored_sequence - secondary structure object provides the sequence and secondary structure data for energy evaluation
   * @param parameters - a vector of string parameters. Three  parameters are required by this constructor:
   *     - energy parameters file; use '-' (dash) character to use the default file that is stored in parameters set
   *     - pseudocounts, e.g "0.001"
   */
  SurpassA13(const systems::surpass::SurpassModel <C> &system, const SecondaryStructure_SP scored_sequence, const std::vector<std::string> &parameters) :
    SurpassA13<C>(system, scored_sequence,parameters[0], utils::from_string<core::real>(parameters[1])) { }

  /** \brief Creates SurpassA12 energy function.
   *
   * @param system - the system whose energy will be evaluated
   * @param scored_sequence - secondary structure object provides the sequence and secondary structure data for energy evaluation
   * @param ff_file - energy parameters file; use '-' (dash) character to use the default file that is stored in parameters set
   * @param pseudocounts - pseudocounts value, e.g 0.001
   */
  SurpassA13(const systems::surpass::SurpassModel <C> &system, const SecondaryStructure_SP scored_sequence, const std::string &ff_file,
             const double pseudocounts = -1) : mf::ShortRangeMFBase<C>(system, (ff_file != "-") ? ff_file : "forcefield/local/A13_surpass.dat",
      simulations::representations::surpass_representation(*scored_sequence), 0, 2, 3, pseudocounts),
    xyz(system.coordinates), logger("SurpassA13") { }

  /// Empty virtual constructor to satisfy the compiler
  virtual ~SurpassA13() {}

  /// Returns the name of this energy term, which is "SurpassA13"
  const std::string & name() const { return name_; }

  /** \brief Calculates energy related to a given residue.
   *
   * The energy is evaluated based on the structural property (\f$R_{15}\f$ distance) and the mean field potential functions.
   */
  virtual inline double calculate_by_residue(const core::index2 which_residue)   {

    core::real en = 0.0;

    if ((which_residue >= 2) && (which_residue < mf::ShortRangeMFBase<C>::n_residues)) { // distance towards the N termini works for residue at least at index 4
      core::real a = core::calc::structural::evaluate_planar_angle(
        xyz[which_residue],
        xyz[which_residue - 1],
        xyz[which_residue - 2]);
      en += mf::ShortRangeMFBase<C>::score_property(which_residue - 2, a);
//      std::cerr<<"-2: "<<which_residue<<" "<<a<<" "<<en<<"\n";/////////////////////////////////////////////////////////
    }
    if (which_residue <= mf::ShortRangeMFBase<C>::last_positions_scored) {
      core::real a = core::calc::structural::evaluate_planar_angle(
        xyz[which_residue],
        xyz[which_residue + 1],
        xyz[which_residue + 2]);
      en += mf::ShortRangeMFBase<C>::score_property(which_residue, a);
//      std::cerr<<"+2: "<<which_residue<<" "<<a<<" "<<en<<"\n";/////////////////////////////////////////////////////////
    }

    return en;
  }

private:
  static const std::string name_;
  const std::unique_ptr<C[]> &xyz;
  utils::Logger logger;
};

template<typename C>
const std::string SurpassA13<C>::name_ = "SurpassA13";

}
}
}

#endif
