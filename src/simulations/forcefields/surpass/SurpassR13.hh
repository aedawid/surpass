#ifndef SIMULATIONS_CARTESIAN_FF_MF_SurpassR13_HH
#define SIMULATIONS_CARTESIAN_FF_MF_SurpassR13_HH

#include <memory>
#include <vector>

#include <core/data/sequence/SecondaryStructure.hh>
#include <core/data/io/ss2_io.hh>

#include <simulations/forcefields/mf/ShortRangeMFBase.hh>
#include <simulations/systems/surpass/SurpassModel.hh>
#include <simulations/representations/surpass_utils.hh>

namespace simulations {
namespace forcefields {
namespace surpass {

using namespace core::calc::numeric;
using namespace core::data::sequence;

/** @brief Knowledge based potential that evaluates energy of a \f$R_{13}\f$ distance.
 */
template<typename C>
class SurpassR13 : public mf::ShortRangeMFBase<C> {
public:

  /** \brief Creates SurpassR13 energy function based on string parameters.
   *
   * @param system - the system whose energy will be evaluated
   * @param scored_sequence - secondary structure object provides the sequence and secondary structure data for energy evaluation
   * @param parameters - a vector of string parameters. Three  parameters are required by this constructor:
   *     - energy parameters file; use '-' (dash) character to use the default file that is stored in parameters set
   *     - pseudocounts, e.g "0.001"
   */
  SurpassR13(const systems::surpass::SurpassModel <C> &system, const SecondaryStructure_SP scored_sequence, const std::vector<std::string> &parameters) :
    SurpassR13<C>(system, scored_sequence,parameters[0], utils::from_string<core::real>(parameters[1])) { }

  /** \brief Creates SurpassR12 energy function.
   *
   * @param system - the system whose energy will be evaluated
   * @param scored_sequence - secondary structure object provides the sequence and secondary structure data for energy evaluation
   * @param ff_file - energy parameters file; use '-' (dash) character to use the default file that is stored in parameters set
   * @param pseudocounts - pseudocounts value, e.g 0.001
   */
  SurpassR13(const systems::surpass::SurpassModel <C> &system, const SecondaryStructure_SP scored_sequence, const std::string &ff_file,
             const double pseudocounts = -1) : mf::ShortRangeMFBase<C>(system, (ff_file != "-") ? ff_file : "forcefield/local/R13_surpass.dat",
      simulations::representations::surpass_representation(*scored_sequence), 0, 2, 3, pseudocounts),
    xyz(system.coordinates), logger("SurpassR13") { }

  /// Empty virtual constructor to satisfy the compiler
  virtual ~SurpassR13() {}

  /// Returns the name of this energy term, which is "SurpassR13"
  const std::string & name() const { return name_; }

  /** \brief Calculates energy of pseudo-bonds connected to a given residue.
   *
   * The energy is evaluated based on the structural property (\f$R_{13}\f$ distance) and the mean field potential functions.
   */
  virtual inline double calculate_by_residue(const residue_index which_residue) {

    core::real en = 0.0;
    if ((which_residue >= 2) && (which_residue < mf::ShortRangeMFBase<C>::n_residues)) { // distance towards the N termini works for residue at least at index 4
      core::real d = xyz[which_residue].distance_to(xyz[which_residue - 2]);
      en += mf::ShortRangeMFBase<C>::score_property(which_residue - 2, d);
    }
    if (which_residue <= mf::ShortRangeMFBase<C>::last_positions_scored) {
      core::real d = xyz[which_residue].distance_to(xyz[which_residue + 2]);
      en += mf::ShortRangeMFBase<C>::score_property(which_residue, d);
    }

    return en;
  }

private:
  const std::unique_ptr<C[]> &xyz;
  static const std::string name_;
  utils::Logger logger;
};

template<typename C>
const std::string SurpassR13<C>::name_ = "SurpassR13";

}
}
}

#endif
