/** @file TotalEnergyByResidue.hh
 * @brief Provides TotalEnergyByResidue type and related ostream operator
 */
#ifndef SIMULATIONS_GENERIC_FF_TotalEnergyByResidue_HH
#define SIMULATIONS_GENERIC_FF_TotalEnergyByResidue_HH

#include <core/index.hh>
#include <core/data/basic/Array2D.hh>

#include <utils/Logger.hh>

#include <simulations/forcefields/TotalEnergy.hh>
#include <simulations/forcefields/ByResidueEnergy.hh>
#include <simulations/evaluators/Evaluator.hh>

namespace simulations {
namespace forcefields {

/** @brief Total energy based on <code>ByResidueEnergy</code> terms.
 * This class itself also implements ByResidueEnergy interface, so one can calculate
 * weighted combination of energy terms for a given residue or a chunk of a modelled chain
 */
class TotalEnergyByResidue : public TotalEnergy<ByResidueEnergy>, public virtual ByResidueEnergy {
public:

  /// Empty c'tor
  TotalEnergyByResidue() : logger("TotalEnergyByResidue") {}

  /// Virtual destructor
  virtual ~TotalEnergyByResidue() {}

  /// Calculates the total energy of a system
  virtual double calculate() { return TotalEnergy<ByResidueEnergy>::calculate(); }

  const std::string & name() const { return name_; }

  /// The returned string is the header line for energy components table (printed by ostream operator)
  virtual std::string header_string() const;

  virtual double calculate_by_residue(const core::index2 which_residue);

  virtual double calculate_by_chunk(const core::index2 chunk_from, const core::index2 chunk_to);

private:
  utils::Logger logger;
  static const std::string name_;

  friend std::ostream &operator<<(std::ostream &out, const TotalEnergyByResidue &e);
};

/** @brief Prints a nice table with energy broken out into components
 * @param out - where to print the row of energy values
 * @param e - energy to be printed
 */
std::ostream &operator<<(std::ostream &out, const TotalEnergyByResidue &e);

} // ~ simulations
} // ~ ff
#endif
