#ifndef SIMULATIONS_GENERIC_FF_CalculateEnergyBase_HH
#define SIMULATIONS_GENERIC_FF_CalculateEnergyBase_HH

#include <string>
#include <core/real.hh>
#include <core/data/basic/Array2D.hh>

namespace simulations {
namespace forcefields {

/// Bare interface for a class that can calculate energy.
class CalculateEnergyBase {
public:
  /// calculates energy
  virtual double calculate() = 0;

  /// Virtual destructor
  virtual ~CalculateEnergyBase() {}

  /// Should return the name of a derived energy term
  virtual const std::string & name() const = 0;
};

/// Defines a new type that represents a shared pointer to a CalculateEnergyBase object
typedef std::shared_ptr<CalculateEnergyBase> CalculateEnergyBase_SP;

} // ~ simulations
} // ~ ff
#endif
