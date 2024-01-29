#ifndef SIMULATIONS_GENERIC_EVALUATORS_CrmsdEvaluator_HH
#define SIMULATIONS_GENERIC_EVALUATORS_CrmsdEvaluator_HH

#include <memory>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <core/data/structural/Structure.hh>
#include <core/calc/structural/transformations/Crmsd.hh>

#include <simulations/systems/CartesianAtomsSimple.hh>
#include <simulations/evaluators/Evaluator.hh>

namespace simulations {
namespace evaluators {
namespace cartesian {

using namespace std::chrono;

template <typename C>
class CrmsdEvaluator : public Evaluator {
public:

  CrmsdEvaluator(core::data::structural::Structure_SP reference,
      const systems::CartesianAtomsSimple<C> & system) :
      n(reference->count_atoms()), xyz(system) {

    ref = std::make_shared<core::data::basic::Coordinates>(n);
    core::data::structural::structure_to_coordinates(reference,*ref);
    structure_to_coordinates(reference, *ref);
  }

  virtual core::real evaluate() { return rms.crmsd(xyz.coordinates, *ref, n); }

  virtual const std::string & name() const { return name_; }

  virtual core::index1 precision() const { return 2; }

  virtual ~CrmsdEvaluator() {}

  virtual core::index2 min_width() const { return 8; }

private:
  const core::index4 n;
  static const std::string name_;
  core::data::basic::Coordinates_SP ref;
  const systems::CartesianAtomsSimple<C> & xyz;
  mutable core::calc::structural::transformations::Crmsd<std::unique_ptr<core::data::basic::Vec3[]>,core::data::basic::Coordinates> rms;
};

template<typename C>
const std::string CrmsdEvaluator<C>::name_ = "CrmsdEvaluator";

} // ~ cartesian
} // ~ evaluators
} // ~ simulations
#endif
