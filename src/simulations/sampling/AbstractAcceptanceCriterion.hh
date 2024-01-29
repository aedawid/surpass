#ifndef SIMULATIONS_GENERIC_SAMPLING_AbstractAcceptanceCriterion_HH
#define SIMULATIONS_GENERIC_SAMPLING_AbstractAcceptanceCriterion_HH

#include <core/real.hh>

namespace simulations {
namespace sampling {


class AbstractAcceptanceCriterion {
public:
	virtual bool test(const core::real old_energy,const core::real new_energy) = 0;
	virtual ~AbstractAcceptanceCriterion() {};
};

} // ~ movers
} // ~ simulations

#endif
