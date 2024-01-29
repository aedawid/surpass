#ifndef SIMULATIONS_GENERIC_EVALUATOR_HH
#define SIMULATIONS_GENERIC_EVALUATOR_HH

#include <memory>
#include <iostream>
#include <iomanip>
#include <string>
#include <core/real.hh>
#include <core/index.hh>

namespace simulations {
namespace evaluators {

class Evaluator {
public:
    virtual core::real evaluate() = 0; // Pure virtual function
    virtual const std::string &name() const = 0; // Pure virtual function
    virtual std::string header() const; // Virtual function with definition
    virtual core::index1 precision() const { return 2; } // Virtual function with default implementation
    virtual core::index2 width() const; // Virtual function with definition
    virtual core::index2 min_width() const { return 5; } // Virtual function with default implementation
    virtual ~Evaluator(); // Virtual destructor with definition in .cc file
};

typedef std::shared_ptr<Evaluator> Evaluator_SP;
std::ostream &operator<<(std::ostream &out, Evaluator &e);

} // namespace evaluators
} // namespace simulations

#endif // SIMULATIONS_GENERIC_EVALUATOR_HH

