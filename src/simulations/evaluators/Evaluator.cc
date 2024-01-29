#include "Evaluator.hh"
#include <algorithm>

namespace simulations {
namespace evaluators {

Evaluator::~Evaluator() {
    // Destructor implementation (if needed)
}

std::string Evaluator::header() const {
    core::index2 l = width() - name().length();
    std::string headr = (l > 1) ? std::string(size_t(l / 2), ' ') : "";
    headr += name();
    if (l > headr.size()) headr += std::string(l - headr.size(), ' ');
    return headr;
}

core::index2 Evaluator::width() const {
    return std::max(min_width(), core::index2(name().size()));
}

std::ostream &operator<<(std::ostream &out, Evaluator &e) {
    out << std::setw(e.width()) << std::setprecision(e.precision()) << std::fixed << std::showpoint << e.evaluate();
    return out;
}

std::ostream &operator<<(std::ostream &out, Evaluator_SP e) {
    return out << *e;
}

} // namespace evaluators
} // namespace simulations

