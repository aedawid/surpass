#include <vector>
#include <iostream>
#include <iomanip>

#include <utils/Logger.hh>

#include <simulations/evaluators/Evaluator.hh>
#include <simulations/observers/ObserveEvaluators.hh>

namespace simulations {
namespace observers {


void ObserveEvaluators::add_evaluator(evaluators::Evaluator_SP evaluator) {

  evaluators.push_back(evaluator);
  sw.push_back(evaluator->width());
  logger << utils::LogLevel::INFO << "observing new evaluator: " << evaluator->name() << "\n";
}

std::string ObserveEvaluators::header_string() const {

  std::stringstream ss;
  for (const auto &e : evaluators) ss << " " << e->header() << " ";

  return ss.str();
}

void ObserveEvaluators::finalize() {

  if(is_file_) {
    std::shared_ptr<std::ofstream> of = std::static_pointer_cast<std::ofstream>(outstream);
    of->close();
  } else outstream->flush();
}

bool ObserveEvaluators::observe() {

  ++cnt;
  if(!trigger->operator()()) return false;

  (*outstream) << std::fixed << std::setw(5) << cnt;
  for (size_t i = 0; i < evaluators.size(); ++i)
    (*outstream) << ' ' << std::setw(sw[i]) << std::setprecision(evaluators[i]->precision()) << evaluators[i]->evaluate();
  (*outstream) << "\n";
  outstream->flush();

  return true;
}

} // ~ simulations
} // ~ observers
