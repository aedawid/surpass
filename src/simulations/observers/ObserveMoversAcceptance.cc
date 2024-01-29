#include <iostream>
#include <iomanip>

#include <simulations/observers/ObserveMoversAcceptance.hh>

namespace simulations {
namespace observers {

bool ObserveMoversAcceptance::observe() {

  ++cnt;
  if(!ObserverInterface::trigger->operator()()) return false;

  *(outstream) << std::setw(5) << cnt << " "<< ms_ << "\n";
  outstream->flush();
  return true;
}

void ObserveMoversAcceptance::finalize() {

  if (is_file_) {
    std::shared_ptr<std::ofstream> of = std::static_pointer_cast<std::ofstream>(outstream);
    of->close();
  } else outstream->flush();
}

} // ~ simulations
} // ~ observers
