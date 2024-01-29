#include <iostream>
#include <utils/Logger.hh>

#include <simulations/systems/AtomTypingVariants.hh>

namespace simulations {
namespace systems {

std::ostream & operator<<(std::ostream & out,const AtomTypingVariants v) {
  switch(v) {
    case AtomTypingVariants::C_TERMINAL : out << "C_TERMINAL"; return out;
    case AtomTypingVariants::N_TERMINAL : out << "N_TERMINAL"; return out;
    case AtomTypingVariants::STANDARD : out << "STANDARD"; return out;
    default: out << "UNKNOWN"; return out;
  }
}

utils::Logger & operator<<(utils::Logger & logs,const AtomTypingVariants v) {

  switch(v) {
    case AtomTypingVariants::C_TERMINAL : logs << "C_TERMINAL"; return logs;
    case AtomTypingVariants::N_TERMINAL : logs << "N_TERMINAL"; return logs;
    case AtomTypingVariants::STANDARD : logs << "STANDARD"; return logs;
    default: logs << "UNKNOWN"; return logs;
  }
}

}
}
