/** @file ToStreamObserver.hh
 *  @brief Observer that can change a stream used to write observations.
 */
#ifndef SIMULATIONS_OBSERVERS_ToStreamObserver_HH
#define SIMULATIONS_OBSERVERS_ToStreamObserver_HH

#include <memory>
#include <ostream>

#include <simulations/observers/ObserverInterface.hh>

namespace simulations {
namespace observers {

/** @brief Observer which can expose the pointer to its output stream.
 *
 * It is also possible to substitute stream pointer with another one, which effectively vhanges the destination
 * where observations are written.
 */
class ToStreamObserver : public ObserverInterface {
public:

  /// Returns a pointer to the stream this observer uses to write observations
  virtual std::shared_ptr<std::ostream> output_stream() = 0;

  /// Change a stream to write observations from this observer
  virtual void output_stream(std::shared_ptr<std::ostream> out) = 0;
};

}
}

#endif
