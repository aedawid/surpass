#ifndef SIMULATIONS_OBSERVERS_EndVectorObserver_HH
#define SIMULATIONS_OBSERVERS_EndVectorObserver_HH

#include <iostream>
#include <fstream>

#include <simulations/systems/CartesianAtomsSimple.hh>
#include <simulations/observers/ObserverInterface.hh>

namespace simulations {
namespace observers {
namespace cartesian {

/** @brief Observes end-to-end vector.
 */
template<typename C>
class EndVectorObserver : public ObserverInterface {
public:

    /** @brief Creates an observer that writes a conformation to a file in PDB format
   * @param observed_object - system whose coordinates will be observed
   * @param out_fname - name of the output file
   */
    EndVectorObserver(const systems::CartesianAtomsSimple<C> &observed_object, const std::string & out_fname)
      : out_fname(out_fname), observed_object_(observed_object) {
    outstream =   std::make_shared<std::ofstream>(out_fname);
  }

  /// Measure the vector and write to the file
  bool observe();

  void finalize() { (std::static_pointer_cast<std::ofstream>(outstream))->close(); }

  virtual std::shared_ptr<std::ostream> output_stream() { return outstream; };

  virtual void output_stream(std::shared_ptr<std::ostream> out) { outstream = out; };

  virtual core::index4 count_observe_calls() const { return cnt; }

private:
  std::string out_fname;
  core::index4  cnt = 0;
  std::shared_ptr<std::ostream> outstream;
  const simulations::systems::CartesianAtomsSimple<C> &observed_object_;
};

}
}
}
#endif
