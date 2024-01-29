#ifndef SIMULATIONS_CARTESIAN_EVALUATORS_CM_HH
#define SIMULATIONS_CARTESIAN_EVALUATORS_CM_HH

#include <memory>

#include <core/real.hh>
#include <core/data/basic/Vec3.hh>
#include <utils/string_utils.hh>
#include <simulations/evaluators/Evaluator.hh>
#include <simulations/systems/CartesianAtomsSimple.hh>

namespace simulations {
namespace evaluators {
namespace cartesian {

/** \brief Evaluator that reports the distance made by the center of mass of a monitored system.
 *
 * This evaluator also reports the vector of CM location
 */
template <typename C>
class CM : public simulations::evaluators::Evaluator {
public:

  /** \brief Create the evaluator for a given system.
   *
   * @param system - system to be monitored
   */
  CM(const systems::CartesianAtomsSimple<C> & system) : xyz(system) {}

  /** \brief Evaluate the monitored property i.e. the distance of the CM from the origin
   *
   * @return the distance of the CM from the origin
   */
  virtual core::real evaluate() {

    cm.set(0);
    for (size_t i=0;i<xyz.count_atoms();++i)
      cm += xyz[i];

    cm /= double(xyz.count_atoms());

    return cm.length();
  }

  /** \brief Evaluate the position of the center of mass (CM)
   *
   * @param result -  the position of the CM will be stored here
   */
  void evaluate_vector(core::data::basic::Vec3 & result) {

    result.set(0.0);
    for (size_t i=0;i<xyz.count_atoms();++i)
      result += xyz[i];
    result /= double(xyz.count_atoms());
  }

  virtual std::string header() const {

    core::index2 l = width() + 2;
    std::string header(size_t(l / 2), ' ');
    header += name_;
    header += std::string(l - header.size(), ' ');

    return header;
 }

  /// Returns the name of this evaluator, which is <code>CM</code>
  virtual const std::string & name() const { return name_; }

  virtual core::index1 precision() const { return 2; }

  virtual core::index2 min_width() const { return 7; }

  virtual ~CM() {}

  /// Returns the X coordinate of the center of mass for a monitored system that has been measured at the very recent <code>evaluate()</code> call
  core::real recent_x() const { return cm.x; }
  /// Returns the Y coordinate of the center of mass for a monitored system that has been measured at the very recent <code>evaluate()</code> call
  core::real recent_y() const { return cm.y; }
  /// Returns the Z coordinate of the center of mass for a monitored system that has been measured at the very recent <code>evaluate()</code> call
  core::real recent_z() const { return cm.z; }


private:
  mutable C cm;
  const systems::CartesianAtomsSimple<C> & xyz;
  static const std::string name_;
};

template <typename C>
const std::string CM<C>::name_ = "CM";

template <typename C>
std::ostream & operator<<(std::ostream &out, CM<C> & e) {

  e.evaluate();
  out << std::setw(e.width()) << std::setprecision(e.precision())<<e.recent_x();
  out << " " << std::setw(e.width()) << std::setprecision(e.precision()) << e.recent_y();
  out << " " << std::setw(e.width()) << std::setprecision(e.precision()) << e.recent_z();
  return out;
}


} // ~ cartesian
} // ~ evaluators
} // ~ simulations
#endif
