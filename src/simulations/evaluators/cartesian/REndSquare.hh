#ifndef SIMULATIONS_CARTESIAN_EVALUATORS_REndSquare_HH
#define SIMULATIONS_CARTESIAN_EVALUATORS_REndSquare_HH

#include <memory>

#include <core/real.hh>
#include <core/index.hh>
#include <core/data/basic/Vec3.hh>
#include <utils/string_utils.hh>
#include <simulations/evaluators/Evaluator.hh>
#include <simulations/systems/CartesianAtomsSimple.hh>

namespace simulations {
namespace evaluators {
namespace cartesian {

/** \brief Evaluator that reports the end-to-end length.
 *
 * This evaluator also reports the vector of end-to-end vector.
 */
template <typename C>
class REndSquare : public simulations::evaluators::Evaluator {
public:

  /** \brief Create the evaluator for a given system.
   *
   * @param system - system to be monitored
   */
  REndSquare(const systems::CartesianAtomsSimple<C> & system) : xyz(system) { n = xyz.count_atoms(); }

  /** \brief Evaluate the length of an end-to-end vector for a single chain
   *
   * @return the distance of the CM from the origin
   */
  virtual core::real evaluate() {

    r_end.set(xyz[n - 1]);
    r_end -= xyz[0];
    return r_end.length();
  }

  /// Returns the X coordinate of the center of mass for a monitored system that has been measured at the very recent <code>evaluate()</code> call
  core::real recent_x() const { return r_end.x; }
  /// Returns the Y coordinate of the center of mass for a monitored system that has been measured at the very recent <code>evaluate()</code> call
  core::real recent_y() const { return r_end.y; }
  /// Returns the Z coordinate of the center of mass for a monitored system that has been measured at the very recent <code>evaluate()</code> call
  core::real recent_z() const { return r_end.z; }

  /** \brief Evaluate the vector between the first and last atom of the monitored system
   *
   * @param result -  the position of the CM will be stored here
   */
  void evaluate_vector(core::data::basic::Vec3 & result) { result.set(xyz[n-1]); result -= xyz[0];}

  /// Returns a header for a column for results of this evaluator
  virtual std::string header() const {

    core::index2 l = width() * 3 + 2;
    std::string header(size_t(l / 2), ' ');
    header += name_;
    header += std::string(l - header.size(), ' ');

    return header;
 }

  /// Returns the name of this evaluator, which is <code>REndSquare</code>
  virtual const std::string & name() const { return name_; }

  virtual core::index1 precision() const { return 2; }

  virtual core::index2 min_width() const { return 7; }

  /// Virtual destructor to satisfy a compiler
  virtual ~REndSquare() {}

private:
  mutable C r_end;
  core::index4 n;
  const systems::CartesianAtomsSimple<C> & xyz;
  static const std::string name_;
};

template <typename C>
const std::string REndSquare<C>::name_ = "REndSquare";

template <typename C>
std::ostream & operator<<(std::ostream &out, REndSquare<C> & e) {

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
