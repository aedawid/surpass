#ifndef CORE_CALC_STATISTICS_PDF_HH
#define CORE_CALC_STATISTICS_PDF_HH

#include <ostream>
#include <vector>
#include <core/real.hh>
#include <core/index.hh>

namespace core {
namespace calc {
namespace statistics {

/** \brief A virtual base class for a probability distribution function (pdf).
 */
class PDF {
public:

  /** @brief Constructor creates an instance of this distribution based on the given parameters
   * @param parameters - parameters of the distribution, e.g. expectation, variance, moments, etc.
   */
  PDF(const std::vector<core::real> & parameters) : parameters_(parameters) {}

   /// Default virtual destructor does nothing
  virtual ~PDF() {}

  /** @brief Estimates parameters of this distribution based on a given sample.
   *
   * @param observations - a sample of random observations used for the estimation
   */
  virtual std::vector<core::real> & estimate(const std::vector<std::vector<core::real>> & observations) = 0;

  /** @brief Evaluate the probability of withdrawing a given value from distribution.
   *
   * @param random_value - PDF function argument
   */
  virtual core::real evaluate(const std::vector<core::real> & random_value) const = 0;

  /** @brief Evaluate the probability of withdrawing a given value from distribution.
   *
   * This operator simply calls <code>evaluate()</code>
   * @param random_value - PDF function argument
   */
  core::real operator()(const std::vector<core::real> & random_value) const { return evaluate(random_value); }

  /** @brief Provide read-only access to parameters of this distribution
   *
   * @return vector of distribution's parameters, such as mean, standard deviation, etc.
   */
  const std::vector<core::real> & parameters() const { return parameters_; }

  /** @brief Writes parameters of a distribution into a stream
   *
   * @param stream - output stream
   * @param pdf - distribution to be written
   * @return reference to the output stream
   */
  friend std::ostream &operator<<(std::ostream &stream, const PDF &pdf) {

    for (real p : pdf.parameters_) stream << p << " ";

    return stream;
  }

  /** @brief Copies parameters of this distribution to the given <code>destination</code> vector.
   * <code>destination</code> vector must be at least of the same size as the number of parameters of this distribution
   * @param destination - vector where parameters are copied to
   */
  void get_parameters(std::vector<core::real> &destination) const {

    for (index1 i = 0; i < parameters_.size(); ++i) destination[i] = parameters_[i];
  }

  /** @brief Copies parameters of this distribution to the given <code>destination</code> vector.
   *
   * <code>destination</code> vector must be at least of the same size as the number of paarmeters of this distribution.
   * This function is the opposite to <code>get_parameters(std::vector<core::real> &)</code>
   *
   * @param destination - vector where parameters are copied to
   */
  virtual void set_parameters(const std::vector<core::real> &source) = 0;

protected:
  ///  parameters of the distribution, e.g. expectation, variance, moments, etc.
  std::vector<core::real> parameters_;
};

}
}
}

#endif
