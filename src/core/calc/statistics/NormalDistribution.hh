#ifndef CORE_CALC_STATISTICS_NormalDistribution_HH
#define CORE_CALC_STATISTICS_NormalDistribution_HH

#include <math.h>
#include <algorithm>    // std::sort
#include <vector>
#include <core/real.hh>
#include <core/calc/statistics/PDF.hh>

namespace core {
namespace calc {
namespace statistics {

/** \brief Normal distribution function.
 *
 * \include ex_NormalDistribution.cc
 */
class NormalDistribution: public PDF {
public:

  /** @brief Constructor creates an instance of this distribution based on the given parameters : mean and sdev
   * @param parameters - parameters of the distribution, e.g. expectation, variance, moments, etc.
   */
  NormalDistribution(const real avg, const real sdev) : NormalDistribution(std::vector<core::real>({avg,sdev})) {}

  /** @brief Constructor creates an instance of this distribution based on the given parameters
   * @param parameters - parameters of the distribution, e.g. expectation, variance, moments, etc.
   */
  NormalDistribution(const std::vector<core::real> & parameters);

   /// Default virtual destructor does nothing
  ~NormalDistribution() {}

  /** @brief Estimates parameters of this distribution based on a given sample.
   *
   * @param observations - PDF function argument
   */
  virtual std::vector<core::real> & estimate(const std::vector<std::vector<core::real>> & observations);

  /** @brief Estimates parameters of this distribution based on a given sample.
   *
   * @param observations - PDF function argument
   * @param col - a colum index with data
   */
  virtual std::vector<core::real> & estimate(const std::vector<std::vector<core::real>> & observations, const core::index1 col);

  /** @brief Estimates parameters of this distribution based on a given sample using the average trimming as robust estimation algorithm
   *
   * @param observations - PDF function argument
   */
  std::vector<core::real> & robust_estimate(std::vector<std::vector<core::real>> & observations);

 /** @brief Estimates parameters of this distribution based on a given sample using the average trimming as robust estimation algorithm
   *
   * @param observations - PDF function argument
   */
  std::vector<core::real> & robust_estimate(std::vector<std::vector<core::real>> & observations, const core::index1 col);


  /** @brief Evaluate the probability of withdrawing a given vector of values from this distribution.
   *
   * @param random_value - vector of PDF function arguments : X for 1D normal
   */
  virtual inline core::real evaluate(const std::vector<core::real> & random_value) const  {
    return evaluate(random_value[0]);
  }

  /** @brief Evaluate the probability of withdrawing a given value from distribution.
   *
   * @param x_1d - x for the 1D normal term
   */
  core::real evaluate(core::real x_1d) const;

  virtual void set_parameters(const std::vector<core::real> &source);

  /** @brief Returns the probability of \f$ [-\inf,x] \f$ of a normal distribution
   *
   * @param x - function argument (e.g. statistics value to be tested)
   * @param avg - expected value (distribution mean)
   * @param sdev - standard deviation (first central moment of the distribution)
   * @return probability of withdrawing  a value <strong>smaller</strong> than <code>x</code> from a given normal distribution
   */
  static double cdf(double x, double avg, double sdev) { return 0.5 * (1 + erf((x - avg) / (sdev * sqrt(2.)))); }

private:
  double const_1; ///< normalization of the 1D normal distribution

  void set_up_constants();
};

/**
 * \example ex_NormalDistribution.cc
 */

}
}
}

#endif
