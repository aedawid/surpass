#include <math.h>
#include <algorithm>    // std::sort
#include <vector>
#include <iostream>
#include <core/index.hh>
#include <core/real.hh>
#include <core/calc/statistics/NormalDistribution.hh>

namespace core {
namespace calc {
namespace statistics {

static bool sort_by_x(const std::vector<core::real>& p1, const std::vector<core::real>& p2) {return p1[0] < p2[0];}             // This makes the sort be according to column 1 (x) and ascending

NormalDistribution::NormalDistribution(const std::vector<core::real> &parameters) : PDF(parameters) {

  if (parameters_.size() == 0) {
    parameters_.push_back(0.0);
    parameters_.push_back(1.0);
  }
  set_up_constants();
}

void NormalDistribution::set_parameters(const std::vector<core::real> &source) {

  for(unsigned int i = 0; i < parameters_.size(); ++i) {
    parameters_[i] = source[i];
  }
  set_up_constants();
}

void NormalDistribution::set_up_constants() { const_1 = (1.0 / (parameters_[1] * (sqrt(M_PI * 2.0)))); }

std::vector<core::real> & NormalDistribution::estimate(const std::vector<std::vector<core::real>> &observations) {
  return estimate(observations, 0);
}

std::vector<core::real> & NormalDistribution::estimate(const std::vector<std::vector<core::real>> &observations, const core::index1 col) {

//    calculate average and sdev
  real sum = 0;
  real sum2 = 0;
  real n_data = observations.size();
  for (core::index4 i = 0; i < observations.size(); ++i) {
    sum += observations[i][col];
    sum2 += observations[i][col] * observations[i][col];
  }

  sum /= n_data;
  sum2 = sum2 / n_data - sum * sum;
  parameters_[0] = sum;                       // average
  parameters_[1] = sqrt(sum2);                // standard deviation x

  set_up_constants();

  return parameters_;
}

std::vector<core::real> & NormalDistribution::robust_estimate(std::vector<std::vector<core::real>> &observations) {
  return robust_estimate(observations, 0);
}

std::vector<core::real> & NormalDistribution::robust_estimate(std::vector<std::vector<core::real>> &observations, const core::index1 col) {

//    calculate average and sdev
  real sum = 0;
  real sum2 = 0;

//    sort_by_x & erase observations before calculate the average trimming
  bool is_begin = true, is_end = true;
  std::sort(observations.begin(), observations.end(), sort_by_x);
  for (unsigned k = 0; k <= 100; ++k) {
    if ((is_begin == false) && (is_end == false)) break;
    if ((observations.size() <= 100) && (k == observations.size() / 4)) break;
    if ((is_begin == true) && (observations[0][col] < parameters_[0] - (2 * parameters_[1])))
      observations.erase(observations.begin());
    if ((is_end == true) &&
        (observations[observations.size() - 1][col] > parameters_[0] + (2 * parameters_[1])))
      observations.erase(observations.end() - 1);
  }
  for (core::index4 i = 0; i < observations.size(); ++i) {
    sum += observations[i][0];
    sum2 += observations[i][0] * observations[i][col];
  }
  real n_data = observations.size();
  sum /= n_data;
  sum2 = sum2 / n_data - sum * sum;
  parameters_[0] = sum;                       // average
  parameters_[1] = sqrt(sum2);                // standard deviation

  set_up_constants();

  return parameters_;
}

/** @brief Evaluate the probability of withdrawing a given value from distribution.
 *
 * @param random_value - PDF function argument
 */
core::real NormalDistribution::evaluate(core::real x) const {

  x -= parameters_[0];

  core::real pdf = const_1 * exp(-(x * x) / (2 * parameters_[1] * parameters_[1]));

  return pdf;
}

}
}
}

