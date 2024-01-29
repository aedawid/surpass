#ifndef CORE_CALC_NUMERIC_basic_math_HH
#define CORE_CALC_NUMERIC_basic_math_HH

#include <vector>

#include <core/real.hh>
#include <core/index.hh>

namespace core {
namespace calc {
namespace numeric {

/** @brief signum function.
 */
template<typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }

double gamma_function(const double x);

/** @brief Evaluates logarithm of gamma function
 *
 */
long double log_gamma_function(const double x);

/** @brief Evaluates incomplete gamma function
 *
 */
double incomplete_gamma_function(const double S, const double Z);

/** @brief Evaluates logarithm of incomplete gamma function
 */
long double log_incomplete_gamma_function(const long double S, const long double Z);

/** @brief evaluates chi-square p-value
 * @param n_dof - number of degrees of freedom
 * @param value - chi-square statistics value
 */
double chi_square_pvalue(const int n_dof, const double value);

/** @brief Evaluates the modified Bessel function of the first kind of order 0
 *
 * @param x the value to compute the Bessel function
 */
double mod_bessel_first_kind_zero(double x);

/** @brief Evaluates the modified Bessel function of the first kind of order 1
 *
 * @param x the value to compute the Bessel function
 */
double mod_bessel_first_kind_one(double x);

/** @brief Evaluates the Erf function
 *
 * @param x the value to compute the erf function
 */
double erf(double x);

/** @brief Calculates mean of angular data
 *
 * @param data - input data
 * @return angular average computed as \f$ \bar{\alpha} = \mathrm{atan2}\left(\frac{1}{n}\cdot\sum_{j=1}^n \sin\alpha_j,\frac{1}{n}\cdot\sum_{j=1}^n \cos\alpha_j\right) \f$
 */
double circular_mean(const std::vector<core::real> & data);

/** Fill the given vector with <code>n</code> equally spaced double values
 *
 * @param min - smallest value; it will be placed as the very first one in the vector
 * @param max - largest value; it will be placed as the very last one in the vector
 * @param n - the total number of values, including the two given ones
 * @param result - where the generated numbers will be stored; will  be cleared before the output values are inserted
 */
std::vector<core::real> & evenly_spaced_values(core::real min, core::real max, core::index4 n, std::vector<core::real> & result);

/** @brief Finds roots of a cubic polynomial
 *
 * @param coeff - four coefficients of the polynomial in the following order: \f$ a_0, a_1, a_2, a_3 \f$
 * @param roots - real roots will be stored in this vector
 * @tparam T - either float or double
 * @return the number of real roots found: 1 or 3
 */
template <typename T>
void find_cubic_roots(const std::vector<T> & coeff, std::vector<T> & roots);

/** @brief Finds roots of a cubic polynomial
 *
 * @param coeff - four coefficients of the polynomial in the following order: \f$ a_0, a_1, a_2\f$
 * i.e. \f$ c, b, a\f
 * @param roots - real roots will be stored in this vector
 * @tparam T - either float or double
 * @return the number of real roots found: 1 or 3
 */
template <typename T>
void find_quadratic_roots(const std::vector<T> & coeff, std::vector<T> & roots);

}
}
}

#endif
