#include <iostream>
#include <math.h>

#include <core/calc/numeric/basic_math.hh>

namespace core {
namespace calc {
namespace numeric {


double gamma_function(const double N) {

  const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;

  long double Z = (long double) N;
  long double Sc = powl((Z + 15), (Z + 0.5));
  Sc *= expl(-1.0 * (Z + 15));
  Sc /= Z;

  long double F = 1.0;
  long double Ck;
  long double Sum = SQRT2PI;


  for (int K = 1; K < 15; K++) {
    Z++;
    Ck = powl(15 - K, K - 0.5);
    Ck *= expl(15 - K);
    Ck /= F;

    Sum += (Ck / Z);

    F *= (-1.0 * K);
  }

  return (double) (Sum * Sc);
}

long double log_gamma_function(const double N) {

  const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;
  long double Z = (long double) N;
  long double Sc;

  Sc = (logl(Z + 15) * (Z + 0.5)) - (Z + 15) - logl(Z);

  long double F = 1.0;
  long double Ck;
  long double Sum = SQRT2PI;

  for (int K = 1; K < 15; K++) {
    Z++;
    Ck = powl(15 - K, K - 0.5);
    Ck *= expl(15 - K);
    Ck /= F;

    Sum += (Ck / Z);

    F *= (-1.0 * K);
  }

  return logl(Sum) + Sc;
}

static long double KM(const long double S,const long double Z) {
  long double Sum = 1.0;
  long double Nom = 1.0;
  long double Denom = 1.0;
  long double s_ = S;
  for (int I = 0; I < 1000; I++) { // Loops for 1000 iterations
    Nom *= Z;
    s_++;
    Denom *= s_;
    Sum += (Nom / Denom);
  }

  return Sum;
}

long double log_incomplete_gamma_function(const long double S,const long double Z) {

  if (Z < 0.0) return 0.0;
  long double Sc, K;
  Sc = (logl(Z) * S) - Z - logl(S);

  K = KM(S, Z);

  return logl(K) + Sc;
}

double incomplete_gamma_function(const double S, const double Z) {

  if (Z < 0.0) return 0.0;
  long double s_ = S;
  double Sc = (1.0 / s_);
  Sc *= pow(Z, s_);
  Sc *= exp(-Z);

  double Sum = 1.0;
  double Nom = 1.0;
  double Denom = 1.0;

  for (int I = 0; I < 200; I++) {
    Nom *= Z;
    s_++;
    Denom *= s_;
    Sum += (Nom / Denom);
  }

  return Sum * Sc;
}

double erf(double x) {
  double y = 1.0 / (1.0 + 0.3275911 * x);
  return 1 - (((((+1.061405429 * y - 1.453152027) * y + 1.421413741) * y - 0.284496736) * y + 0.254829592) * y) *
             exp(-x * x);
}

double chi_square_pvalue(const int Dof, const double Cv) {

  if (Cv < 0 || Dof < 1) return 0.0;
  if (Dof == 2) return exp(-1.0 * Cv * 0.5);

  long double PValue = 1.0L - expl(log_incomplete_gamma_function(((double) Dof) * 0.5, Cv * 0.5)) / gamma_function(((double) Dof) * 0.5);

  return (double) PValue;
}

/**
 * Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */
static constexpr double A_i0[] = { -4.41534164647933937950E-18, 3.33079451882223809783E-17,
                                       -2.43127984654795469359E-16, 1.71539128555513303061E-15, -1.16853328779934516808E-14, 7.67618549860493561688E-14,
                                       -4.85644678311192946090E-13, 2.95505266312963983461E-12, -1.72682629144155570723E-11, 9.67580903537323691224E-11,
                                       -5.18979560163526290666E-10, 2.65982372468238665035E-9, -1.30002500998624804212E-8, 6.04699502254191894932E-8,
                                       -2.67079385394061173391E-7, 1.11738753912010371815E-6, -4.41673835845875056359E-6, 1.64484480707288970893E-5,
                                       -5.75419501008210370398E-5, 1.88502885095841655729E-4, -5.76375574538582365885E-4, 1.63947561694133579842E-3,
                                       -4.32430999505057594430E-3, 1.05464603945949983183E-2, -2.37374148058994688156E-2, 4.93052842396707084878E-2,
                                       -9.49010970480476444210E-2, 1.71620901522208775349E-1, -3.04682672343198398683E-1, 6.76795274409476084995E-1 };

/**
 * Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */
static constexpr double B_i0[] = {
  -7.23318048787475395456E-18, -4.83050448594418207126E-18,
                                       4.46562142029675999901E-17, 3.46122286769746109310E-17, -2.82762398051658348494E-16, -3.42548561967721913462E-16,
                                       1.77256013305652638360E-15, 3.81168066935262242075E-15, -9.55484669882830764870E-15, -4.15056934728722208663E-14,
                                       1.54008621752140982691E-14, 3.85277838274214270114E-13, 7.18012445138366623367E-13, -1.79417853150680611778E-12,
                                       -1.32158118404477131188E-11, -3.14991652796324136454E-11, 1.18891471078464383424E-11, 4.94060238822496958910E-10,
                                       3.39623202570838634515E-9, 2.26666899049817806459E-8, 2.04891858946906374183E-7, 2.89137052083475648297E-6,
                                       6.88975834691682398426E-5, 3.36911647825569408990E-3, 8.04490411014108831608E-1 };

static constexpr double B_i1[] = {
  7.51729631084210481353E-18,
  4.41434832307170791151E-18,
  -4.65030536848935832153E-17,
  -3.20952592199342395980E-17,
  2.96262899764595013876E-16,
  3.30820231092092828324E-16,
  -1.88035477551078244854E-15,
  -3.81440307243700780478E-15,
  1.04202769841288027642E-14,
  4.27244001671195135429E-14,
  -2.10154184277266431302E-14,
  -4.08355111109219731823E-13,
  -7.19855177624590851209E-13,
  2.03562854414708950722E-12,
  1.41258074366137813316E-11,
  3.25260358301548823856E-11,
  -1.89749581235054123450E-11,
  -5.58974346219658380687E-10,
  -3.83538038596423702205E-9,
  -2.63146884688951950684E-8,
  -2.51223623787020892529E-7,
  -3.88256480887769039346E-6,
  -1.10588938762623716291E-4,
  -9.76109749136146840777E-3,
  7.78576235018280120474E-1
};

static constexpr double A_i1[] {
  2.77791411276104639959E-18,
  -2.11142121435816608115E-17,
  1.55363195773620046921E-16,
  -1.10559694773538630805E-15,
  7.60068429473540693410E-15,
  -5.04218550472791168711E-14,
  3.22379336594557470981E-13,
  -1.98397439776494371520E-12,
  1.17361862988909016308E-11,
  -6.66348972350202774223E-11,
  3.62559028155211703701E-10,
  -1.88724975172282928790E-9,
  9.38153738649577178388E-9,
  -4.44505912879632808065E-8,
  2.00329475355213526229E-7,
  -8.56872026469545474066E-7,
  3.47025130813767847674E-6,
  -1.32731636560394358279E-5,
  4.78156510755005422638E-5,
  -1.61760815825896745588E-4,
  5.12285956168575772895E-4,
  -1.51357245063125314899E-3,
  4.15642294431288815669E-3,
  -1.05640848946261981558E-2,
  2.47264490306265168283E-2,
  -5.29459812080949914269E-2,
  1.02643658689847095384E-1,
  -1.76416518357834055153E-1,
  2.52587186443633654823E-1
};

static double chbevl(double x, const double* coef,const int N) {

  double b0, b1, b2;
  int p = 0;
  int i;
  b0 = coef[p++];
  b1 = 0.0;
  i = N - 1;
  do {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + coef[p++];
  } while (--i > 0);

  return (0.5 * (b0 - b2));
}

double mod_bessel_first_kind_zero(double x) {

  double y;
  if (x < 0)
    x = -x;
  if (x <= 8.0) {
    y = (x / 2.0) - 2.0;
    return (exp(x) * chbevl(y, A_i0, 30));
  }
  return (exp(x) * chbevl(32.0 / x - 2.0, B_i0, 25) / sqrt(x));
}

double mod_bessel_first_kind_one(double x) {
  double y, z;

  z = abs(x);
  if (z <= 8.0) {
    y = (z / 2.0) - 2.0;
    z = chbevl(y, A_i1, 29) * z * exp(z);
  } else {
    z = exp(z) * chbevl(32.0 / z - 2.0, B_i1, 25) / sqrt(z);
  }
  if (x < 0.0)
    z = -z;
  return (z);
}

double circular_mean(const std::vector<core::real> & data) {

  double s = 0, c = 0;
  for (real a:data) {
    s += sin(a);
    c += cos(a);
  }
  s /= double(data.size());
  c /= double(data.size());

  return atan2(s,c);
}

std::vector<core::real> & evenly_spaced_values(core::real min, core::real max, core::index4 n, std::vector<core::real> & result) {

  result.clear();
  double delta = (max-min)/double(n-1);
  for(core::index4 i=0; i<n; i++) result.push_back(min + i*delta);

  return result;
}

/* ax^2 + bx + c = 0
   x = ( -b +- sqrt(b^2 - 4ac) ) / 2a */
template<typename T>
void find_quadratic_roots(const std::vector<T> & coeff, std::vector<T> & roots) {

  if (roots.size() > 0) roots.clear();

  if (coeff[0] == 0.0) {
    roots.push_back(coeff[2] / coeff[1]);
  } else {
    const double delta = (coeff[1] * coeff[1]) - (4 * coeff[0] * coeff[2]);

    if (delta == 0.0) {
      roots.push_back(-coeff[1] / (2 * coeff[0]));
    } else {
      const double sqrt_delta = sqrt(delta);
      const double two_a = 2 * coeff[0];
      roots.push_back((-coeff[1] + sqrt_delta) / two_a);
      roots.push_back((-coeff[1] - sqrt_delta) / two_a);
    }
  }
}

/* Copyright (C) 1997-2001 Ken Turkowski. <turk@computer.org>
 *
 * All rights reserved.
 *
 * Warranty Information
 *  Even though I have reviewed this software, I make no warranty
 *  or representation, either express or implied, with respect to this
 *  software, its quality, accuracy, merchantability, or fitness for a
 *  particular purpose.  As a result, this software is provided "as is,"
 *  and you, its user, are assuming the entire risk as to its quality
 *  and accuracy.
 *
 * This code may be used and freely distributed as long as it includes
 * this copyright notice and the above warranty information.
 */

/*******************************************************************************
 * FindCubicRoots
 *
 *  Solve:
 *      coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 *
 *  returns:
 *      3 - 3 real roots
 *      1 - 1 real root (2 complex conjugate)
 *******************************************************************************/
template<typename T>
void find_cubic_roots(const std::vector<T> & coeff, std::vector<T> & roots) {

  if (roots.size() > 0) roots.clear();

  if (coeff[3] == 0.0) { // --- Actually it's a simple quadratic equation
    std::vector<T> coeff2{coeff[1], coeff[2], coeff[3]};
    find_quadratic_roots(coeff2, roots);
    return;
  }

  if (coeff[0] == 0.0) { // --- One of the roots is zero, solve a quadratic equation to get the other two roots
    find_quadratic_roots(coeff, roots);
    roots.push_back(0.0);
    return;
  }

  double a1 = coeff[2] / coeff[3];
  double a2 = coeff[1] / coeff[3];
  double a3 = coeff[0] / coeff[3];

  double Q = (a1 * a1 - 3 * a2) / 9;
  double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
  double Qcubed = Q * Q * Q;
  double d = Qcubed - R * R;

  if (d >= 0) { /* Three real roots */
    double theta = acos(R / sqrt(Qcubed));
    double sqrtQ = sqrt(Q);
    roots.push_back( -2 * sqrtQ * cos( theta             / 3) - a1 / 3 );
    roots.push_back( -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3 );
    roots.push_back( -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3 );
  }  else { /* One real root */
    double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
    if (R > 0) e = -e;
    roots.push_back(  (e + Q / e) - a1 / 3. );
  }
}

template
void find_quadratic_roots(const std::vector<double> & coeff, std::vector<double> & roots);
template
void find_quadratic_roots(const std::vector<float> & coeff, std::vector<float> & roots);

template
void find_cubic_roots(const std::vector<double> & coeff, std::vector<double> & roots);
template
void find_cubic_roots(const std::vector<float> & coeff, std::vector<float> & roots);

}
}
}

