#ifndef CORE_CALC_STRUCTURAL_TRANSFORMATIONS_Crmsd_H
#define CORE_CALC_STRUCTURAL_TRANSFORMATIONS_Crmsd_H

#include <memory>
#include <exception>

#include <core/real.hh>
#include <utils/string_utils.hh>
#include <utils/exit.hh>
#include <core/data/basic/Vec3.hh>
#include <core/calc/structural/transformations/Rototranslation.hh>

#include <utils/string_utils.hh>

using core::data::basic::Vec3;

namespace core {
namespace calc {
namespace structural {
namespace transformations {

/** @brief Calculates crsmd between structures.
 * @tparam T1 - type of the vector of coordinates; something like <code>std::unique_ptr<Vec3[]></code> or <code>std::vector<Vec3></code>.
 *    In general - something which allows index coordinates in atoms as: <code>t[i].x</code>
 * @tparam T2 - yet another type providing coordinates; in general might be different that T1
 *
 * @include ex_Crmsd.cc
 */
template<typename T1, typename T2>
class Crmsd : public Rototranslation {
public:

  /** @brief Calculates crmsd value.
   *
   * The calculations will be performed based on the rotation-translation transformation
   * computed during the most recent call of <code>crmsd()</code> method.
   * Neiter <code>templt</code> nor <code>query</code> atoms will be moved
   * @param query - the set of query atoms which remain in their positions
   * @param templt - the set of template atoms which remain in their positions
   * @param n_atoms - the number of atoms in each of the two sets being compared
   */
  double calculate_crmsd_value(const T1 & query, const T2 & templt, const size_t n_atoms) const {

    double ret_val = 0.0;
    for (size_t i = 0; i < n_atoms; i++) {
      core::real x = query[i].x - tr_before().x;
      core::real y = query[i].y - tr_before().y;
      core::real z = query[i].z - tr_before().z;
      const core::real tmpx = x * rot_x().x + y * rot_x().y + z * rot_x().z;
      const core::real tmpy = x * rot_y().x + y * rot_y().y + z * rot_y().z;
      const core::real tmpz = x * rot_z().x + y * rot_z().y + z * rot_z().z;
      x = tmpx + tr_after().x - templt[i].x;
      y = tmpy + tr_after().y - templt[i].y;
      z = tmpz + tr_after().z - templt[i].z;
      ret_val += x * x + y * y + z * z;
    }
    return sqrt(ret_val / (core::real) n_atoms);
  }

  /** @brief Calculates crmsd value <strong>rotate the query set</strong>
   *
   * The set of <code>query</code> atoms will be superimposed on the template set
   * @param query - the set of query atoms which <strong>will be moved</strong>
   * @param templt - the set of template atoms which remain in their positions
   * @param n_atoms - the number of atoms in each of the two sets being compared
   */
  double calculate_crmsd_value_rotate(T1 & query, const T2 & templt, const size_t n_atoms) const {

    double ret_val = 0.0;
    for (size_t i = 0; i < n_atoms; i++) {
      core::real x = query[i].x - tr_before().x;
      core::real y = query[i].y - tr_before().y;
      core::real z = query[i].z - tr_before().z;
      query[i].x = x * rot_x().x + y * rot_x().y + z * rot_x().z + tr_after().x;
      query[i].y = x * rot_y().x + y * rot_y().y + z * rot_y().z + tr_after().y;
      query[i].z = x * rot_z().x + y * rot_z().y + z * rot_z().z + tr_after().z;
      x = query[i].x - templt[i].x;
      y = query[i].y - templt[i].y;
      z = query[i].z - templt[i].z;
      ret_val += x * x + y * y + z * z;
    }
    return sqrt(ret_val / (core::real) n_atoms);
  }

  /** @brief Finds the optimal transformation (rotation and translation) to superimpose the  <code>query</code> atoms on <code>templt</code> atoms.
   *
   * The set of <code>query</code> atoms will be superimposed on the template set
   * @param query - the set of query atoms which remain in their positions
   * @param templt - the set of template atoms which remain in their positions
   * @param n_atoms - the number of atoms in each of the two sets being compared
   * @param store_transformation - if true, the data for the optimal transformation (rotation matrix) will be computed and stored within this object.
   */
  double crmsd(const T1 & query, const T2 & templt, const size_t n_atoms, const bool store_transformation = false) {

    ++call_counter;

#ifdef DEBUG
    runtime_assert(n_atoms > 3, "superimposition needs at least 4 atoms !");
#endif

    const double R = center_covariance_rg2(query, templt, n_atoms);
    covariance2();

    eigenvalues(); // eigenvalues of the covariance_square
    const double determinant = recent_cov[0] * (recent_cov[4] * recent_cov[8] - recent_cov[5] * recent_cov[7])
        - recent_cov[1] * (recent_cov[3] * recent_cov[8] - recent_cov[5] * recent_cov[6])
        + recent_cov[2] * (recent_cov[3] * recent_cov[7] - recent_cov[4] * recent_cov[6]);
    const double d =
        (determinant > 0) ?
            sqrt(recent_eigenvalues[0]) + sqrt(recent_eigenvalues[1]) + sqrt(fabs(recent_eigenvalues[2])) :
            sqrt(recent_eigenvalues[0]) + sqrt(recent_eigenvalues[1]) - sqrt(fabs(recent_eigenvalues[2]));
    double rms = sqrt(R - d - d);

    rms = (rms < 1e-6) ? 0.0 : rms;
    recent_crmsd_ = rms;

    if (store_transformation) {
      eigenvectors();
      rotationMatrix();
      update_mm();
#ifdef DEBUG
      real rms_val = calculate_crmsd_value(query,templt,n_atoms);
      if(fabs(rms_val-rms)>0.01) {
        std::cerr << utils::string_format("Crmsd from covariance and determinant differs from the value after rotation (%f != %f)!\n",rms,rms_val);
        std::cerr << "eigenvalues: "<<recent_eigenvalues[0]<<" "<<recent_eigenvalues[1]<<" "<<recent_eigenvalues[2]<<"\n";
        std::cerr << "determinant: "<<determinant<<"\n";
        std::cerr << "R, d: "<<R<<" "<<d<<" "<<(R-d-d)<<" "<<rms<<"\n";
        throw std::runtime_error(utils::string_format("Crmsd from covariance and determinant differs from the value after rotation (%f != %f)!\n",rms,rms_val));
      }
#endif
    }
    return rms;
  }

  /** @name Methods for benchmarking code that uses the crmsd calculator
   */
  ///@{
  /** @brief Says how many times the <code>crmsd()</code> method was called since the last reset of the internal counter
   */
  inline index4 crmsd_calls_counter() const {
    return call_counter;
  }

  /** @brief Resets the internal counter for <code>crmsd()</code> calls
   */
  inline void reset_crmsd_calls_counter() { call_counter = 0; }
  ///@}

  /// Returns the crmsd value from the most recent crmsd() method call
  inline double recent_crmsd() { return recent_crmsd_; }

private:
  core::index4 call_counter = 0;
  double recent_crmsd_ = 0.0; ///< crmsd value found in the most recent crsmd() method call
  // Covariance matrix packed row-wise
  double recent_cov[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  // The square of the covariance matrix packed row-wise
  double recent_cov2[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  // Eigenvectors packed row-wise
  double recent_eigenvectors[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  double recent_eigenvalues[3] = { 0, 0, 0 };

  inline double center_covariance_rg2(const T1 & query, const T2 & templt, const size_t n_atoms) {

    double recent_n_atoms_inv = 1.0 / double(n_atoms);
    // ---------- Find centers for both structures
    Vec3 tb;
    Vec3 ta;

    for (size_t i = 0; i < n_atoms; i++)
      tb+= query[i];
    for (size_t i = 0; i < n_atoms; i++)
      ta += templt[i];

    tb *= recent_n_atoms_inv;
    ta *= recent_n_atoms_inv;
    Rototranslation::tr_before(tb);
    Rototranslation::tr_after(ta);

    // ---------- Find the value of S2
    // ---------- Also, compute the covariance matrix
    for (core::index2 i = 0; i < 9; i++)
      recent_cov[i] = 0;
    double sum = 0.0;
    Vec3 t, q;
    for (size_t i = 0; i < n_atoms; i++) {
      t.set(templt[i]);
      t -= Rototranslation::tr_after();
      q.set(query[i]);
      q -= Rototranslation::tr_before();
      sum += t.x * t.x + q.x * q.x;
      sum += t.y * t.y + q.y * q.y;
      sum += t.z * t.z + q.z * q.z;

      recent_cov[0] += t.x * q.x;
      recent_cov[1] += t.x * q.y;
      recent_cov[2] += t.x * q.z;

      recent_cov[3] += t.y * q.x;
      recent_cov[4] += t.y * q.y;
      recent_cov[5] += t.y * q.z;

      recent_cov[6] += t.z * q.x;
      recent_cov[7] += t.z * q.y;
      recent_cov[8] += t.z * q.z;
    }
    for (size_t i = 0; i < 9; i++)
      recent_cov[i] *= recent_n_atoms_inv;
    sum *= recent_n_atoms_inv;

    return sum;
  }

  inline void covariance2() {

    recent_cov2[0] = recent_cov[0] * recent_cov[0] + recent_cov[3] * recent_cov[3] + recent_cov[6] * recent_cov[6];
    recent_cov2[1] = recent_cov[0] * recent_cov[1] + recent_cov[3] * recent_cov[4] + recent_cov[6] * recent_cov[7];
    recent_cov2[2] = recent_cov[0] * recent_cov[2] + recent_cov[3] * recent_cov[5] + recent_cov[6] * recent_cov[8];

    recent_cov2[3] = recent_cov[1] * recent_cov[0] + recent_cov[4] * recent_cov[3] + recent_cov[7] * recent_cov[6];
    recent_cov2[4] = recent_cov[1] * recent_cov[1] + recent_cov[4] * recent_cov[4] + recent_cov[7] * recent_cov[7];
    recent_cov2[5] = recent_cov[1] * recent_cov[2] + recent_cov[4] * recent_cov[5] + recent_cov[7] * recent_cov[8];

    recent_cov2[6] = recent_cov[2] * recent_cov[0] + recent_cov[5] * recent_cov[3] + recent_cov[8] * recent_cov[6];
    recent_cov2[7] = recent_cov[2] * recent_cov[1] + recent_cov[5] * recent_cov[4] + recent_cov[8] * recent_cov[7];
    recent_cov2[8] = recent_cov[2] * recent_cov[2] + recent_cov[5] * recent_cov[5] + recent_cov[8] * recent_cov[8];
  }

  inline void eigenvalues() {
    const double SQRT3 = 1.7320508075688772;
    double b, c, d, f, g, h, i, j, K, P;
    double ma, mb, mc, me, mf, mi;
    ma = recent_cov2[0];
    mb = recent_cov2[1];
    mc = recent_cov2[2];
    me = recent_cov2[4];
    mf = recent_cov2[5];
    mi = recent_cov2[8];

    double mcmc = mc * mc;
    double mbmb = mb * mb;
    double mfmf = mf * mf;
    double mame = ma * me;
    b = ma + me + mi;
    c = mbmb + mcmc + mfmf - mame - mi * (ma + me);
    d = mi * (mame - mbmb) - ma * mfmf + 2 * (mb * mc * mf) - mcmc * me;
    double bb = b * b;
    f = -c - bb / 3;
    g = ((-2 * bb * b) - (9 * b * c)) / 27 - d;
    double gg4 = g * g / 4;
    h = gg4 + f * f * f / 27;
    i = sqrt(gg4 - h);
    j = -cbrt(i);
    K = acos(-1.0 * (g / (2 * i))) / 3;
    P = b / 3;
    double jM = j * cos(K);
    double jN = j * SQRT3 * sin(K);
    double e1 = -2 * j * cos(K) + P;
    double e2 = jM + jN + P;
    double e3 = jM - jN + P;

    recent_eigenvalues[2] = std::min(e1, std::min(e2, e3));
    recent_eigenvalues[0] = std::max(e1, std::max(e2, e3));
    recent_eigenvalues[1] =
        (recent_eigenvalues[0] == e1) ?
            (recent_eigenvalues[2] == e3 ? e2 : e3) : (recent_eigenvalues[2] == e3 ? e2 : e1);
  }

  void eigenvectors() {

    double x, y, z = 1;
    double b, c, e, f, l, li;
    int co = 0, j;

    b = recent_cov2[1];
    c = recent_cov2[2];
    e = recent_cov2[4];
    f = recent_cov2[5];
    double i = recent_cov2[8];
    double len;

    double bf = b * f - c * e;
    for (j = 0; j < 2; j++) {
      l = recent_eigenvalues[j];
      //symmetric : g=c,d=b,h=f
      double clbf = c * l + bf;
      li = l - i;
      x = (l * (li - e) + e * i - f * f) / clbf;
      y = (b * li + c * f) / clbf;
      recent_eigenvectors[co++] = x;
      recent_eigenvectors[co++] = y;
      recent_eigenvectors[co++] = z;

    }
    len = 1.0
        / sqrt(
            recent_eigenvectors[0] * recent_eigenvectors[0] + recent_eigenvectors[1] * recent_eigenvectors[1]
                + recent_eigenvectors[2] * recent_eigenvectors[2]);
    recent_eigenvectors[0] *= len;
    recent_eigenvectors[1] *= len;
    recent_eigenvectors[2] *= len;

    len = 1.0
        / sqrt(
            recent_eigenvectors[3] * recent_eigenvectors[3] + recent_eigenvectors[4] * recent_eigenvectors[4]
                + recent_eigenvectors[5] * recent_eigenvectors[5]);
    recent_eigenvectors[3] *= len;
    recent_eigenvectors[4] *= len;
    recent_eigenvectors[5] *= len;

    recent_eigenvectors[6] = recent_cov2[1] * recent_cov2[5]
        - recent_cov2[2] * (recent_cov2[4] - recent_eigenvalues[2]);
    recent_eigenvectors[7] = recent_cov2[3] * recent_cov2[2]
        - (recent_cov2[0] - recent_eigenvalues[2]) * recent_cov2[5];
    recent_eigenvectors[8] = (recent_cov2[0] - recent_eigenvalues[2]) * (recent_cov2[4] - recent_eigenvalues[2])
        - recent_cov2[3] * recent_cov2[1];

    len = 1.0
        / sqrt(
            recent_eigenvectors[6] * recent_eigenvectors[6] + recent_eigenvectors[7] * recent_eigenvectors[7]
                + recent_eigenvectors[8] * recent_eigenvectors[8]);
    recent_eigenvectors[6] *= len;
    recent_eigenvectors[7] *= len;
    recent_eigenvectors[8] *= len;

    double determinant = recent_eigenvectors[0]
        * (recent_eigenvectors[4] * recent_eigenvectors[8] - recent_eigenvectors[5] * recent_eigenvectors[7])
        - recent_eigenvectors[1]
            * (recent_eigenvectors[3] * recent_eigenvectors[8] - recent_eigenvectors[5] * recent_eigenvectors[6])
        + recent_eigenvectors[2]
            * (recent_eigenvectors[3] * recent_eigenvectors[7] - recent_eigenvectors[4] * recent_eigenvectors[6]);

    if (determinant < 0) {
      recent_eigenvectors[3] *= -1;
      recent_eigenvectors[4] *= -1;
      recent_eigenvectors[5] *= -1;
    }
  }

  void rotationMatrix() {

    double len;
    // licze macierz B: cov_matrix * wektory wlasne
    double B[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double* c = recent_cov;
    const double* e = recent_eigenvectors;

    B[0] = c[0] * e[0] + c[1] * e[1] + c[2] * e[2];
    B[1] = c[3] * e[0] + c[4] * e[1] + c[5] * e[2];
    B[2] = c[6] * e[0] + c[7] * e[1] + c[8] * e[2];

    B[3] = c[0] * e[3] + c[1] * e[4] + c[2] * e[5];
    B[4] = c[3] * e[3] + c[4] * e[4] + c[5] * e[5];
    B[5] = c[6] * e[3] + c[7] * e[4] + c[8] * e[5];

    len = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    B[0] /= len;
    B[1] /= len;
    B[2] /= len; // normalize
    len = sqrt(B[3] * B[3] + B[4] * B[4] + B[5] * B[5]);
    B[3] /= len;
    B[4] /= len;
    B[5] /= len;
    // trzeci wiersz macierzy B
    B[6] = B[1] * B[5] - B[4] * B[2];
    B[7] = -B[0] * B[5] + B[3] * B[2];
    B[8] = B[0] * B[4] - B[3] * B[1];
    // teraz macierz U = B^T*cov

    Rototranslation::rot_x(B[0] * e[0] + B[3] * e[3] + B[6] * e[6], B[0] * e[1] + B[3] * e[4] + B[6] * e[7],
        B[0] * e[2] + B[3] * e[5] + B[6] * e[8]);

    Rototranslation::rot_y(B[1] * e[0] + B[4] * e[3] + B[7] * e[6], B[1] * e[1] + B[4] * e[4] + B[7] * e[7],
        B[1] * e[2] + B[4] * e[5] + B[7] * e[8]);

    Rototranslation::rot_z(B[2] * e[0] + B[5] * e[3] + B[8] * e[6], B[2] * e[1] + B[5] * e[4] + B[8] * e[7],
        B[2] * e[2] + B[5] * e[5] + B[8] * e[8]);
  }
};

/**
 * @example ex_Crmsd.cc
 */
}
}
}
}

#endif
