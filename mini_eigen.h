// mini_eigen.h -- self-contained replacement for the slice of Eigen this
// project uses: dynamic MatrixXd / VectorXd and least-squares solves on
// small dense systems (nothing here exceeds 5x4).
//
// API-compatible with the Eigen spellings used in the .cpp files:
//   MatrixXd m(r,c); m(i,j); m*x
//   VectorXd b(n); b(i); b.size(); b.begin(); b.end()
//   m.colPivHouseholderQr().solve(b)
//   m.fullPivHouseholderQr().solve(b)
//   m.lu().solve(b)
// All three solvers are backed by one column-pivoted Householder QR:
// exact solution for square full-rank systems (what lu() is used for),
// least-squares solution for overdetermined ones, and a rank-revealing
// basic solution when rank-deficient (callers verify with a residual
// check, e.g. b_prime = m*x). Zlib-style freedom: no external deps.
#pragma once

#include <cmath>
#include <vector>

namespace Eigen {

class VectorXd {
    std::vector<double> d;
public:
    VectorXd() {}
    explicit VectorXd(int n) : d(n, 0.0) {}
    double &operator()(int i) { return d[i]; }
    const double &operator()(int i) const { return d[i]; }
    int size() const { return (int)d.size(); }
    double *begin() { return d.data(); }
    double *end() { return d.data() + d.size(); }
    const double *begin() const { return d.data(); }
    const double *end() const { return d.data() + d.size(); }
};

class MatrixXd;

// Column-pivoted Householder QR of an m x n matrix (m >= n for the
// systems in this codebase). solve() returns the least-squares /
// basic solution, matching Eigen's ColPivHouseholderQR behavior.
class ColPivHouseholderQR {
    int m, n, rank;
    std::vector<double> qr;      // column-major; Householder vectors below R
    std::vector<double> beta;    // Householder scalars
    std::vector<int> perm;       // column permutation
    double &at(int i, int j) { return qr[j * m + i]; }
    double atc(int i, int j) const { return qr[j * m + i]; }
public:
    ColPivHouseholderQR(const MatrixXd &A);

    VectorXd solve(const VectorXd &b) const {
        // c = Q^T b via stored Householder reflectors
        std::vector<double> c(b.begin(), b.end());
        for (int k = 0; k < rank; k++) {
            double dot = c[k];
            for (int i = k + 1; i < m; i++) dot += atc(i, k) * c[i];
            dot *= beta[k];
            c[k] -= dot;
            for (int i = k + 1; i < m; i++) c[i] -= atc(i, k) * dot;
        }
        // back-substitute R(0:rank,0:rank) y = c(0:rank)
        std::vector<double> y(n, 0.0);
        for (int k = rank - 1; k >= 0; k--) {
            double s = c[k];
            for (int j = k + 1; j < rank; j++) s -= atc(k, j) * y[j];
            y[k] = s / atc(k, k);
        }
        VectorXd x(n);
        for (int j = 0; j < n; j++) x(perm[j]) = y[j];
        return x;
    }
};

class MatrixXd {
    int r_, c_;
    std::vector<double> d;       // column-major
public:
    MatrixXd() : r_(0), c_(0) {}
    MatrixXd(int r, int c) : r_(r), c_(c), d((size_t)r * c, 0.0) {}
    double &operator()(int i, int j) { return d[(size_t)j * r_ + i]; }
    const double &operator()(int i, int j) const { return d[(size_t)j * r_ + i]; }
    int rows() const { return r_; }
    int cols() const { return c_; }

    VectorXd operator*(const VectorXd &x) const {
        VectorXd out(r_);
        for (int j = 0; j < c_; j++)
            for (int i = 0; i < r_; i++)
                out(i) += (*this)(i, j) * x(j);
        return out;
    }

    ColPivHouseholderQR colPivHouseholderQr() const { return ColPivHouseholderQR(*this); }
    ColPivHouseholderQR fullPivHouseholderQr() const { return ColPivHouseholderQR(*this); }
    ColPivHouseholderQR lu() const { return ColPivHouseholderQR(*this); }
};

inline ColPivHouseholderQR::ColPivHouseholderQR(const MatrixXd &A)
    : m(A.rows()), n(A.cols()), rank(0), qr((size_t)A.rows() * A.cols()),
      beta(std::min(A.rows(), A.cols()), 0.0), perm(A.cols()) {
    for (int j = 0; j < n; j++) {
        perm[j] = j;
        for (int i = 0; i < m; i++) at(i, j) = A(i, j);
    }
    int steps = std::min(m, n);
    double max_pivot = 0.0;
    for (int k = 0; k < steps; k++) {
        // pivot: bring the column with the largest remaining norm to k
        int best = k;
        double best_norm = -1.0;
        for (int j = k; j < n; j++) {
            double s = 0.0;
            for (int i = k; i < m; i++) s += atc(i, j) * atc(i, j);
            if (s > best_norm) { best_norm = s; best = j; }
        }
        if (best != k) {
            for (int i = 0; i < m; i++) std::swap(at(i, k), at(i, best));
            std::swap(perm[k], perm[best]);
        }
        // Householder reflector zeroing column k below the diagonal
        double norm = std::sqrt(best_norm);
        if (norm == 0.0) { beta[k] = 0.0; continue; }
        double alpha = atc(k, k) >= 0 ? -norm : norm;
        double v0 = atc(k, k) - alpha;
        beta[k] = -v0 / alpha;
        for (int i = k + 1; i < m; i++) at(i, k) /= v0;
        at(k, k) = alpha;
        for (int j = k + 1; j < n; j++) {
            double dot = atc(k, j);
            for (int i = k + 1; i < m; i++) dot += atc(i, k) * atc(i, j);
            dot *= beta[k];
            at(k, j) -= dot;
            for (int i = k + 1; i < m; i++) at(i, j) -= atc(i, k) * dot;
        }
        if (std::fabs(alpha) > max_pivot) max_pivot = std::fabs(alpha);
    }
    // rank: pivots above Eigen's default threshold (eps * diag size,
    // relative to the largest pivot)
    double thresh = max_pivot * steps * 2.220446049250313e-16;
    for (int k = 0; k < steps; k++)
        if (std::fabs(atc(k, k)) > thresh) rank++;
        else break;
}

}  // namespace Eigen
