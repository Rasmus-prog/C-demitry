#include "lsfit.h"
#include "matrix.h"
#include "qr.h"
#include <cassert>
#include <cmath>

namespace pp {

std::tuple<pp::vector, pp::matrix> lsfit(
    const std::vector<std::function<double(double)>>& fs,
    const pp::vector& x,
    const pp::vector& y,
    const pp::vector& dy)
{
    int n = x.size();
    int m = (int)fs.size();
    assert((int)y.size() == n && (int)dy.size() == n && n >= m);

    // Build weighted design matrix A[i,k] = f_k(x_i)/dy_i  and  b[i] = y_i/dy_i
    pp::matrix A(n, m);
    pp::vector b(n);
    for (int i = 0; i < n; i++) {
        b[i] = y[i] / dy[i];
        for (int k = 0; k < m; k++) {
            A[i, k] = fs[k](x[i]) / dy[i];
        }
    }

    // QR-decompose A (tall matrix: n x m, n >= m) and solve normal equations
    pp::QR qr(A);
    pp::vector c = qr.solve(b);

    // Covariance matrix C = (A^T A)^{-1} = R^{-1} (R^{-1})^T
    // Compute R_inv column by column via back-substitution: R * x = e_j
    const pp::matrix& R = qr.R;
    pp::matrix R_inv(m, m);
    for (int j = 0; j < m; j++) {
        for (int i = m - 1; i >= 0; i--) {
            double s = (i == j) ? 1.0 : 0.0;
            for (int k = i + 1; k < m; k++) s -= R[i, k] * R_inv[k, j];
            R_inv[i, j] = s / R[i, i];
        }
    }
    // C = R_inv * R_inv^T
    pp::matrix C(m, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            double s = 0;
            for (int k = 0; k < m; k++) s += R_inv[i, k] * R_inv[j, k];
            C[i, j] = s;
        }

    return {c, C};
}

} // namespace pp
