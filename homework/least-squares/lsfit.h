#pragma once
#include <vector>
#include <functional>
#include <tuple>
#include "vector.h"
#include "matrix.h"

namespace pp {

// Least-squares fit of data {x_i, y_i, dy_i} with a linear combination
//   F_c(x) = sum_k c_k f_k(x)
// Returns {c, C} where c is the best-fit coefficient vector and
// C = (A^T A)^{-1} = R^{-1} R^{-T} is the covariance matrix.
// Diagonal elements give sigma_k^2 = C[k,k].
//
// The weighted design matrix is A[i,k] = f_k(x_i)/dy_i and b[i] = y_i/dy_i,
// reducing chi^2 minimisation to the overdetermined system A*c = b.
std::tuple<pp::vector, pp::matrix> lsfit(
    const std::vector<std::function<double(double)>>& fs,
    const pp::vector& x,
    const pp::vector& y,
    const pp::vector& dy
);

} // namespace pp
