#pragma once

#include <functional>
#include <limits>
#include <utility>

using Function = std::function<double(double)>;

// Adaptive 1D integrator with random (stratified) abscissas.
//
// At each recursion level N new points are drawn uniformly at random within
// the subinterval.  Statistics inherited from the parent (sum, count) are
// reused so no re-evaluation is needed.
//
// Error estimate: |mean_new - mean_inherited| * (b-a)
// When there is no inherited mean (first call) the sample variance is used.
//
// N is the number of NEW random points thrown at each level (default 4).

std::pair<double, double> integrate_random(
    const Function& f,
    double a,
    double b,
    double acc   = 1e-3,
    double eps   = 1e-3,
    int    N     = 4,
    // inherited statistics from the parent split (nullopt = first call)
    double inherited_mean = std::numeric_limits<double>::quiet_NaN(),
    double inherited_count = 0.0
);

// Convenience wrapper: Clenshaw-Curtis transform applied before random integration.
// Handles endpoint singularities much better than plain random abscissas.
std::pair<double, double> integrate_random_cc(
    const Function& f,
    double a,
    double b,
    double acc = 1e-3,
    double eps = 1e-3,
    int    N   = 4
);
