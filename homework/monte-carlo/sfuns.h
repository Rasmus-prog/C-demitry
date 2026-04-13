#pragma once
#include <functional>
#include <utility>
#include <vector>
using Function = std::function<double(const std::vector<double>&)>;

// Plain Monte Carlo multi-dimensional integration.
// Returns {integral estimate, estimated standard error}.
std::pair<double, double> monte_carlo_integrate(
    const Function& f,
    const std::vector<std::pair<double, double>>& bounds,
    int num_samples);