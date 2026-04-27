#pragma once

#include <functional>
#include <vector>

using Function = std::function<double(const std::vector<double>&)>;

struct NewtonResult {
    std::vector<double> x;
    int steps = 0;
    bool converged = false;
};

std::vector<double> gradient(const Function& phi, const std::vector<double>& x);
std::vector<std::vector<double>> hessian(const Function& phi, const std::vector<double>& x);
NewtonResult newton_minimize(
    const Function& phi,
    const std::vector<double>& x0,
    double acc = 1e-3,
    int max_steps = 1000
);
