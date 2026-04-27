#include "sfuns.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

double rosenbrock(const std::vector<double>& x) {
    const double xx = x[0];
    const double yy = x[1];
    return (1.0 - xx) * (1.0 - xx) + 100.0 * (yy - xx * xx) * (yy - xx * xx);
}

double himmelblau(const std::vector<double>& x) {
    const double xx = x[0];
    const double yy = x[1];
    const double a = xx * xx + yy - 11.0;
    const double b = xx + yy * yy - 7.0;
    return a * a + b * b;
}

void print_vector(const std::vector<double>& x) {
    std::cout << '(';
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::cout << x[i];
        if (i + 1 != x.size()) {
            std::cout << ", ";
        }
    }
    std::cout << ')';
}

void run_case(const std::string& name, const Function& f, const std::vector<double>& x0, double acc = 1e-6) {
    const NewtonResult res = newton_minimize(f, x0, acc, 1000);

    std::cout << name << '\n';
    std::cout << "  start        = ";
    print_vector(x0);
    std::cout << '\n';
    std::cout << "  minimum      = ";
    print_vector(res.x);
    std::cout << '\n';
    std::cout << "  f(minimum)   = " << f(res.x) << '\n';
    std::cout << "  steps        = " << res.steps << '\n';
    std::cout << "  converged    = " << (res.converged ? "yes" : "no") << "\n\n";
}

} // namespace

int main() {
    std::cout << std::setprecision(12);

    run_case("Rosenbrock valley", rosenbrock, {-1.2, 1.0});
    run_case("Himmelblau function from (-2,3)", himmelblau, {-2.0, 3.0});

    return 0;
}
