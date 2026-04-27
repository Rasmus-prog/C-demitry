#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include "sfuns.h"


int main() {
    // Compare pseudo-random and quasi-random scaling on a smooth benchmark.
    auto smooth = [](const std::vector<double>& x) {
        double sum = 0.0;
        for (double value : x) {
            sum += value * value;
        }
        return std::exp(-sum);
    };
    std::vector<std::pair<double, double>> unit_cube3 = {{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
    std::vector<int> samples = {10, 100, 1000, 10000, 100000};
    const double smooth_exact = std::pow(0.5 * std::sqrt(M_PI) * std::erf(1.0), 3);
    std::ofstream scaling_out("scaling.data");
    scaling_out << "# N std_error qmc_error lcg_error std_actual qmc_actual lcg_actual\n";
    for (int num_samples : samples) {
        auto [std_integral, std_error] = monte_carlo_integrate_mt19937(smooth, unit_cube3, num_samples);
        auto [qmc_integral, qmc_error] = quasi_monte_carlo_integrate(smooth, unit_cube3, num_samples);
        auto [lcg_integral, lcg_error] = monte_carlo_integrate_lcg(smooth, unit_cube3, num_samples);
        scaling_out << num_samples << ' ' << std_error << ' ' << qmc_error << ' ' << lcg_error << ' '
                    << std::abs(std_integral - smooth_exact) << ' '
                    << std::abs(qmc_integral - smooth_exact) << ' '
                    << std::abs(lcg_integral - smooth_exact) << '\n';
    }
    scaling_out.close();

    // monte_carlo_integrate test for unit circle area as a function of the number of samples
    auto f = [](const std::vector<double>& x) {
        double r2 = x[0] * x[0] + x[1] * x[1];
        return r2 <= 1.0 ? 1.0 : 0.0;
    };
    std::vector<std::pair<double, double>> bounds = {{-1.0, 1.0}, {-1.0, 1.0}};
    samples = {1, 10, 100, 1000, 10000, 100000, 1000000};
    std::ofstream out("unit_circle.data");
    out << "# Samples Estimated_Area Error Actual" << std::endl;
    const double exact = M_PI;
    for (int num_samples : samples)
    {
    auto [integral, error] = monte_carlo_integrate(f, bounds, num_samples);
    std::cout << "Estimated area of unit circle: " << integral << " ± " << error << " with sample size: " << num_samples << std::endl;
    const double actual = std::abs(integral - exact);
    out << num_samples << " " << integral << " " << error << " " << actual << std::endl;
    }
    out.close();

    // Calculate the volume of a three-dimensional ellipsoid with semi-axes a=1, b=2, and c=3, defined by the equation
    //    x²/a²+y²/b²+z²/c² ≤ 1 
    auto f1 = [](const std::vector<double>& x) {
        double r2 = (x[0] * x[0]) / (1.0 * 1.0) + (x[1] * x[1]) / (2.0 * 2.0) + (x[2] * x[2]) / (3.0 * 3.0);
        return r2 <= 1.0 ? 1.0 : 0.0;
    };
    std::vector<std::pair<double, double>> bounds1 = {{-1.0, 1.0}, {-2.0, 2.0}, {-3.0, 3.0}};
    auto [integral1, error1] = monte_carlo_integrate(f1, bounds1, 1000000);
    std::cout << "\n";
    std::cout << "Estimated volume of the ellipsoid: " << integral1 << " ± " << error1 << std::endl;
    std::cout << "Actual volume of the ellipsoid: " << (4.0 / 3.0) * M_PI * 1.0 * 2.0 * 3.0 << std::endl;

    // Compare the requested singular integral using three generators.
    std::cout << "\n";
    std::cout << "Singular integral: ∫0π dx/π ∫0π dy/π ∫0π dz/π [1-cos(x)cos(y)cos(z)]^-1" << std::endl;
    auto f2 = [](const std::vector<double>& x) {
        return (1.0 / (M_PI * M_PI * M_PI)) /
               (1.0 - std::cos(x[0]) * std::cos(x[1]) * std::cos(x[2]));
    };
    std::vector<std::pair<double, double>> bounds2 = {{0.0, M_PI}, {0.0, M_PI}, {0.0, M_PI}};
    auto [integral_std, error_std] = monte_carlo_integrate_mt19937(f2, bounds2, 1000000);
    auto [integral_lcg, error_lcg] = monte_carlo_integrate_lcg(f2, bounds2, 1000000);
    auto [integral_qmc, error_qmc] = quasi_monte_carlo_integrate(f2, bounds2, 1000000);
    std::cout << "std::mt19937: " << integral_std << " ± " << error_std << std::endl;
    std::cout << "LCG:          " << integral_lcg << " ± " << error_lcg << std::endl;
    std::cout << "QMC:          " << integral_qmc << " ± " << error_qmc << std::endl;
    std::cout << "Reference:    " << 1.3932039296856769 << std::endl;

    
    return 0;
}
