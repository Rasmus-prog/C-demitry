#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include "sfuns.h"


int main() {
    // monte_carlo_integrate test for unit circle area as a function of the number of samples
    auto f = [](const std::vector<double>& x) {
        double r2 = x[0] * x[0] + x[1] * x[1];
        return r2 <= 1.0 ? 1.0 : 0.0;
    };
    std::vector<std::pair<double, double>> bounds = {{-1.0, 1.0}, {-1.0, 1.0}};
    std::vector<int> samples = {1, 10, 100, 1000, 10000, 100000, 1000000};
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

    // ∫0π  dx/π ∫0π  dy/π ∫0π  dz/π [1-cos(x)cos(y)cos(z)]-1
    std::cout << "\n";
    std::cout << "Estimated value of the integral: ∫0π  dx/π ∫0π  dy/π ∫0π  dz/π [1-cos(x)cos(y)cos(z)]-1" << std::endl;
    auto f2 = [](const std::vector<double>& x) {
        return (1.0 / (M_PI * M_PI * M_PI)) /
               (1.0 - std::cos(x[0]) * std::cos(x[1]) * std::cos(x[2]));
    };
    std::vector<std::pair<double, double>> bounds2 = {{0.0, M_PI}, {0.0, M_PI}, {0.0, M_PI}};
    auto [integral2, error2] = monte_carlo_integrate(f2, bounds2, 1000000);
    std::cout << "Estimated value of the integral: " << integral2 << " ± " << error2 << std::endl;
	
    auto [integral_qmc, error_qmc] = quasi_monte_carlo_integrate(f2, bounds2, 1000000);
    std::cout << "Estimated value of the integral (QMC): " << integral_qmc << " ± " << error_qmc << std::endl;

    
    return 0;
}
