
#include "sfuns.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main() {
	std::cout << std::setprecision(12);

	// Simple 1D debug: root of x^2-2 = 0.
	const auto f1 = [](const std::vector<double>& x) {
		return std::vector<double>{x[0] * x[0] - 2.0};
	};
	reportSolution("Simple 1D test: x^2-2", f1, {2.0}, {1e-8});

	// Simple 2D debug: circle-line intersection.
	const auto f2 = [](const std::vector<double>& x) {
		return std::vector<double>{
			x[0] * x[0] + x[1] * x[1] - 1.0,
			x[0] - x[1]
		};
	};
	reportSolution("Simple 2D test: x^2+y^2=1, x=y", f2, {1.0, 0.5});

	// Rosenbrock: solve gradient = 0.
	const auto rosenbrockGrad = [](const std::vector<double>& x) {
		const double xx = x[0];
		const double yy = x[1];
		return std::vector<double>{
			2.0 * (xx - 1.0) - 400.0 * xx * (yy - xx * xx),
			200.0 * (yy - xx * xx)
		};
	};

	const std::vector<double> rb = newton(rosenbrockGrad, {-1.2, 1.0}, 1e-10, 1e-3, 200);
	std::cout << "Rosenbrock stationary point from (-1.2,1.0)\n  point = ";
	printVector(rb);
	std::cout << "\n  ||grad|| = " << norm(rosenbrockGrad(rb));
	std::cout << "\n  f(point) = " << rosenbrock(rb) << "\n\n";

	// Himmelblau: solve gradient = 0 from different starting points.
	const auto himmelblauGrad = [](const std::vector<double>& x) {
		const double xx = x[0];
		const double yy = x[1];
		const double a = xx * xx + yy - 11.0;
		const double b = xx + yy * yy - 7.0;
		return std::vector<double>{
			4.0 * xx * a + 2.0 * b,
			2.0 * a + 4.0 * yy * b
		};
	};

	const std::vector<std::vector<double>> starts = {
		{3.0, 2.0},
		{-2.0, 3.0},
		{-4.0, -4.0},
		{4.0, -2.0}
	};

	for (const auto& x0 : starts) {
		const std::vector<double> p = newton(himmelblauGrad, x0, 1e-10, 1e-3, 200);
		std::cout << "Himmelblau stationary point from ";
		printVector(x0);
		std::cout << "\n  point = ";
		printVector(p);
		std::cout << "\n  ||grad|| = " << norm(himmelblauGrad(p));
		std::cout << "\n  f(point) = " << himmelblau(p) << "\n\n";
	}

	return 0;
}
