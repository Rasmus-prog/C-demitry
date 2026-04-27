
#include "sfuns.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

double norm(const std::vector<double>& v) {
	double sum = 0.0;
	for (double value : v) {
		sum += value * value;
	}
	return std::sqrt(sum);
}

void printVector(const std::vector<double>& x) {
	std::cout << '(';
	for (std::size_t i = 0; i < x.size(); ++i) {
		std::cout << x[i];
		if (i + 1 != x.size()) {
			std::cout << ", ";
		}
	}
	std::cout << ')';
}

void reportSolution(
	const std::string& name,
	const std::function<std::vector<double>(const std::vector<double>&)>& f,
	const std::vector<double>& x0,
	const std::vector<double>& dx = {}
) {
	const std::vector<double> root = newton(f, x0, 1e-10, 1e-3, 200, dx);
	std::cout << name << "\n  start = ";
	printVector(x0);
	std::cout << "\n  root  = ";
	printVector(root);
	std::cout << "\n  ||f(root)|| = " << norm(f(root)) << "\n\n";
}

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

} // namespace

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
