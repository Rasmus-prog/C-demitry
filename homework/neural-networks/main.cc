#include "sfuns.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
	constexpr int sample_count = 21;
	constexpr double left = -1.0;
	constexpr double right = 1.0;

	std::vector<double> x(sample_count);
	std::vector<double> y(sample_count);

	for (int i = 0; i < sample_count; ++i) {
		const double t = static_cast<double>(i) / (sample_count - 1);
		x[i] = left + t * (right - left);
		y[i] = tabulated_function(x[i]);
	}

	ann network(8, gaussian_wavelet());
	network.train(x, y);

	std::cout << std::fixed << std::setprecision(8);
	std::cout << "# x target response\n";
	for (int i = 0; i <= 100; ++i) {
		const double t = static_cast<double>(i) / 100.0;
		const double xi = left + t * (right - left);
		std::cout << xi << ' ' << tabulated_function(xi) << ' ' << network.response(xi) << '\n';
	}

	std::cout << "# cost " << network.cost(x, y) << '\n';
	return 0;
}
