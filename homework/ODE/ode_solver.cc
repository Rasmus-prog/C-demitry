#include "ode_solver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ode {

std::tuple<Vec, Vec> rkstep12(const OdeFunction& f, double x, const Vec& y, double h) {
	const Vec k0 = f(x, y);
	const Vec y_mid = y + k0 * (h / 2.0);
	const Vec k1 = f(x + h / 2.0, y_mid);

	const Vec yh = y + k1 * h;
	const Vec dy = (k1 - k0) * h;
	return {yh, dy};
}

std::tuple<std::vector<double>, std::vector<Vec>> driver(
	const OdeFunction& f,
	double a,
	double b,
	const Vec& yinit,
	double h,
	double acc,
	double eps
) {
	if (a == b) {
		return {{a}, {yinit}};
	}

	const double direction = (b > a) ? 1.0 : -1.0;
	if (h == 0.0) {
		h = 0.125 * direction;
	}
	if (h * direction < 0.0) {
		h = -h;
	}

	double x = a;
	Vec y = yinit;
	std::vector<double> xlist{x};
	std::vector<Vec> ylist{y};

	int steps_taken = 0;

	while ((b - x) * direction > 0.0) {
		if ((x + h - b) * direction > 0.0) {
			h = b - x;
		}

		auto [yh, dy] = rkstep12(f, x, y, h);
		const double tol = (acc + eps * yh.norm()) * std::sqrt(std::abs(h / (b - a)));
		const double err = dy.norm();

		if (err <= tol) {
			x += h;
			y = yh;
			xlist.push_back(x);
			ylist.push_back(y);
		}

		if (err > 0.0) {
			const double factor = std::min(0.95 * std::pow(tol / err, 0.25), 2.0);
			h *= std::max(0.1, factor);
		} else {
			h *= 2.0;
		}

		++steps_taken;
	}

	return {xlist, ylist};
}

} // namespace ode
