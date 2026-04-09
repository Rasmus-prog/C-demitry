#include "sfuns.h"

#include <cmath>

double integrate_finite(const Function& f,
						double a,
						double b,
						double acc = 1e-3,
						double eps = 1e-3,
						std::optional<double> f2 = std::nullopt,
						std::optional<double> f3 = std::nullopt) {
	const double h = b - a;

	const double x1 = a + h / 6.0;
	const double x2 = a + 2.0 * h / 6.0;
	const double x3 = a + 4.0 * h / 6.0;
	const double x4 = a + 5.0 * h / 6.0;

	const double v1 = f(x1);
	const double v2 = f2.has_value() ? *f2 : f(x2);
	const double v3 = f3.has_value() ? *f3 : f(x3);
	const double v4 = f(x4);

	const double Q = (2 * v1 + v2 + v3 + 2 * v4) / 6 * (b - a);
	const double q = (v1 + v2 + v3 + v4) / 4 * (b - a);

	const double err = std::abs(Q - q);
	const double tol = acc + eps * std::abs(Q);

	if (err < tol) {
		return Q;
	}

	const double m = (a + b) / 2.0;
	const double nextAcc = acc / std::sqrt(2.0);

	return integrate_finite(f, a, m, nextAcc, eps, v1, v2) +
		   integrate_finite(f, m, b, nextAcc, eps, v3, v4);
}

double integrate(const Function& f,
				 double a,
				 double b,
				 double acc,
				 double eps,
				 std::optional<double> f2,
				 std::optional<double> f3) {
	if (std::isfinite(a) && std::isfinite(b)) {
		return integrate_finite(f, a, b, acc, eps, f2, f3);
	}

	if (std::isinf(a) && std::isinf(b)) {
		auto right_half = [f](double t) {
			const double x = t / (1.0 - t);
			const double dxdt = 1.0 / ((1.0 - t) * (1.0 - t));
			return f(x) * dxdt;
		};
		auto left_half = [f](double t) {
			const double x = -t / (1.0 - t);
			const double dxdt = 1.0 / ((1.0 - t) * (1.0 - t));
			return f(x) * dxdt;
		};
		return integrate_clenshaw_curtis(left_half, 0.0, 1.0, acc / std::sqrt(2.0), eps) +
			   integrate_clenshaw_curtis(right_half, 0.0, 1.0, acc / std::sqrt(2.0), eps);
	}

	if (std::isinf(b)) {
		auto transformed = [f, a](double t) {
			const double x = a + t / (1.0 - t);
			const double dxdt = 1.0 / ((1.0 - t) * (1.0 - t));
			return f(x) * dxdt;
		};
		return integrate_clenshaw_curtis(transformed, 0.0, 1.0, acc, eps);
	}

	auto transformed = [f, b](double t) {
		const double x = b - t / (1.0 - t);
		const double dxdt = 1.0 / ((1.0 - t) * (1.0 - t));
		return f(x) * dxdt;
	};
	return integrate_clenshaw_curtis(transformed, 0.0, 1.0, acc, eps);
}

double errf_func(double x, double acc, double eps) {
	if (x < 0.0) {
		return -errf_func(-x, acc, eps);
	} else if (x <= 1.0) {
		auto f = [](double t) { return std::exp(-t * t); };
		return 2.0 / std::sqrt(M_PI) * integrate(f, 0.0, x, acc, eps);
	} else {
		auto f = [x](double t) { return std::exp(-(x + (1.0 - t) / t) * (x + (1.0 - t) / t)) / t / t; };
		return 1.0 - 2.0 / std::sqrt(M_PI) * integrate(f, 0.0, 1.0, acc, eps);
	}
}

// Inplement an (open quandrature) adaptive integrator with the Clenshaw-Curtis variable transformation,
double integrate_clenshaw_curtis(const Function& f,
								 double a,
								 double b,
								 double acc,
								 double eps) {
	auto transformed_f = [f, a, b](double t) {
		const double x = (a + b) / 2.0 + (b - a) / 2.0 * std::cos(M_PI * t);
		return f(x) * (b - a) / 2.0 * M_PI * std::sin(M_PI * t);
	};
	return integrate(transformed_f, 0.0, 1.0, acc, eps);
}
