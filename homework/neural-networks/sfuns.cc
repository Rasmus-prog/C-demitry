#include "sfuns.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace {

double clamp_scale(double value) {
	constexpr double minimum = 1e-3;
	if (std::abs(value) < minimum) {
		return value < 0 ? -minimum : minimum;
	}
	return value;
}

std::mt19937& rng() {
	static std::mt19937 generator{std::random_device{}()};
	return generator;
}

double uniform(double min, double max) {
	std::uniform_real_distribution<double> dist(min, max);
	return dist(rng());
}

std::vector<double> add_scaled(const std::vector<double>& a, const std::vector<double>& b, double scale) {
	std::vector<double> result(a.size());
	for (std::size_t i = 0; i < a.size(); ++i) {
		result[i] = a[i] + scale * b[i];
	}
	return result;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
	double sum = 0;
	for (std::size_t i = 0; i < a.size(); ++i) {
		sum += a[i] * b[i];
	}
	return sum;
}

double norm(const std::vector<double>& v) {
	return std::sqrt(dot(v, v));
}

} // namespace

activation_function gaussian() {
	return {
		[](double x) { return std::exp(-x * x); },
		[](double x) { return -2.0 * x * std::exp(-x * x); }
	};
}

activation_function gaussian_wavelet() {
	return {
		[](double x) { return x * std::exp(-x * x); },
		[](double x) { return (1.0 - 2.0 * x * x) * std::exp(-x * x); }
	};
}

activation_function wavelet() {
	return {
		[](double x) { return std::cos(5.0 * x) * std::exp(-x * x); },
		[](double x) {
			return (-5.0 * std::sin(5.0 * x) - 2.0 * x * std::cos(5.0 * x)) * std::exp(-x * x);
		}
	};
}

ann::ann(int n, activation_function act)
	: n(n), act(std::move(act)), p(3 * n) {
	for (int i = 0; i < n; ++i) {
		p[3 * i + 0] = uniform(-1.0, 1.0);
		p[3 * i + 1] = uniform(0.2, 1.2);
		p[3 * i + 1] = clamp_scale(p[3 * i + 1]);
		p[3 * i + 2] = uniform(-1.0, 1.0);
	}
}

double ann::response(double x) const {
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		const double a = p[3 * i + 0];
		const double b = clamp_scale(p[3 * i + 1]);
		const double w = p[3 * i + 2];
		const double z = (x - a) / b;
		sum += w * act.f(z);
	}
	return sum;
}

double ann::cost(const std::vector<double>& x, const std::vector<double>& y) const {
	double sum = 0.0;
	for (std::size_t i = 0; i < x.size(); ++i) {
		const double d = response(x[i]) - y[i];
		sum += d * d;
	}
	return sum;
}

std::vector<double> ann::gradient(const std::vector<double>& x, const std::vector<double>& y) const {
	std::vector<double> g(3 * n, 0.0);
	for (std::size_t k = 0; k < x.size(); ++k) {
		double fx = response(x[k]);
		double err = fx - y[k];

		for (int i = 0; i < n; ++i) {
			const double a = p[3 * i + 0];
			const double b = clamp_scale(p[3 * i + 1]);
			const double w = p[3 * i + 2];
			const double z = (x[k] - a) / b;
			const double fz = act.f(z);
			const double dfz = act.df(z);

			g[3 * i + 0] += 2.0 * err * w * dfz * (-1.0 / b);
			g[3 * i + 1] += 2.0 * err * w * dfz * (-(x[k] - a) / (b * b));
			g[3 * i + 2] += 2.0 * err * fz;
		}
	}
	return g;
}

// GAI start
void ann::train(const std::vector<double>& x, const std::vector<double>& y) {
	if (x.empty()) {
		return;
	}

	double current_cost = cost(x, y);
	for (int iteration = 0; iteration < 10000; ++iteration) {
		std::vector<double> g = gradient(x, y);
		const double gnorm = norm(g);
		if (gnorm < 1e-8) {
			break;
		}

		double step = 0.1;
		bool accepted = false;
		while (step > 1e-10) {
			std::vector<double> candidate = add_scaled(p, g, -step);
			for (int i = 0; i < n; ++i) {
				candidate[3 * i + 1] = clamp_scale(candidate[3 * i + 1]);
			}

			ann trial = *this;
			trial.p = std::move(candidate);
			const double trial_cost = trial.cost(x, y);
			if (trial_cost < current_cost) {
				p = std::move(trial.p);
				current_cost = trial_cost;
				accepted = true;
				break;
			}
			step *= 0.5;
		}

		if (!accepted) {
			break;
		}
	}
}
// GAI end

double tabulated_function(double x) {
	return std::cos(5.0 * x - 1.0) * std::exp(-x * x);
}
