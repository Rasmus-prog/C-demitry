#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "vector.h"

using Vec = pp::vector;
using OdeFunction = std::function<Vec(double, const Vec&)>;

double vec_norm(const Vec& v) {
	return v.norm();
}

Vec vec_add(const Vec& a, const Vec& b) {
	if (a.size() != b.size()) {
		throw std::runtime_error("Vector sizes do not match in vec_add");
	}
	Vec c(a.size());
	for (int i = 0; i < a.size(); ++i) {
		c[i] = a[i] + b[i];
	}
	return c;
}

Vec vec_sub(const Vec& a, const Vec& b) {
	if (a.size() != b.size()) {
		throw std::runtime_error("Vector sizes do not match in vec_sub");
	}
	Vec c(a.size());
	for (int i = 0; i < a.size(); ++i) {
		c[i] = a[i] - b[i];
	}
	return c;
}

Vec vec_scale(const Vec& v, double s) {
	Vec out(v.size());
	for (int i = 0; i < v.size(); ++i) {
		out[i] = v[i] * s;
	}
	return out;
}

std::tuple<Vec, Vec> rkstep12(
	const OdeFunction& f,
	double x,
	const Vec& y,
	double h
) {
	Vec k0 = f(x, y);
	Vec y_mid = vec_add(y, vec_scale(k0, h / 2.0));
	Vec k1 = f(x + h / 2.0, y_mid);

	Vec yh = vec_add(y, vec_scale(k1, h));
	Vec dy = vec_scale(vec_sub(k1, k0), h);
	return {yh, dy};
}

std::tuple<std::vector<double>, std::vector<Vec>> driver(
	const OdeFunction& f,
	double a,
	double b,
	const Vec& yinit,
	double h = 0.125,
	double acc = 0.01,
	double eps = 0.01
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
	std::vector<double> xlist{ x };
	std::vector<Vec> ylist{ y };

	int steps_taken = 0;
	const int max_steps = 1000000;

	while ((b - x) * direction > 0.0) {
		if ((x + h - b) * direction > 0.0) {
			h = b - x;
		}

		auto [yh, dy] = rkstep12(f, x, y, h);
		const double tol = (acc + eps * vec_norm(yh)) * std::sqrt(std::abs(h / (b - a)));
		const double err = vec_norm(dy);

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
		if (steps_taken > max_steps) {
			throw std::runtime_error("Too many integration steps in driver");
		}
	}

	return {xlist, ylist};
}

void write_solution(const std::string& file_name, const std::vector<double>& x, const std::vector<Vec>& y) {
	std::ofstream out(file_name);
	out << std::setprecision(16);
	for (std::size_t i = 0; i < x.size(); ++i) {
		out << x[i];
		for (int j = 0; j < y[i].size(); ++j) {
			double yi = y[i][j];
			out << ' ' << yi;
		}
		out << '\n';
	}
}

void run_simple_harmonic_oscillator() {
	OdeFunction f = [](double /*x*/, const Vec& y) {
		return Vec{y[1], -y[0]};
	};

	auto [x, y] = driver(f, 0.0, 10.0, Vec{0.0, 1.0}, 0.1, 1e-3, 1e-3);
	write_solution("out_sho.txt", x, y);

	std::cout << "Simple harmonic oscillator (u'' = -u):\n";
	std::cout << "  Saved " << x.size() << " points to out_sho.txt\n";
	std::cout << "  Sample: x, u_numeric, u_exact=sin(x), abs_error\n";
	const std::size_t stride = std::max<std::size_t>(1, x.size() / 8);
	for (std::size_t i = 0; i < x.size(); i += stride) {
		const double u_exact = std::sin(x[i]);
		const double abs_error = std::abs(y[i][0] - u_exact);
		std::cout << "  " << std::setw(9) << x[i]
				  << "  " << std::setw(12) << y[i][0]
				  << "  " << std::setw(12) << u_exact
				  << "  " << std::setw(12) << abs_error << '\n';
	}
	std::cout << '\n';
}

void run_damped_pendulum() {
	const double b = 0.25;
	const double c = 5.0;

	OdeFunction f = [b, c](double /*t*/, const Vec& y) {
		const double theta = y[0];
		const double omega = y[1];
		return Vec{omega, -b * omega - c * std::sin(theta)};
	};

	auto [t, y] = driver(f, 0.0, 10.0, Vec{std::numbers::pi_v<double> - 0.1, 0.0}, 0.05, 1e-3, 1e-3);
	write_solution("out_damped_oscillator.txt", t, y);

	std::cout << "Damped oscillator from scipy.integrate.odeint manual:\n";
	std::cout << "  Saved " << t.size() << " points to out_damped_oscillator.txt\n";
	std::cout << "  Final state: theta=" << y.back()[0] << ", omega=" << y.back()[1] << "\n\n";
}

void run_lotka_volterra() {
	const double alpha = 1.5;
	const double beta = 1.0;
	const double delta = 1.0;
	const double gamma = 3.0;

	OdeFunction f = [alpha, beta, delta, gamma](double /*t*/, const Vec& z) {
		const double x = z[0];
		const double y = z[1];
		return Vec{
			alpha * x - beta * x * y,
			delta * x * y - gamma * y
		};
	};

	auto [t, z] = driver(f, 0.0, 15.0, Vec{10.0, 5.0}, 0.05, 1e-3, 1e-3);
	write_solution("out_lotka_volterra.txt", t, z);

	std::cout << "Lotka-Volterra from scipy.integrate.solve_ivp examples:\n";
	std::cout << "  Saved " << t.size() << " points to out_lotka_volterra.txt\n";
	std::cout << "  Final state: prey=" << z.back()[0] << ", predator=" << z.back()[1] << "\n\n";
}

int main() {
	try {
		run_simple_harmonic_oscillator();
		run_damped_pendulum();
		run_lotka_volterra();
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return 1;
	}

	return 0;
}
