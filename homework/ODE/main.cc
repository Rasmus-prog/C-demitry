#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ode_solver.h"
#include "vector.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

void write_solution(const std::string& file_name, const std::vector<double>& x, const std::vector<ode::Vec>& y) {
	std::ofstream out(file_name);
	out << std::setprecision(16);
	for (std::size_t i = 0; i < x.size(); ++i) {
		out << x[i];
		for (int j = 0; j < y[i].size(); ++j) {
			out << ' ' << y[i][j];
		}
		out << '\n';
	}
}

void run_simple_harmonic_oscillator() {
	ode::OdeFunction f = [](double /*x*/, const ode::Vec& y) {
		return ode::Vec{y[1], -y[0]};
	};

	auto [x, y] = ode::driver(f, 0.0, 10.0, ode::Vec{0.0, 1.0}, 0.1, 1e-3, 1e-3);
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

	ode::OdeFunction f = [b, c](double /*t*/, const ode::Vec& y) {
		const double theta = y[0];
		const double omega = y[1];
		return ode::Vec{omega, -b * omega - c * std::sin(theta)};
	};

	auto [t, y] = ode::driver(f, 0.0, 10.0, ode::Vec{kPi - 0.1, 0.0}, 0.05, 1e-3, 1e-3);
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

	ode::OdeFunction f = [alpha, beta, delta, gamma](double /*t*/, const ode::Vec& z) {
		const double prey = z[0];
		const double predator = z[1];
		return ode::Vec{
			alpha * prey - beta * prey * predator,
			delta * prey * predator - gamma * predator
		};
	};

	auto [t, z] = ode::driver(f, 0.0, 15.0, ode::Vec{10.0, 5.0}, 0.05, 1e-3, 1e-3);
	write_solution("out_lotka_volterra.txt", t, z);

	std::cout << "Lotka-Volterra from scipy.integrate.solve_ivp examples:\n";
	std::cout << "  Saved " << t.size() << " points to out_lotka_volterra.txt\n";
	std::cout << "  Final state: prey=" << z.back()[0] << ", predator=" << z.back()[1] << "\n\n";
}

void run_relativistic_orbits() {
	auto simulate_orbit = [](double epsilon, double u0, double up0, const std::string& file_name) {
		ode::OdeFunction f = [epsilon](double /*phi*/, const ode::Vec& y) {
			const double u = y[0];
			const double up = y[1];
			return ode::Vec{up, 1.0 - u + epsilon * u * u};
		};

		const double rotations = 10.0;
		const double phi_max = rotations * 2.0 * kPi;
		auto [phi, y] = ode::driver(f, 0.0, phi_max, ode::Vec{u0, up0}, 0.02, 1e-5, 1e-5);
		write_solution(file_name, phi, y);

		std::cout << "  Saved " << phi.size() << " points to " << file_name << "\n";
		std::cout << "  Final state: u=" << y.back()[0] << ", u'=" << y.back()[1] << "\n";
	};

	std::cout << "Relativistic planetary orbit (u'' + u = 1 + eps*u^2):\n";
	std::cout << "  Case 1: Newtonian circular orbit (eps=0, u(0)=1, u'(0)=0)\n";
	simulate_orbit(0.0, 1.0, 0.0, "out_orbit_circular.txt");

	std::cout << "  Case 2: Newtonian elliptical orbit (eps=0, u(0)=1, u'(0)=-0.5)\n";
	simulate_orbit(0.0, 1.0, -0.5, "out_orbit_newtonian_elliptic.txt");

	std::cout << "  Case 3: Relativistic precession (eps=0.01, u(0)=1, u'(0)=-0.5)\n";
	simulate_orbit(0.01, 1.0, -0.5, "out_orbit_relativistic_precession.txt");
	std::cout << '\n';
}

} // namespace

int main() {
		run_simple_harmonic_oscillator();
		run_damped_pendulum();
		run_lotka_volterra();
		run_relativistic_orbits();
	return 0;
}
