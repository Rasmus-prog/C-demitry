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
	write_solution("sho.data", x, y);

	std::cout << "Simple harmonic oscillator (u'' = -u):\n";
	std::cout << "  Saved " << x.size() << " points to sho.data\n";
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
	write_solution("damped_oscillator.data", t, y);

	std::cout << "Damped oscillator from scipy.integrate.odeint manual:\n";
	std::cout << "  Saved " << t.size() << " points to damped_oscillator.data\n";
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
	write_solution("lotka_volterra.data", t, z);

	std::cout << "Lotka-Volterra from scipy.integrate.solve_ivp examples:\n";
	std::cout << "  Saved " << t.size() << " points to lotka_volterra.data\n";
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
	simulate_orbit(0.0, 1.0, 0.0, "orbit_circular.data");

	std::cout << "  Case 2: Newtonian elliptical orbit (eps=0, u(0)=1, u'(0)=-0.5)\n";
	simulate_orbit(0.0, 1.0, -0.5, "orbit_newtonian_elliptic.data");

	std::cout << "  Case 3: Relativistic precession (eps=0.01, u(0)=1, u'(0)=-0.5)\n";
	simulate_orbit(0.01, 1.0, -0.5, "orbit_relativistic_precession.data");
	std::cout << '\n';
}

void run_three_body_figure8() {
	// Equal-mass Newtonian three-body figure-8 initial conditions.
	const ode::Vec z0{
		 0.4662036850,  0.4323657300,
		 0.4662036850,  0.4323657300,
		-0.9324073700, -0.8647314600,
		 0.9700043600, -0.2430875300,
		-0.9700043600,  0.2430875300,
		 0.0,           0.0
	};

	ode::OdeFunction f = [](double /*t*/, const ode::Vec& z) {
		auto ax_from = [](double xi, double yi, double xj, double yj) {
			const double dx = xj - xi;
			const double dy = yj - yi;
			const double r2 = dx * dx + dy * dy;
			const double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
			return std::pair<double, double>{dx * inv_r3, dy * inv_r3};
		};

		const double x1 = z[6], y1 = z[7];
		const double x2 = z[8], y2 = z[9];
		const double x3 = z[10], y3 = z[11];

		const auto [a12x, a12y] = ax_from(x1, y1, x2, y2);
		const auto [a13x, a13y] = ax_from(x1, y1, x3, y3);
		const auto [a21x, a21y] = ax_from(x2, y2, x1, y1);
		const auto [a23x, a23y] = ax_from(x2, y2, x3, y3);
		const auto [a31x, a31y] = ax_from(x3, y3, x1, y1);
		const auto [a32x, a32y] = ax_from(x3, y3, x2, y2);

		return ode::Vec{
			a12x + a13x, a12y + a13y,
			a21x + a23x, a21y + a23y,
			a31x + a32x, a31y + a32y,
			z[0], z[1],
			z[2], z[3],
			z[4], z[5]
		};
	};

	const double period = 6.32591398;

	auto [t, z] = ode::driver(f, 0.0, 3.0 * period, z0, 1e-2, 1e-4, 1e-4);
	write_solution("three_body_figure8.data", t, z);

	std::cout << "Newtonian three-body figure-8 orbit (m_i=1, G=1):\n";
	std::cout << "  Saved " << t.size() << " points to three_body_figure8.data\n";
	std::cout << "  Final positions:"
	          << " r1=(" << z.back()[6] << ", " << z.back()[7] << ")"
	          << " r2=(" << z.back()[8] << ", " << z.back()[9] << ")"
	          << " r3=(" << z.back()[10] << ", " << z.back()[11] << ")\n\n";
}

} // namespace

int main() {
		run_simple_harmonic_oscillator();
		run_damped_pendulum();
		run_lotka_volterra();
		run_relativistic_orbits();
		run_three_body_figure8();
	return 0;
}
