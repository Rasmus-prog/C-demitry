#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "jacobi.h"

struct Eigensystem {
	Vector w;
	Matrix V;
};

Matrix build_hydrogen_hamiltonian(double rmax, double dr, Vector& r_grid) {
	const int npoints = static_cast<int>(rmax / dr) - 1;
	if (npoints < 2) throw std::runtime_error("Too few grid points: increase rmax or decrease dr.");

	r_grid.assign(static_cast<size_t>(npoints), 0.0);
	for (int i = 0; i < npoints; ++i) r_grid[static_cast<size_t>(i)] = dr * (i + 1);

	Matrix H(npoints, npoints);
	const double tdiag = -2.0 * (-0.5 / (dr * dr));
	const double toff = 1.0 * (-0.5 / (dr * dr));

	for (int i = 0; i < npoints - 1; ++i) {
		H(i, i) = tdiag;
		H(i, i + 1) = toff;
		H(i + 1, i) = toff;
	}
	H(npoints - 1, npoints - 1) = tdiag;

	for (int i = 0; i < npoints; ++i) H(i, i) += -1.0 / r_grid[static_cast<size_t>(i)];
	return H;
}

Eigensystem sorted_eigensystem(const Matrix& H) {
	auto [w, V] = jacobi_cyclic(H);
	const int n = static_cast<int>(w.size());
	std::vector<int> idx(static_cast<size_t>(n));
	std::iota(idx.begin(), idx.end(), 0);
	std::sort(idx.begin(), idx.end(), [&](int a, int b) { return w[static_cast<size_t>(a)] < w[static_cast<size_t>(b)]; });

	Vector w_sorted(static_cast<size_t>(n));
	Matrix V_sorted(n, n);
	for (int k = 0; k < n; ++k) {
		const int old_col = idx[static_cast<size_t>(k)];
		w_sorted[static_cast<size_t>(k)] = w[static_cast<size_t>(old_col)];
		for (int i = 0; i < n; ++i) V_sorted(i, k) = V(i, old_col);
	}

	return {w_sorted, V_sorted};
}

double exact_energy_s(int n) {
	return -0.5 / (static_cast<double>(n) * static_cast<double>(n));
}

double assoc_laguerre_alpha1(int m, double x) {
	if (m == 0) return 1.0;
	if (m == 1) return 2.0 - x;
	double Lm2 = 1.0;
	double Lm1 = 2.0 - x;
	for (int k = 2; k <= m; ++k) {
		const double kk = static_cast<double>(k);
		const double Lk = ((2.0 * kk - x) * Lm1 - kk * Lm2) / kk;
		Lm2 = Lm1;
		Lm1 = Lk;
	}
	return Lm1;
}

double exact_reduced_radial_s(int n, double r) {
	const double nn = static_cast<double>(n);
	const double x = 2.0 * r / nn;
	const double pref = 2.0 / std::pow(nn, 2.5);
	return r * pref * std::exp(-r / nn) * assoc_laguerre_alpha1(n - 1, x);
}

void write_energies(const std::string& filename, const Vector& w, int nlevels) {
	std::ofstream out(filename);
	out << "# n  E_numerical  E_exact  abs_error\n";
	for (int k = 0; k < nlevels; ++k) {
		const int n = k + 1;
		const double e_exact = exact_energy_s(n);
		const double e_num = w[static_cast<size_t>(k)];
		out << n << ' ' << std::setprecision(16) << e_num << ' ' << e_exact << ' ' << std::abs(e_num - e_exact) << '\n';
	}
}

void write_convergence_dr(const std::string& filename, double rmax) {
	const std::vector<double> dr_values = {0.60, 0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10};
	std::ofstream out(filename);
	out << "# dr  E0_numerical  E0_exact  abs_error\n";
	for (double dr : dr_values) {
		if (static_cast<int>(rmax / dr) - 1 < 2) continue;
		Vector r;
		const Matrix H = build_hydrogen_hamiltonian(rmax, dr, r);
		const Eigensystem evd = sorted_eigensystem(H);
		const double e0 = evd.w[0];
		const double e_exact = exact_energy_s(1);
		out << std::setprecision(16) << dr << ' ' << e0 << ' ' << e_exact << ' ' << std::abs(e0 - e_exact) << '\n';
	}
}

void write_convergence_rmax(const std::string& filename, double dr) {
	const std::vector<double> rmax_values = {4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 22.0, 26.0, 30.0};
	std::ofstream out(filename);
	out << "# rmax  E0_numerical  E0_exact  abs_error\n";
	for (double rmax : rmax_values) {
		if (static_cast<int>(rmax / dr) - 1 < 2) continue;
		Vector r;
		const Matrix H = build_hydrogen_hamiltonian(rmax, dr, r);
		const Eigensystem evd = sorted_eigensystem(H);
		const double e0 = evd.w[0];
		const double e_exact = exact_energy_s(1);
		out << std::setprecision(16) << rmax << ' ' << e0 << ' ' << e_exact << ' ' << std::abs(e0 - e_exact) << '\n';
	}
}

void write_wavefunctions(const std::string& filename, const Vector& r, const Matrix& V, double dr, int nstates) {
	const int npoints = static_cast<int>(r.size());
	const double scale = 1.0 / std::sqrt(dr);
	std::ofstream out(filename);
	out << "# r";
	for (int n = 1; n <= nstates; ++n) out << " f_num_n=" << n << " f_exact_n=" << n;
	out << "\n";

	for (int i = 0; i < npoints; ++i) {
		const double rr = r[static_cast<size_t>(i)];
		out << std::setprecision(16) << rr;
		for (int n = 1; n <= nstates; ++n) {
			double f_num = scale * V(i, n - 1);
			const double f_exact = exact_reduced_radial_s(n, rr);
			if (f_num * f_exact < 0) f_num = -f_num;
			out << ' ' << f_num << ' ' << f_exact;
		}
		out << '\n';
	}
}

Matrix random_symmetric_matrix(int n, std::mt19937_64& rng) {
	std::uniform_real_distribution<double> dist(-1.0, 1.0);
	Matrix A(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			const double x = dist(rng);
			A(i, j) = x;
			A(j, i) = x;
		}
	}
	return A;
}

int main(int argc, char** argv) {
	double rmax = 10.0;
	double dr = 0.30;
	int benchmark_n = -1;
	bool ground_only = false;

	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "-rmax" && i + 1 < argc) {
			rmax = std::stod(argv[++i]);
		} else if (arg == "-dr" && i + 1 < argc) {
			dr = std::stod(argv[++i]);
		} else if (arg == "-benchmark" && i + 1 < argc) {
			benchmark_n = std::stoi(argv[++i]);
		} else if (arg == "-ground") {
			ground_only = true;
		} else {
			std::cerr << "Usage: ./main [-rmax value] [-dr value] [-ground] [-benchmark N]\n";
			return 1;
		}
	}

	if (benchmark_n > 0) {
		std::mt19937_64 rng(20260220 + static_cast<unsigned long long>(benchmark_n));
		const Matrix A = random_symmetric_matrix(benchmark_n, rng);
		const auto t0 = std::chrono::steady_clock::now();
		const auto [w, V] = jacobi_cyclic(A);
		const auto t1 = std::chrono::steady_clock::now();
		const std::chrono::duration<double> dt = t1 - t0;
		double trace_sum = 0.0;
		for (double x : w) trace_sum += x;
		std::cout << benchmark_n << ' ' << std::setprecision(16) << dt.count() << ' ' << trace_sum << '\n';
		(void)V;
		return 0;
	}

	if (rmax <= 0 || dr <= 0) {
		std::cerr << "Both rmax and dr must be positive.\n";
		return 1;
	}

	const int npoints = static_cast<int>(rmax / dr) - 1;
	if (npoints < 2) {
		std::cerr << "Need at least 2 grid points. Try larger rmax or smaller dr.\n";
		return 1;
	}

	Vector r;
	const Matrix H = build_hydrogen_hamiltonian(rmax, dr, r);
	const Eigensystem evd = sorted_eigensystem(H);
	if (ground_only) {
		const double e0 = evd.w[0];
		const double e_exact = exact_energy_s(1);
		std::cout << std::setprecision(16) << rmax << ' ' << dr << ' ' << e0 << ' ' << e_exact << ' ' << std::abs(e0 - e_exact)
			      << '\n';
		return 0;
	}

	const int nlevels = std::min(5, static_cast<int>(evd.w.size()));
	write_energies("hydrogen_energies.txt", evd.w, nlevels);
	write_convergence_dr("hydrogen_convergence_dr.txt", rmax);
	write_convergence_rmax("hydrogen_convergence_rmax.txt", dr);
	write_wavefunctions("hydrogen_wavefunctions.txt", r, evd.V, dr, std::min(3, nlevels));

	std::cout << std::setprecision(8) << std::fixed;
	std::cout << "Hydrogen s-wave radial equation solved on grid\n";
	std::cout << "rmax=" << rmax << ", dr=" << dr << ", npoints=" << npoints << "\n";
	std::cout << "Lowest energies (numerical vs exact):\n";
	for (int k = 0; k < nlevels; ++k) {
		const int n = k + 1;
		const double e_num = evd.w[static_cast<size_t>(k)];
		const double e_exact = exact_energy_s(n);
		std::cout << "n=" << n << "  E_num=" << e_num << "  E_exact=" << e_exact
			      << "  |err|=" << std::abs(e_num - e_exact) << "\n";
	}
	std::cout << "Wrote files: hydrogen_energies.txt, hydrogen_convergence_dr.txt, hydrogen_convergence_rmax.txt, hydrogen_wavefunctions.txt\n";

	return 0;
}
