#include "evd.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

static double exact_energy(int n) {
	return -0.5 / (n * n);
}

static double assoc_laguerre_n_minus_1_alpha_1(int n, double x) {
	int kmax = n - 1;
	if (kmax == 0) return 1.0;
	double Lkm1 = 1.0;
	double Lk = 2.0 - x;
	if (kmax == 1) return Lk;
	for (int k = 1; k < kmax; k++) {
		double Lkp1 = ((2.0 * k + 2.0 - x) * Lk - (k + 1.0) * Lkm1) / (k + 1.0);
		Lkm1 = Lk;
		Lk = Lkp1;
	}
	return Lk;
}

static double exact_reduced_radial(int n, double r) {
	double x = 2.0 * r / n;
	double Ln = assoc_laguerre_n_minus_1_alpha_1(n, x);
	double pref = 2.0 * r / std::pow(n, 1.5);
	return pref * std::exp(-r / n) * Ln;
}

static double discrete_norm(const pp::vector& f, double dr) {
	double s = 0;
	for (int i = 0; i < f.size(); i++) s += f[i] * f[i];
	return std::sqrt(s * dr);
}

static void normalize(pp::vector& f, double dr) {
	double nrm = discrete_norm(f, dr);
	if (nrm > 0) f /= nrm;
}

static pp::matrix random_symmetric_matrix(int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(-2.0, 2.0);
	pp::matrix A(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			double x = dist(gen);
			A[i, j] = x;
			A[j, i] = x;
		}
	}
	return A;
}

static pp::matrix diag_from_vector(const pp::vector& d) {
	int n = d.size();
	pp::matrix D(n, n);
	for (int i = 0; i < n; i++) D[i, i] = d[i];
	return D;
}

static pp::matrix identity(int n) {
	pp::matrix I(n, n);
	I.setid();
	return I;
}

static int run_jacobi_selftest() {
	const int n = 5;
	const double tol = 1e-7;

	std::cout << std::setprecision(10);
	pp::matrix A = random_symmetric_matrix(n);
	A.print("A:");

	pp::EVD evd(A);
	pp::matrix V = evd.V;
	pp::matrix D = diag_from_vector(evd.w);
	pp::matrix VT = V.T();

	V.print("V:");
	D.print("D:");

	pp::matrix VTAV = VT * A * V;
	pp::matrix VDVT = V * D * VT;
	pp::matrix VTV = VT * V;
	pp::matrix VVT = V * VT;
	pp::matrix I = identity(n);

	VTAV.print("V^TAV:");
	std::cout << "V^TAV approx D? " << (pp::mat_approx(VTAV, D, tol) ? "yes" : "no") << "\n";
	std::cout << "\n";

	A.print("A:");
	VDVT.print("VDV^T:");
	std::cout << "VDV^T approx A? " << (pp::mat_approx(VDVT, A, tol) ? "yes" : "no") << "\n";
	std::cout << "\n";

	VTV.print("V^TV:");
	std::cout << "V^TV approx I? " << (pp::mat_approx(VTV, I, tol) ? "yes" : "no") << "\n";
	std::cout << "\n";

	VVT.print("VV^T:");
	std::cout << "VV^T approx I? " << (pp::mat_approx(VVT, I, tol) ? "yes" : "no") << "\n";

	return 0;
}

static bool build_hamiltonian(double rmax, double dr, pp::matrix& H, pp::vector& r) {
	if (!(rmax > 0 && dr > 0)) return false;
	int npoints = static_cast<int>(rmax / dr) - 1;
	if (npoints < 2) return false;

	r.resize(npoints);
	for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);

	H.resize(npoints, npoints);
	double kin = -0.5 / (dr * dr);
	for (int i = 0; i < npoints - 1; i++) {
		H[i, i] = -2 * kin;
		H[i, i + 1] = kin;
		H[i + 1, i] = kin;
	}
	H[npoints - 1, npoints - 1] = -2 * kin;
	for (int i = 0; i < npoints; i++) H[i, i] += -1.0 / r[i];
	return true;
}

static bool lowest_energy(double rmax, double dr, double& e0) {
	pp::matrix H;
	pp::vector r;
	if (!build_hamiltonian(rmax, dr, H, r)) return false;

	pp::EVD evd(H);
	e0 = evd.w[0];
	for (int i = 1; i < evd.w.size(); i++) {
		if (evd.w[i] < e0) e0 = evd.w[i];
	}
	return true;
}

static int run_convergence(double fixed_rmax, double fixed_dr) {
	std::ofstream drout("convergence_dr.dat");
	std::ofstream rmaxout("convergence_rmax.dat");
	if (!drout || !rmaxout) {
		std::cerr << "Failed to open convergence output files.\n";
		return 1;
	}

	drout << std::setprecision(10);
	rmaxout << std::setprecision(10);
	drout << "# fixed_rmax=" << fixed_rmax << "\n";
	drout << "# dr epsilon0 exact abs_error\n";
	rmaxout << "# fixed_dr=" << fixed_dr << "\n";
	rmaxout << "# rmax epsilon0 exact abs_error\n";

	const double exact0 = exact_energy(1);

	std::cout << "Convergence sweep 1: epsilon0 vs dr (fixed rmax=" << fixed_rmax << ")\n";
	for (double sweep_dr : std::vector<double>{0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.10, 0.08, 0.06, 0.05}) {
		double e0 = std::numeric_limits<double>::quiet_NaN();
		if (lowest_energy(fixed_rmax, sweep_dr, e0)) {
			drout << sweep_dr << " " << e0 << " " << exact0 << " " << std::abs(e0 - exact0) << "\n";
			std::cout << "  dr=" << std::setw(6) << sweep_dr << "  epsilon0=" << std::setw(14) << e0
			          << "  abs.error=" << std::abs(e0 - exact0) << "\n";
		}
	}

	std::cout << "\nConvergence sweep 2: epsilon0 vs rmax (fixed dr=" << fixed_dr << ")\n";
	for (double sweep_rmax : std::vector<double>{4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0}) {
		double e0 = std::numeric_limits<double>::quiet_NaN();
		if (lowest_energy(sweep_rmax, fixed_dr, e0)) {
			rmaxout << sweep_rmax << " " << e0 << " " << exact0 << " " << std::abs(e0 - exact0) << "\n";
			std::cout << "  rmax=" << std::setw(6) << sweep_rmax << "  epsilon0=" << std::setw(14) << e0
			          << "  abs.error=" << std::abs(e0 - exact0) << "\n";
		}
	}

	std::cout << "\nSaved convergence data to convergence_dr.dat and convergence_rmax.dat\n";
	return 0;
}

static int run_wavefunctions(double rmax, double dr, int requested_states) {
	pp::matrix H;
	pp::vector r;
	if (!build_hamiltonian(rmax, dr, H, r)) {
		std::cerr << "Invalid grid for wavefunctions task: ensure rmax>0, dr>0 and at least 2 grid points.\n";
		return 1;
	}

	pp::EVD evd(H);
	int npoints = r.size();
	std::vector<int> order(npoints);
	for (int i = 0; i < npoints; i++) order[i] = i;
	std::sort(order.begin(), order.end(), [&](int a, int b) { return evd.w[a] < evd.w[b]; });

	int nstates = std::min(std::max(requested_states, 1), npoints);
	double c = 1.0 / std::sqrt(dr);

	std::ofstream summary("wavefunctions_summary.txt");
	if (!summary) {
		std::cerr << "Failed to open wavefunctions_summary.txt\n";
		return 1;
	}
	summary << std::setprecision(10);
	summary << "# rmax=" << rmax << " dr=" << dr << " Const=" << c << "\n";
	summary << "# n epsilon_numeric epsilon_exact abs_error overlap\n";

	std::cout << "Wave-functions task: compare lowest " << nstates
	          << " eigen-functions with analytical results\n";
	std::cout << "Using Const = 1/sqrt(dr) = " << c << "\n";

	for (int s = 0; s < nstates; s++) {
		int n = s + 1;
		int col = order[s];
		double en = evd.w[col];
		double ex = exact_energy(n);

		pp::vector f_num = evd.V[col];
		f_num *= c;

		pp::vector f_ex(r.size());
		for (int i = 0; i < r.size(); i++) f_ex[i] = exact_reduced_radial(n, r[i]);
		normalize(f_ex, dr);

		double overlap = 0;
		for (int i = 0; i < r.size(); i++) overlap += f_num[i] * f_ex[i];
		overlap *= dr;
		if (overlap < 0) {
			for (int i = 0; i < r.size(); i++) f_num[i] = -f_num[i];
			overlap = -overlap;
		}

		std::string filename = "wave_n" + std::to_string(n) + ".dat";
		std::ofstream fout(filename);
		if (!fout) {
			std::cerr << "Failed to open output file: " << filename << "\n";
			return 1;
		}
		fout << std::setprecision(10);
		fout << "# n=" << n << " rmax=" << rmax << " dr=" << dr << " Const=" << c << "\n";
		fout << "# r f_num f_exact\n";
		for (int i = 0; i < r.size(); i++) fout << r[i] << " " << f_num[i] << " " << f_ex[i] << "\n";

		summary << n << " " << en << " " << ex << " " << std::abs(en - ex) << " " << overlap << "\n";
		std::cout << "  n=" << n << "  overlap=" << overlap << "  saved " << filename << "\n";
	}

	std::cout << "Saved summary to wavefunctions_summary.txt\n";
	return 0;
}

int main(int argc, char** argv) {
	double rmax = 10.0;
	double dr = 0.3;
	bool jacobi_selftest = false;
	bool do_convergence = false;
	bool do_wavefunctions = false;
	double conv_fixed_rmax = 10.0;
	double conv_fixed_dr = 0.1;
	int wf_states = 3;

	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-rmax" && i + 1 < argc) {
			rmax = std::stod(argv[++i]);
		} else if (arg == "-dr" && i + 1 < argc) {
			dr = std::stod(argv[++i]);
		} else if (arg == "--jacobi-test") {
			jacobi_selftest = true;
		} else if (arg == "--convergence") {
			do_convergence = true;
		} else if (arg == "--wavefunctions") {
			do_wavefunctions = true;
		} else if (arg == "--wf-states" && i + 1 < argc) {
			wf_states = std::stoi(argv[++i]);
		} else if (arg == "--conv-fixed-rmax" && i + 1 < argc) {
			conv_fixed_rmax = std::stod(argv[++i]);
		} else if (arg == "--conv-fixed-dr" && i + 1 < argc) {
			conv_fixed_dr = std::stod(argv[++i]);
		} else {
			std::cerr << "Unknown or incomplete option: " << arg << "\n";
			std::cerr << "Usage: ./main -rmax <value> -dr <value> [--jacobi-test] [--convergence]"
			          << " [--conv-fixed-rmax <value>] [--conv-fixed-dr <value>]"
			          << " [--wavefunctions] [--wf-states <int>]\n";
			return 1;
		}
	}

	if (jacobi_selftest) return run_jacobi_selftest();
	if (do_convergence) return run_convergence(conv_fixed_rmax, conv_fixed_dr);
	if (do_wavefunctions) return run_wavefunctions(rmax, dr, wf_states);

	std::cout << "No explicit mode selected; running hydrogen solver with -rmax and -dr.\n";
	return run_wavefunctions(rmax, dr, wf_states);
}
