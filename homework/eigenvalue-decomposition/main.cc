#include "evd.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

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

int main(int argc, char** argv) {
	double rmax = 10.0;
	double dr = 0.3;
	bool jacobi_selftest = false;

	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-rmax" && i + 1 < argc) {
			rmax = std::stod(argv[++i]);
		} else if (arg == "-dr" && i + 1 < argc) {
			dr = std::stod(argv[++i]);
		} else if (arg == "--jacobi-test") {
			jacobi_selftest = true;
		} else {
			std::cerr << "Unknown or incomplete option: " << arg << "\n";
			std::cerr << "Usage: ./main -rmax <value> -dr <value> [--jacobi-test]\n";
			return 1;
		}
	}

	if (jacobi_selftest) return run_jacobi_selftest();

	if (!(rmax > 0 && dr > 0)) {
		std::cerr << "Both rmax and dr must be positive.\n";
		return 1;
	}

	int npoints = static_cast<int>(rmax / dr) - 1;
	if (npoints < 2) {
		std::cerr << "Need at least 2 grid points. Increase rmax or decrease dr.\n";
		return 1;
	}

	pp::vector r(npoints);
	for (int i = 0; i < npoints; i++) r[i] = dr * (i + 1);

	pp::matrix H(npoints, npoints);
	double kin = -0.5 / (dr * dr);
	for (int i = 0; i < npoints - 1; i++) {
		H[i, i] = -2 * kin;
		H[i, i + 1] = 1 * kin;
		H[i + 1, i] = 1 * kin;
	}
	H[npoints - 1, npoints - 1] = -2 * kin;
	for (int i = 0; i < npoints; i++) H[i, i] += -1.0 / r[i];

	pp::EVD evd(H);

	std::vector<int> order(npoints);
	for (int i = 0; i < npoints; i++) order[i] = i;
	std::sort(order.begin(), order.end(), [&](int a, int b) { return evd.w[a] < evd.w[b]; });

	int nshow = std::min(5, npoints);
	std::cout << std::setprecision(10);
	std::cout << "Hydrogen s-wave on grid: rmax=" << rmax << ", dr=" << dr << ", npoints=" << npoints
	          << "\n\n";
	std::cout << "Lowest eigenvalues (Hartree):\n";
	std::cout << " state   numeric            exact              abs.error\n";
	for (int s = 0; s < nshow; s++) {
		int n = s + 1;
		double en = evd.w[order[s]];
		double ex = exact_energy(n);
		std::cout << "  " << std::setw(2) << n << "    " << std::setw(14) << en << "  " << std::setw(14) << ex
		          << "  " << std::setw(14) << std::abs(en - ex) << "\n";
	}

	int wf_states = std::min(3, nshow);
	std::cout << "\nRadial reduced wavefunctions f_n(r):\n";
	std::cout << "Each table compares normalized numerical and exact f_n for n=1.." << wf_states << "\n";
	for (int s = 0; s < wf_states; s++) {
		int n = s + 1;
		pp::vector f_num = evd.V[order[s]];
		normalize(f_num, dr);

		pp::vector f_ex(npoints);
		for (int i = 0; i < npoints; i++) f_ex[i] = exact_reduced_radial(n, r[i]);
		normalize(f_ex, dr);

		double overlap = 0;
		for (int i = 0; i < npoints; i++) overlap += f_num[i] * f_ex[i];
		overlap *= dr;
		if (overlap < 0) {
			for (int i = 0; i < npoints; i++) f_num[i] = -f_num[i];
			overlap = -overlap;
		}

		std::cout << "\nState n=" << n << ", overlap=<f_num|f_exact>=" << overlap << "\n";
		std::string filename = "wavefunction_n" + std::to_string(n) + ".dat";
		std::ofstream fout(filename);
		if (!fout) {
			std::cerr << "Failed to open output file: " << filename << "\n";
			return 1;
		}
		fout << std::setprecision(10);
		fout << "# n=" << n << " overlap=" << overlap << "\n";
		fout << "# r f_num f_exact\n";
		for (int i = 0; i < npoints; i++) {
			fout << r[i] << " " << f_num[i] << " " << f_ex[i] << "\n";
		}
		std::cout << "Saved table to " << filename << "\n";
	}

	return 0;
}
