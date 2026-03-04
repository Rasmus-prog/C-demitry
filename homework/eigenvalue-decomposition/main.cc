#include "evd.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

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

int main() {
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
