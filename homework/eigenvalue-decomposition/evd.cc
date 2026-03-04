#include "evd.h"
#include <cassert>
#include <cmath>

namespace pp {

void EVD::timesJ(pp::matrix& A, int p, int q, double theta) {
	double c = std::cos(theta), s = std::sin(theta);
	for (int i = 0; i < A.size1(); i++) {
		double aip = A[i, p], aiq = A[i, q];
		A[i, p] = c * aip - s * aiq;
		A[i, q] = s * aip + c * aiq;
	}
}

void EVD::Jtimes(pp::matrix& A, int p, int q, double theta) {
	double c = std::cos(theta), s = std::sin(theta);
	for (int j = 0; j < A.size2(); j++) {
		double apj = A[p, j], aqj = A[q, j];
		A[p, j] = c * apj + s * aqj;
		A[q, j] = -s * apj + c * aqj;
	}
}

EVD::EVD(pp::matrix A) : w(A.size1()), V(A.size1(), A.size1()) {
	assert(A.size1() == A.size2());
	int n = A.size1();
	V.setid();

	bool changed;
	do {
		changed = false;
		for (int p = 0; p < n - 1; p++) {
			for (int q = p + 1; q < n; q++) {
				double apq = A[p, q], app = A[p, p], aqq = A[q, q];
				double theta = 0.5 * std::atan2(2 * apq, aqq - app);
				double c = std::cos(theta), s = std::sin(theta);
				double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
				double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
				if (new_app != app || new_aqq != aqq) {
					changed = true;
					timesJ(A, p, q, theta);
					Jtimes(A, p, q, -theta);
					timesJ(V, p, q, theta);
				}
			}
		}
	} while (changed);

	for (int i = 0; i < n; i++) w[i] = A[i, i];
}

}
