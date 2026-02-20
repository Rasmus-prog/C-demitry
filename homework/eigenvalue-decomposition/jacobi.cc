#include "jacobi.h"
#include <cmath>

void timesJ(Matrix& a, int p, int q, double theta) {
	const double c = std::cos(theta);
	const double s = std::sin(theta);
	for (int i = 0; i < a.rows; ++i) {
		const double aip = a(i, p);
		const double aiq = a(i, q);
		a(i, p) = c * aip - s * aiq;
		a(i, q) = s * aip + c * aiq;
	}
}

void Jtimes(Matrix& a, int p, int q, double theta) {
	const double c = std::cos(theta);
	const double s = std::sin(theta);
	for (int j = 0; j < a.cols; ++j) {
		const double apj = a(p, j);
		const double aqj = a(q, j);
		a(p, j) = c * apj + s * aqj;
		a(q, j) = -s * apj + c * aqj;
	}
}

std::pair<Vector, Matrix> jacobi_cyclic(const Matrix& m) {
	Matrix a = m;
	Matrix v = Matrix::identity(m.rows);
	const int n = m.rows;
	bool changed;

	do {
		changed = false;
		for (int p = 0; p < n - 1; ++p) {
			for (int q = p + 1; q < n; ++q) {
				const double apq = a(p, q);
				const double app = a(p, p);
				const double aqq = a(q, q);
				const double theta = 0.5 * std::atan2(2.0 * apq, aqq - app);
				const double c = std::cos(theta);
				const double s = std::sin(theta);
				const double new_app = c * c * app - 2.0 * s * c * apq + s * s * aqq;
				const double new_aqq = s * s * app + 2.0 * s * c * apq + c * c * aqq;

				if (new_app != app || new_aqq != aqq) {
					changed = true;
					timesJ(a, p, q, theta);
					Jtimes(a, p, q, -theta);
					timesJ(v, p, q, theta);
				}
			}
		}
	} while (changed);

	Vector w(static_cast<size_t>(n));
	for (int i = 0; i < n; ++i) w[static_cast<size_t>(i)] = a(i, i);
	return {w, v};
}