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

// Helper functions for optimized version - only update upper triangle
static void timesJ_opt(Matrix& a, int p, int q, double theta) {
	const double c = std::cos(theta);
	const double s = std::sin(theta);
	for (int i = 0; i < a.rows; ++i) {
		if (i == p || i == q) continue;
		const double aip = (i < p) ? a(i, p) : a(p, i);
		const double aiq = (i < q) ? a(i, q) : a(q, i);
		const double new_ip = c * aip - s * aiq;
		const double new_iq = s * aip + c * aiq;
		if (i < p) a(i, p) = new_ip; else a(p, i) = new_ip;
		if (i < q) a(i, q) = new_iq; else a(q, i) = new_iq;
	}
	// Update diagonal elements
	const double app = a(p, p);
	const double aqq = a(q, q);
	const double apq = (p < q) ? a(p, q) : a(q, p);
	a(p, p) = c * c * app - 2.0 * s * c * apq + s * s * aqq;
	a(q, q) = s * s * app + 2.0 * s * c * apq + c * c * aqq;
	a(p, q) = (c * c - s * s) * apq + s * c * (app - aqq);
	if (p > q) std::swap(a(p, q), a(q, p)); // ensure upper triangle stays upper
}

static void Jtimes_opt(Matrix& a, int p, int q, double theta) {
	const double c = std::cos(theta);
	const double s = std::sin(theta);
	for (int j = 0; j < a.cols; ++j) {
		if (j == p || j == q) continue;
		const double apj = (p < j) ? a(p, j) : a(j, p);
		const double aqj = (q < j) ? a(q, j) : a(j, q);
		const double new_pj = c * apj + s * aqj;
		const double new_qj = -s * apj + c * aqj;
		if (p < j) a(p, j) = new_pj; else a(j, p) = new_pj;
		if (q < j) a(q, j) = new_qj; else a(j, q) = new_qj;
	}
	// Diagonal already updated in timesJ_opt
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

std::pair<Vector, Matrix> jacobi_cyclic_optimized(const Matrix& m) {
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
					timesJ_opt(a, p, q, theta);
					Jtimes_opt(a, p, q, -theta);
					timesJ(v, p, q, theta);  // V is not symmetric, use standard version
				}
			}
		}
	} while (changed);

	Vector w(static_cast<size_t>(n));
	for (int i = 0; i < n; ++i) w[static_cast<size_t>(i)] = a(i, i);
	return {w, v};
}