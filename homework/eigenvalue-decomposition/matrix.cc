#include "matrix.h"
#include <algorithm>
#include <cmath>

Matrix::Matrix(int r, int c) : rows(r), cols(c), data(static_cast<size_t>(r * c), 0.0) {}

double& Matrix::operator()(int i, int j) { return data[static_cast<size_t>(i * cols + j)]; }

double Matrix::operator()(int i, int j) const { return data[static_cast<size_t>(i * cols + j)]; }

Matrix Matrix::identity(int n) {
	Matrix id(n, n);
	for (int i = 0; i < n; ++i) id(i, i) = 1.0;
	return id;
}

Matrix transpose(const Matrix& a) {
	Matrix t(a.cols, a.rows);
	for (int i = 0; i < a.rows; ++i)
		for (int j = 0; j < a.cols; ++j) t(j, i) = a(i, j);
	return t;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
	Matrix c(a.rows, b.cols);
	for (int i = 0; i < a.rows; ++i) {
		for (int k = 0; k < a.cols; ++k) {
			const double aik = a(i, k);
			for (int j = 0; j < b.cols; ++j) c(i, j) += aik * b(k, j);
		}
	}
	return c;
}

Matrix diagonal_matrix(const Vector& w) {
	Matrix d(static_cast<int>(w.size()), static_cast<int>(w.size()));
	for (int i = 0; i < static_cast<int>(w.size()); ++i) d(i, i) = w[static_cast<size_t>(i)];
	return d;
}

double max_abs(const Matrix& a) {
	double m = 0.0;
	for (double v : a.data) m = std::max(m, std::abs(v));
	return m;
}

double max_abs_diff(const Matrix& a, const Matrix& b) {
	double m = 0.0;
	for (size_t i = 0; i < a.data.size(); ++i) m = std::max(m, std::abs(a.data[i] - b.data[i]));
	return m;
}

double max_offdiag_abs(const Matrix& a) {
	double m = 0.0;
	for (int i = 0; i < a.rows; ++i)
		for (int j = 0; j < a.cols; ++j)
			if (i != j) m = std::max(m, std::abs(a(i, j)));
	return m;
}