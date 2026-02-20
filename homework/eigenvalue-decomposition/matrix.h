#pragma once
#include <cstddef>
#include <vector>

struct Matrix {
	int rows = 0;
	int cols = 0;
	std::vector<double> data;

	Matrix() = default;
	Matrix(int r, int c);

	double& operator()(int i, int j);
	double operator()(int i, int j) const;

	static Matrix identity(int n);
};

using Vector = std::vector<double>;

Matrix transpose(const Matrix& a);

Matrix multiply(const Matrix& a, const Matrix& b);

Matrix diagonal_matrix(const Vector& w);

double max_abs(const Matrix& a);

double max_abs_diff(const Matrix& a, const Matrix& b);

double max_offdiag_abs(const Matrix& a);