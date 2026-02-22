#pragma once
#include <utility>
#include "matrix.h"

void timesJ(Matrix& a, int p, int q, double theta);

void Jtimes(Matrix& a, int p, int q, double theta);

std::pair<Vector, Matrix> jacobi_cyclic(const Matrix& m);

std::pair<Vector, Matrix> jacobi_cyclic_optimized(const Matrix& m);