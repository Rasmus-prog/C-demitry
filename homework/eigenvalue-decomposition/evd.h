#pragma once
#include "matrix.h"
#include "vector.h"

namespace pp {

struct EVD {
	pp::vector w;
	pp::matrix V;

	static void timesJ(pp::matrix& A, int p, int q, double theta);
	static void Jtimes(pp::matrix& A, int p, int q, double theta);

	EVD(pp::matrix A);
};

}
