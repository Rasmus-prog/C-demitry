#pragma once

#include "vector.h"

using pp::vector;

struct cspline {
	const int n;
	vector x, y, b, c, d, I;

	cspline(const vector& xs, const vector& ys);
	int binsearch(double z) const;
	double eval(double z) const;
	double deriv(double z) const;
	double integ(double z) const;
};
