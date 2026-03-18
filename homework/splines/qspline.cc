#include <cassert>
#include "qspline.h"

qspline::qspline(const vector& xs, const vector& ys)
	: n(xs.size()), x(xs), y(ys), b(n - 1), c(n - 1), I(n)
{
	assert(n >= 2);
	assert(y.size() == n);

	vector dx(n - 1), p(n - 1);
	for (int i = 0; i < n - 1; i++) {
		dx[i] = x[i + 1] - x[i];
		assert(dx[i] > 0);
		p[i] = (y[i + 1] - y[i]) / dx[i];
	}

	c[0] = 0;
	for (int i = 0; i < n - 2; i++)
		c[i + 1] = (p[i + 1] - p[i] - c[i] * dx[i]) / dx[i + 1];

	c[n - 2] /= 2;
	for (int i = n - 3; i >= 0; i--)
		c[i] = (p[i + 1] - p[i] - c[i + 1] * dx[i + 1]) / dx[i];

	for (int i = 0; i < n - 1; i++)
		b[i] = p[i] - c[i] * dx[i];

	I[0] = 0;
	for (int i = 0; i < n - 1; i++) {
		double h = dx[i];
		I[i + 1] = I[i] + y[i] * h + b[i] * h * h / 2 + c[i] * h * h * h / 3;
	}
}

int qspline::binsearch(double z) const
{
	assert(z >= x[0] && z <= x[n - 1]);
	int i = 0, j = n - 1;
	while (j - i > 1) {
		int mid = (i + j) / 2;
		if (z > x[mid])
			i = mid;
		else
			j = mid;
	}
	return i;
}

double qspline::eval(double z) const
{
	int i = binsearch(z);
	double h = z - x[i];
	return y[i] + b[i] * h + c[i] * h * h;
}

double qspline::deriv(double z) const
{
	int i = binsearch(z);
	double h = z - x[i];
	return b[i] + 2 * c[i] * h;
}

double qspline::integ(double z) const
{
	int i = binsearch(z);
	double h = z - x[i];
	return I[i] + y[i] * h + b[i] * h * h / 2 + c[i] * h * h * h / 3;
}
