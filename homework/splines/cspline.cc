#include <cassert>
#include "cspline.h"

cspline::cspline(const vector& xs, const vector& ys)
	: n(xs.size()), x(xs), y(ys), b(n - 1), c(n - 1), d(n - 1), I(n)
{
	assert(n >= 2);
	assert(y.size() == n);

	vector h(n - 1);
	for (int i = 0; i < n - 1; i++) {
		h[i] = x[i + 1] - x[i];
		assert(h[i] > 0);
	}

	vector alpha(n), l(n), mu(n), z(n), cfull(n);
	alpha[0] = 0;
	alpha[n - 1] = 0;
	for (int i = 1; i < n - 1; i++)
		alpha[i] = 3 * (y[i + 1] - y[i]) / h[i] - 3 * (y[i] - y[i - 1]) / h[i - 1];

	l[0] = 1;
	mu[0] = 0;
	z[0] = 0;
	for (int i = 1; i < n - 1; i++) {
		l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
		mu[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
	}
	l[n - 1] = 1;
	z[n - 1] = 0;
	cfull[n - 1] = 0;

	for (int i = n - 2; i >= 0; i--) {
		cfull[i] = z[i] - mu[i] * cfull[i + 1];
		b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (cfull[i + 1] + 2 * cfull[i]) / 3;
		c[i] = cfull[i];
		d[i] = (cfull[i + 1] - cfull[i]) / (3 * h[i]);
	}

	I[0] = 0;
	for (int i = 0; i < n - 1; i++) {
		double hi = h[i];
		I[i + 1] = I[i] + y[i] * hi + b[i] * hi * hi / 2 + c[i] * hi * hi * hi / 3 + d[i] * hi * hi * hi * hi / 4;
	}
}

int cspline::binsearch(double z) const
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

double cspline::eval(double z) const
{
	int i = binsearch(z);
	double h = z - x[i];
	return y[i] + b[i] * h + c[i] * h * h + d[i] * h * h * h;
}

double cspline::deriv(double z) const
{
	int i = binsearch(z);
	double h = z - x[i];
	return b[i] + 2 * c[i] * h + 3 * d[i] * h * h;
}

double cspline::integ(double z) const
{
	int i = binsearch(z);
	double h = z - x[i];
	return I[i] + y[i] * h + b[i] * h * h / 2 + c[i] * h * h * h / 3 + d[i] * h * h * h * h / 4;
}
