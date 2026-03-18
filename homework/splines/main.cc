#include <cassert>
#include <cmath>
#include <iostream>
#include "vector.h"

using pp::vector;

int binsearch(const vector& x, double z)
{
	assert(x.size() >= 2);
	assert(z >= x[0] && z <= x[x.size() - 1]);
	int i = 0, j = x.size() - 1;
	while (j - i > 1) {
		int mid = (i + j) / 2;
		if (z > x[mid])
			i = mid;
		else
			j = mid;
	}
	return i;
}

double linterp(const vector& x, const vector& y, double z)
{
	assert(x.size() == y.size());
	int i = binsearch(x, z);
	double dx = x[i + 1] - x[i];
	assert(dx > 0);
	double dy = y[i + 1] - y[i];
	return y[i] + dy / dx * (z - x[i]);
}

double linterpInteg(const vector& x, const vector& y, double z)
{
	assert(x.size() == y.size());
	assert(x.size() >= 2);
	assert(z >= x[0] && z <= x[x.size() - 1]);

	int i = binsearch(x, z);
	double sum = 0;

	for (int k = 0; k < i; k++) {
		double dx = x[k + 1] - x[k];
		assert(dx > 0);
		double slope = (y[k + 1] - y[k]) / dx;
		sum += y[k] * dx + 0.5 * slope * dx * dx;
	}

	{
		double dx = x[i + 1] - x[i];
		assert(dx > 0);
		double slope = (y[i + 1] - y[i]) / dx;
		double h = z - x[i];
		sum += y[i] * h + 0.5 * slope * h * h;
	}

	return sum;
}

int main()
{
	vector x, y;
	for (double xi = 0; xi <= 9; xi += 0.5) {
		x.push_back(xi);
		y.push_back(std::cos(xi));
	}

	std::cout << "# z  spline(z)  cos(z)  integ_spline(z)  sin(z)\n";
	for (double z = x[0]; z <= x[x.size() - 1]; z += 0.02) {
		double s = linterp(x, y, z);
		double S = linterpInteg(x, y, z);
		std::cout << z << " " << s << " " << std::cos(z) << " " << S << " " << std::sin(z) << "\n";
	}

	std::cout << "\n\n# knots: x_i  y_i\n";
	for (int i = 0; i < x.size(); i++)
		std::cout << x[i] << " " << y[i] << "\n";

	return 0;
}
