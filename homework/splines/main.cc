#include <cassert>
#include <cmath>
#include <iostream>
#include "qspline.h"
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

double linterpDeriv(const vector& x, const vector& y, double z)
{
	assert(x.size() == y.size());
	int i = binsearch(x, z);
	double dx = x[i + 1] - x[i];
	assert(dx > 0);
	return (y[i + 1] - y[i]) / dx;
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

void run_debug_cases()
{
	vector x;
	for (int i = 1; i <= 5; i++)
		x.push_back(i);

	auto print_bc = [](const qspline& s, const char* name) {
		std::cerr << "# Debug case: " << name << "\n";
		std::cerr << "# i  b_i  c_i\n";
		for (int i = 0; i < s.b.size(); i++)
			std::cerr << i << " " << s.b[i] << " " << s.c[i] << "\n";
		std::cerr << "\n";
	};

	{
		vector y;
		for (int i = 1; i <= 5; i++)
			y.push_back(1.0);
		qspline s(x, y);
		print_bc(s, "y=1");
		for (int i = 0; i < s.b.size(); i++) {
			assert(pp::approx(s.b[i], 0.0));
			assert(pp::approx(s.c[i], 0.0));
		}
	}

	{
		vector y;
		for (int i = 1; i <= 5; i++)
			y.push_back(i);
		qspline s(x, y);
		print_bc(s, "y=x");
		for (int i = 0; i < s.b.size(); i++) {
			assert(pp::approx(s.b[i], 1.0));
			assert(pp::approx(s.c[i], 0.0));
		}
	}

	{
		vector y;
		for (int i = 1; i <= 5; i++)
			y.push_back(i * i);
		qspline s(x, y);
		print_bc(s, "y=x^2");
		for (int i = 0; i < s.b.size(); i++) {
			double expected_b = 2.0 * x[i];
			assert(pp::approx(s.b[i], expected_b));
			assert(pp::approx(s.c[i], 1.0));
		}
	}
}

int main()
{
	run_debug_cases();

	vector x, y;
	for (double xi = 0; xi <= 9; xi += 0.5) {
		x.push_back(xi);
		y.push_back(std::cos(xi));
	}
	qspline s(x, y);

	std::cout << "# z  linear(z)  quadratic(z)  cos(z)  linear'(z)  quadratic'(z)  -sin(z)  integ_linear(z)  integ_quadratic(z)  sin(z)\n";
	for (double z = x[0]; z <= x[x.size() - 1]; z += 0.02) {
		double sl = linterp(x, y, z);
		double sq = s.eval(z);
		double dl = linterpDeriv(x, y, z);
		double dq = s.deriv(z);
		double Il = linterpInteg(x, y, z);
		double Iq = s.integ(z);
		std::cout << z << " " << sl << " " << sq << " " << std::cos(z) << " " << dl << " " << dq << " " << -std::sin(z) << " " << Il << " " << Iq << " " << std::sin(z) << "\n";
	}

	std::cout << "\n\n# knots: x_i  y_i\n";
	for (int i = 0; i < x.size(); i++)
		std::cout << x[i] << " " << y[i] << "\n";

	return 0;
}
