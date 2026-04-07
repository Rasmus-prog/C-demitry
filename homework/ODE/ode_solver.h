#pragma once

#include <functional>
#include <tuple>
#include <vector>
#include "vector.h"

namespace ode {

using Vec = pp::vector;
using OdeFunction = std::function<Vec(double, const Vec&)>;

std::tuple<Vec, Vec> rkstep12(const OdeFunction& f, double x, const Vec& y, double h);

std::tuple<std::vector<double>, std::vector<Vec>> driver(
	const OdeFunction& f,
	double a,
	double b,
	const Vec& yinit,
	double h = 0.125,
	double acc = 0.01,
	double eps = 0.01
);

} // namespace ode
