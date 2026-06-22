#pragma once

#include <functional>
#include <optional>
#include <utility>

using Function = std::function<double(double)>;

std::pair<double, double> integrate(const Function& f,
							   double a,
							   double b,
							   double acc = 1e-3,
							   double eps = 1e-3,
							   std::optional<double> f2 = std::nullopt,
							   std::optional<double> f3 = std::nullopt);

double errf_func(double x, double acc = 1e-3, double eps = 1e-3);

std::pair<double, double> integrate_clenshaw_curtis(const Function& f,
											   double a,
											   double b,
											   double acc = 1e-3,
											   double eps = 1e-3);
