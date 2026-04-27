#pragma once
#include <functional>
#include <vector>

using Function = std::function<double(const std::vector<double>&)>;

std::vector<double> newton(
	std::function<std::vector<double>(const std::vector<double>&)> f,
	const std::vector<double>& x0,
	double acc = 1e-2,
	double alpha_min = 1e-3,
	int max_iter = 100,
	const std::vector<double>& jacobian_dx = {}
);
