#pragma once
#include <functional>
#include <vector>
#include <string>


using Function = std::function<double(const std::vector<double>&)>;

std::vector<double> newton(
	std::function<std::vector<double>(const std::vector<double>&)> f,
	const std::vector<double>& x0,
	double acc = 1e-2,
	double alpha_min = 1e-3,
	int max_iter = 100,
	const std::vector<double>& jacobian_dx = {}
);

// norm func
double norm(const std::vector<double>& v);

// print vector func
void printVector(const std::vector<double>& x);

// rosenbrock func
double rosenbrock(const std::vector<double>& x);

// himmelblau func
double himmelblau(const std::vector<double>& x);

// report solution func
void reportSolution(
	const std::string& name,
	const std::function<std::vector<double>(const std::vector<double>&)>& f,
	const std::vector<double>& x0,
	const std::vector<double>& dx = {}
);

// negate vector func
std::vector<double> negate(const std::vector<double>& v);

// add vector func
std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);

// scale vector func
std::vector<double> scale(const std::vector<double>& v, double alpha);

// default jacobian dx func
std::vector<double> defaultJacobianDx(const std::vector<double>& x);

// compute jacobian func
std::vector<std::vector<double>> computeJacobian(
	std::function<std::vector<double>(const std::vector<double>&)> f,
	const std::vector<double>& x,
	const std::vector<double>& fx,
	const std::vector<double>& dx
);

// solve linear system func
std::vector<double> solveLinearSystem(
	const std::vector<std::vector<double>>& A,
	const std::vector<double>& b
);