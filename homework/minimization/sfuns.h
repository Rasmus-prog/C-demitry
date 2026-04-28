#pragma once
#include <functional>
#include <vector>
#include <string>


using Function = std::function<double(const std::vector<double>&)>;

struct NewtonResult {
	std::vector<double> x;
	int steps;
	double value;
};

NewtonResult newton(
	const Function& phi,
	const std::vector<double>& x0,
	double acc = 1e-3,
	double alpha_min = 1.0 / 1024.0,
	int max_iter = 1000
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
	const Function& phi,
	const std::vector<double>& x0,
	double acc = 1e-3,
	double alpha_min = 1.0 / 1024.0,
	int max_iter = 1000
);

// solve linear system func
std::vector<double> solveLinearSystem(
	const std::vector<std::vector<double>>& A,
	const std::vector<double>& b
);