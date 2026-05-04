#pragma once

#include <functional>
#include <vector>

struct activation_function {
	std::function<double(double)> f;
	std::function<double(double)> df;
};

activation_function gaussian();
activation_function gaussian_wavelet();
activation_function wavelet();

struct ann {
	int n;
	activation_function act;
	std::vector<double> p;

	ann(int n, activation_function act = gaussian_wavelet());

	double response(double x) const;
	double cost(const std::vector<double>& x, const std::vector<double>& y) const;
	std::vector<double> gradient(const std::vector<double>& x, const std::vector<double>& y) const;
	void train(const std::vector<double>& x, const std::vector<double>& y);
};

double tabulated_function(double x);
