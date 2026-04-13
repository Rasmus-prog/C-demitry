#include "sfuns.h"

#include <cmath>
#include <random>
#include <vector>


std::pair<double, double> monte_carlo_integrate(
    const Function& f,
    const std::vector<std::pair<double, double>>& bounds,
    int num_samples) {
    int dim = bounds.size();
    double volume = 1.0;
    for (const auto& b : bounds) {
        volume *= (b.second - b.first);
    }
    std::vector<double> x(dim);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double sum = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        for (int k = 0; k < dim; ++k) {
            x[k] = bounds[k].first + dis(gen) * (bounds[k].second - bounds[k].first);
        }
        double fx = f(x);
        sum += fx;
        sum2 += fx * fx;
    }
    double mean = sum / num_samples;
    double variance = sum2 / num_samples - mean * mean;
    if (variance < 0.0) variance = 0.0;
    double sigma = std::sqrt(variance);

    return {mean * volume, sigma * volume / std::sqrt(num_samples)};
}