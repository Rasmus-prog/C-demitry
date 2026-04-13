#include "sfuns.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>


// GAI start
namespace {

bool is_prime(int n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (int divisor = 3; divisor * divisor <= n; divisor += 2) {
        if (n % divisor == 0) return false;
    }
    return true;
}

std::vector<int> first_primes(int count) {
    std::vector<int> primes;
    primes.reserve(count);
    for (int candidate = 2; static_cast<int>(primes.size()) < count; ++candidate) {
        if (is_prime(candidate)) primes.push_back(candidate);
    }
    return primes;
}

double radical_inverse(int index, int base) {
    double value = 0.0;
    double factor = 1.0 / base;
    while (index > 0) {
        value += factor * (index % base);
        index /= base;
        factor /= base;
    }
    return value;
}

double map_to_bounds(const std::pair<double, double>& bound, double u) {
    return bound.first + u * (bound.second - bound.first);
}

} // namespace

// GAI end

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

std::pair<double, double> quasi_monte_carlo_integrate(
    const Function& f,
    const std::vector<std::pair<double, double>>& bounds,
    int num_samples) {
    int dim = static_cast<int>(bounds.size());
    if (num_samples <= 0 || dim == 0) return {0.0, 0.0};

    double volume = 1.0;
    for (const auto& bound : bounds) {
        volume *= (bound.second - bound.first);
    }

    std::vector<int> bases = first_primes(dim);
    std::vector<double> x(dim);

    int sequence_a_size = (num_samples + 1) / 2;
    int sequence_b_size = num_samples / 2;

    double sum_a = 0.0;
    double sum_b = 0.0;

    for (int i = 0; i < sequence_a_size; ++i) {
        int index = 2 * i + 1;
        for (int k = 0; k < dim; ++k) {
            x[k] = map_to_bounds(bounds[k], radical_inverse(index, bases[k]));
        }
        sum_a += f(x);
    }

    for (int i = 0; i < sequence_b_size; ++i) {
        int index = 2 * i + 2;
        for (int k = 0; k < dim; ++k) {
            x[k] = map_to_bounds(bounds[k], radical_inverse(index, bases[k]));
        }
        sum_b += f(x);
    }

    double mean_a = sum_a / sequence_a_size;
    double mean_b = sequence_b_size > 0 ? sum_b / sequence_b_size : mean_a;
    double estimate = (sequence_a_size * mean_a + sequence_b_size * mean_b) / num_samples;
    double error = 0.5 * std::abs(mean_a - mean_b) * volume;

    return {estimate * volume, error};
};