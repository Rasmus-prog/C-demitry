#include "sfuns.h"
#include <cmath>
#include <vector>
#include <iostream>

namespace {

constexpr double kGradientStep = 1.0 / 67108864.0;
constexpr double kHessianStep = 1.0 / 8192.0;

std::vector<double> stepSizes(const std::vector<double>& x, double base) {
    std::vector<double> dx(x.size(), base);
    for (std::size_t i = 0; i < x.size(); ++i) {
        dx[i] = (1.0 + std::abs(x[i])) * base;
    }
    return dx;
}

std::vector<double> negate(const std::vector<double>& v) {
    std::vector<double> result(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        result[i] = -v[i];
    }
    return result;
}

std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> scale(const std::vector<double>& v, double alpha) {
    std::vector<double> result(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] * alpha;
    }
    return result;
}

std::vector<double> gradient(const Function& phi, const std::vector<double>& x) {
    const double phi_x = phi(x);
    const std::vector<double> dx = stepSizes(x, kGradientStep);
    std::vector<double> g(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::vector<double> xp = x;
        xp[i] += dx[i];
        g[i] = (phi(xp) - phi_x) / dx[i];
    }
    return g;
}

std::vector<std::vector<double>> hessian(const Function& phi, const std::vector<double>& x) {
    const std::vector<double> g0 = gradient(phi, x);
    const std::vector<double> dx = stepSizes(x, kHessianStep);
    std::vector<std::vector<double>> H(x.size(), std::vector<double>(x.size()));
    for (std::size_t j = 0; j < x.size(); ++j) {
        std::vector<double> xp = x;
        xp[j] += dx[j];
        const std::vector<double> gj = gradient(phi, xp);
        for (std::size_t i = 0; i < x.size(); ++i) {
            H[i][j] = (gj[i] - g0[i]) / dx[j];
        }
    }
    return H;
}

} // namespace


std::vector<double> solveLinearSystem(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    const std::size_t m = A.size();
    if (m == 0) {
        return {};
    }

    const std::size_t n = A.front().size();

    std::vector<std::vector<double>> Q(m, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));

    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> v(m);
        for (std::size_t i = 0; i < m; ++i) {
            v[i] = A[i][j];
        }

        for (std::size_t i = 0; i < j; ++i) {
            double rij = 0.0;
            for (std::size_t k = 0; k < m; ++k) {
                rij += Q[k][i] * v[k];
            }
            R[i][j] = rij;
            for (std::size_t k = 0; k < m; ++k) {
                v[k] -= rij * Q[k][i];
            }
        }

        const double rjj = norm(v);

        R[j][j] = rjj;
        for (std::size_t k = 0; k < m; ++k) {
            Q[k][j] = v[k] / rjj;
        }
    }

    std::vector<double> y(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < m; ++k) {
            y[i] += Q[k][i] * b[k];
        }
    }

    std::vector<double> x(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = y[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= R[i][j] * x[j];
        }
        x[i] = sum / R[i][i];
    }

    return x;
}


double norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double value : v) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

void printVector(const std::vector<double>& x) {
	std::cout << '(';
	for (std::size_t i = 0; i < x.size(); ++i) {
		std::cout << x[i];
		if (i + 1 != x.size()) {
			std::cout << ", ";
		}
	}
	std::cout << ')';
}

double rosenbrock(const std::vector<double>& x) {
	const double xx = x[0];
	const double yy = x[1];
	return (1.0 - xx) * (1.0 - xx) + 100.0 * (yy - xx * xx) * (yy - xx * xx);
}

double himmelblau(const std::vector<double>& x) {
	const double xx = x[0];
	const double yy = x[1];
	const double a = xx * xx + yy - 11.0;
	const double b = xx + yy * yy - 7.0;
	return a * a + b * b;
}

void reportSolution(
    const std::string& name,
    const Function& phi,
    const std::vector<double>& x0,
    double acc,
    double alpha_min,
    int max_iter
) {
    const NewtonResult result = newton(phi, x0, acc, alpha_min, max_iter);
	std::cout << name << "\n  start = ";
	printVector(x0);
    std::cout << "\n  minimum  = ";
    printVector(result.x);
    std::cout << "\n  f(minimum) = " << result.value;
    std::cout << "\n  steps = " << result.steps << "\n\n";
}


NewtonResult newton(
    const Function& phi,
    const std::vector<double>& x0,
    double acc,
    double alpha_min,
    int max_iter
) {
    std::vector<double> x = x0;
    int steps = 0;

    for (; steps < max_iter; ++steps) {
        const std::vector<double> g = gradient(phi, x);
        if (norm(g) < acc) {
            break;
        }

        const std::vector<std::vector<double>> H = hessian(phi, x);
        const std::vector<double> step = solveLinearSystem(H, negate(g));

        double alpha = 1.0;
        std::vector<double> candidate = x;
        double candidate_value = phi(x);
        while (alpha >= alpha_min) {
            candidate = add(x, scale(step, alpha));
            candidate_value = phi(candidate);
            if (candidate_value < phi(x)) {
                break;
            }
            alpha /= 2.0;
        }

        if (candidate_value >= phi(x)) {
            break;
        }

        x = candidate;
    }

    return {x, steps, phi(x)};
}