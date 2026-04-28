#include "sfuns.h"
#include <cmath>
#include <limits>
#include <vector>
#include <iostream>


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

std::vector<double> defaultJacobianDx(const std::vector<double>& x) {
    constexpr double eps_sqrt = 1.0 / 67108864.0; // 2^-26
    std::vector<double> dx(x.size(), eps_sqrt);
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double scale = std::max(std::abs(x[i]), 1.0);
        dx[i] = scale * eps_sqrt;
    }
    return dx;
}

std::vector<std::vector<double>> computeJacobian(
    std::function<std::vector<double>(const std::vector<double>&)> f,
    const std::vector<double>& x,
    const std::vector<double>& fx,
    const std::vector<double>& dx
) {
    const std::size_t n = x.size();
    const std::size_t m = fx.size();
    std::vector<std::vector<double>> J(m, std::vector<double>(n));

    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> x_plus = x;
        const double h = dx[j];
        x_plus[j] += h;
        const std::vector<double> fx_plus = f(x_plus);

        for (std::size_t i = 0; i < m; ++i) {
            J[i][j] = (fx_plus[i] - fx[i]) / h;
        }
    }

    return J;
}

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
    const std::function<std::vector<double>(const std::vector<double>&)>& f,
    const std::vector<double>& x0,
    const std::vector<double>& dx
) {
	const std::vector<double> root = newton(f, x0, 1e-10, 1e-3, 200, dx);
	std::cout << name << "\n  start = ";
	printVector(x0);
	std::cout << "\n  root  = ";
	printVector(root);
	std::cout << "\n  ||f(root)|| = " << norm(f(root)) << "\n\n";
}


std::vector<double> newton(
    std::function<std::vector<double>(const std::vector<double>&)> f,
    const std::vector<double>& x0,
    double acc,
    double alpha_min,
    int max_iter,
    const std::vector<double>& jacobian_dx
) {
    std::vector<double> x = x0;
    const bool has_user_dx = !jacobian_dx.empty();
    std::vector<double> user_dx = jacobian_dx;
    for (double& d : user_dx) {
        d = std::abs(d);
        if (d == 0.0) {
            d = 1.0 / 67108864.0;
        }
    }

    std::vector<double> fx = f(x);

    for (int iter = 0; iter < max_iter; ++iter) {
        if (norm(fx) < acc) {
            break;
        }

        std::vector<double> dx = has_user_dx ? user_dx : defaultJacobianDx(x);
        for (double& d : dx) {
            d = std::abs(d);
            if (d == 0.0) {
                d = 1.0 / 67108864.0;
            }
        }

        const std::vector<std::vector<double>> J = computeJacobian(f, x, fx, dx);
        const std::vector<double> Dx = solveLinearSystem(J, negate(fx));

        double alpha = 1.0;
        std::vector<double> z = x;
        std::vector<double> fz = fx;
        while (true) {
            z = add(x, scale(Dx, alpha));
            fz = f(z);
            if (norm(fz) < norm(fx)) {
                break;
            }
            if (alpha < alpha_min) {
                break;
            }
            alpha /= 2.0;
        }

        const std::vector<double> step = scale(Dx, alpha);
        if (norm(step) <= norm(dx)) {
            x = z;
            fx = fz;
            break;
        }

        x = z;
        fx = fz;
    }

    return x;
}