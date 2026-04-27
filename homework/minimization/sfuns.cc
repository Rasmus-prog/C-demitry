#include "sfuns.h"

#include <cmath>

namespace {

double norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double value : v) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        out[i] = a[i] + b[i];
    }
    return out;
}

std::vector<double> scale(const std::vector<double>& v, double alpha) {
    std::vector<double> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        out[i] = alpha * v[i];
    }
    return out;
}

std::vector<double> negate(const std::vector<double>& v) {
    std::vector<double> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        out[i] = -v[i];
    }
    return out;
}

std::vector<double> solve_linear_system_qr(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b
) {
    const std::size_t m = A.size();
    if (m == 0) {
        return {};
    }
    const std::size_t n = A.front().size();

    std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));

    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> v(n);
        for (std::size_t i = 0; i < n; ++i) {
            v[i] = A[i][j];
        }

        for (std::size_t i = 0; i < j; ++i) {
            double rij = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                rij += Q[k][i] * v[k];
            }
            R[i][j] = rij;
            for (std::size_t k = 0; k < n; ++k) {
                v[k] -= rij * Q[k][i];
            }
        }

        const double rjj = norm(v);
        R[j][j] = rjj;

        for (std::size_t k = 0; k < n; ++k) {
            Q[k][j] = v[k] / rjj;
        }
    }

    std::vector<double> y(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < n; ++k) {
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

} // namespace

std::vector<double> gradient(const Function& phi, const std::vector<double>& x) {
    const double phi_x = phi(x);
    std::vector<double> g(x.size());
    std::vector<double> x_work = x;

    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dxi = (1.0 + std::abs(x_work[i])) * std::pow(2.0, -26);
        x_work[i] += dxi;
        g[i] = (phi(x_work) - phi_x) / dxi;
        x_work[i] = x[i];
    }

    return g;
}

std::vector<std::vector<double>> hessian(const Function& phi, const std::vector<double>& x) {
    const std::vector<double> g_phi_x = gradient(phi, x);
    std::vector<std::vector<double>> H(x.size(), std::vector<double>(x.size(), 0.0));
    std::vector<double> x_work = x;

    for (std::size_t j = 0; j < x.size(); ++j) {
        const double dxj = (1.0 + std::abs(x_work[j])) * std::pow(2.0, -13);
        x_work[j] += dxj;
        const std::vector<double> g_shift = gradient(phi, x_work);
        for (std::size_t i = 0; i < x.size(); ++i) {
            H[i][j] = (g_shift[i] - g_phi_x[i]) / dxj;
        }
        x_work[j] = x[j];
    }

    return H;
}

NewtonResult newton_minimize(
    const Function& phi,
    const std::vector<double>& x0,
    double acc,
    int max_steps
) {
    std::vector<double> x = x0;
    NewtonResult result;

    for (int step = 0; step < max_steps; ++step) {
        const std::vector<double> g = gradient(phi, x);
        if (norm(g) < acc) {
            result.x = x;
            result.steps = step;
            result.converged = true;
            return result;
        }

        const std::vector<std::vector<double>> H = hessian(phi, x);
        std::vector<double> dx;
        try {
            dx = solve_linear_system_qr(H, negate(g));
        } catch (const std::exception&) {
            result.x = x;
            result.steps = step;
            result.converged = false;
            return result;
        }

        const double phi_x = phi(x);
        double lambda = 1.0;
        std::vector<double> x_trial = x;
        bool accepted = false;
        while (lambda >= 1.0 / 1024.0) {
            x_trial = add(x, scale(dx, lambda));
            if (phi(x_trial) < phi_x) {
                accepted = true;
                break;
            }
            lambda /= 2.0;
        }

        if (!accepted) {
            result.x = x;
            result.steps = step;
            result.converged = false;
            return result;
        }

        x = x_trial;
    }

    result.x = x;
    result.steps = max_steps;
    result.converged = false;
    return result;
}
