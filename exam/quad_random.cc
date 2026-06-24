#include "quad_random.h"

#include <cmath>
#include <limits>
#include <random>
// Fixed seed for testing. In production, use std::random_device{}().
// static thread_local std::mt19937_64 rng{42};

std::random_device rd;
static thread_local std::mt19937_64 rng{rd()};
// ---------------------------------------------------------------------------
// Core recursive function.
//
// The key challenge for log(x)/sqrt(x) near x=0 even with CC transform:
// the transformed integrand g(t) = f(cos(pi*t)) * sin(pi*t) * (b-a)/2*pi
// still has large variance near t~1 (where cos(pi*t)~-1 -> x~0).
// We address this with:
//   1. Open interval (eps_t away from endpoints) in integrate_random_cc.
//   2. A hard cap at depth 52 (~double precision limit).
//   3. At the cap, we return the best estimate with a large error so the
//      caller's error accumulation stays finite.
// ---------------------------------------------------------------------------
static std::pair<double,double> recurse(
    const Function& f,
    double a, double b,
    double acc, double eps,
    int N,
    double inherited_mean,
    int depth)
{
    const double width = b - a;

    // Degenerate or depth cap: return best estimate, signal large error.
    if (width <= 0.0
        || !std::isfinite(width)
        || depth > 52)
    {
        // Don't return inf/nan — return a small contribution with large error.
        const double Q = (std::isnan(inherited_mean) ? 0.0 : inherited_mean) * width;
        const double safe_Q = std::isfinite(Q) ? Q : 0.0;
        // Mark as "not converged" by returning tol as error (safe, finite).
        return {safe_Q, acc + eps * std::abs(safe_Q) + 1e-30};
    }

    std::uniform_real_distribution<double> uniform(a, b);

    double sum  = 0.0;
    double sum2 = 0.0;
    bool   has_inf = false;
    for (int i = 0; i < N; ++i) {
        const double v = f(uniform(rng));
        if (!std::isfinite(v)) { has_inf = true; continue; }
        sum  += v;
        sum2 += v * v;
    }

    // If any sample was non-finite, recurse immediately (halve the interval
    // to escape the bad region) without checking tolerance.
    if (has_inf) {
        const double m       = (a + b) / 2.0;
        const double nextAcc = acc / std::sqrt(2.0);
        auto [Ql, el] = recurse(f, a, m, nextAcc, eps, N,
                                std::numeric_limits<double>::quiet_NaN(), depth + 1);
        auto [Qr, er] = recurse(f, m, b, nextAcc, eps, N,
                                std::numeric_limits<double>::quiet_NaN(), depth + 1);
        return {Ql + Qr, std::sqrt(el*el + er*er)};
    }

    const double mean_new = sum / N;
    const double Q        = mean_new * width;

    double err;
    if (std::isnan(inherited_mean)) {
        const double var = std::max(0.0, sum2 / N - mean_new * mean_new);
        err = 3.0 * std::sqrt(var / N) * width;
        if (err == 0.0) return {Q, 0.0};
    } else {
        err = std::abs(mean_new - inherited_mean) * width;
    }

    const double tol = acc + eps * std::abs(Q);
    if (err < tol) return {Q, err};

    const double m       = (a + b) / 2.0;
    const double nextAcc = acc / std::sqrt(2.0);

    auto [Ql, el] = recurse(f, a, m, nextAcc, eps, N, mean_new, depth + 1);
    auto [Qr, er] = recurse(f, m, b, nextAcc, eps, N, mean_new, depth + 1);

    return {Ql + Qr, std::sqrt(el*el + er*er)};
}

std::pair<double,double> integrate_random(
    const Function& f,
    double a, double b,
    double acc, double eps,
    int N,
    double, double)
{
    return recurse(f, a, b, acc, eps, N,
                   std::numeric_limits<double>::quiet_NaN(), 0);
}

// Clenshaw-Curtis wrapper.
// Open the t-interval by eps_t on each side to avoid the Jacobian zero
// (sin(pi*t)=0 at t=0,1) hitting a singularity of f at an endpoint.
std::pair<double,double> integrate_random_cc(
    const Function& f,
    double a, double b,
    double acc, double eps,
    int N)
{
    const double mid  = (a + b) / 2.0;
    const double half = (b - a) / 2.0;
    constexpr double eps_t = 1e-9;

    auto g = [&](double t) -> double {
        const double x    = mid + half * std::cos(M_PI * t);
        const double dxdt = half * M_PI * std::sin(M_PI * t);
        const double v    = f(x) * dxdt;
        return std::isfinite(v) ? v : 0.0;  // treat isolated singularity as 0
    };

    return recurse(g, eps_t, 1.0 - eps_t, acc, eps, N,
                   std::numeric_limits<double>::quiet_NaN(), 0);
}
