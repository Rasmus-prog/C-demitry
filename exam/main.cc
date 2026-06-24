#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "sfuns.h"
#include "quad_random.h"

// Helpers

struct Result {
    double value;
    double err_est;
    double act_err;
    int    ncalls;
};

Result run_fixed(const Function& base, double a, double b, double acc = 1e-3, double eps = 1e-3, double exact = 0.0) {
    int n = 0;
    auto counted = [&](double x) { ++n; return base(x); };
    auto [v, e] = integrate(counted, a, b, acc, eps);
    return {v, e, std::abs(v - exact), n};
}

Result run_fixed_cc(const Function& base, double a, double b, double acc = 1e-3, double eps = 1e-3, double exact = 0.0) {
    int n = 0;
    auto counted = [&](double x) { ++n; return base(x); };
    auto [v, e] = integrate_clenshaw_curtis(counted, a, b, acc, eps);
    return {v, e, std::abs(v - exact), n};
}

Result run_random(const Function& base, double a, double b, int N,
                  double acc = 1e-3, double eps = 1e-3, double exact = 0.0) {
    int n = 0;
    auto counted = [&](double x) { ++n; return base(x); };
    auto [v, e] = integrate_random(counted, a, b, acc, eps, N);
    return {v, e, std::abs(v - exact), n};
}

Result run_random_cc(const Function& base, double a, double b, int N,
                     double acc = 1e-3, double eps = 1e-3, double exact = 0.0) {
    int n = 0;
    auto counted = [&](double x) { ++n; return base(x); };
    auto [v, e] = integrate_random_cc(counted, a, b, acc, eps, N);
    return {v, e, std::abs(v - exact), n};
}



// Write one comparison .data file for a single integral.
// Columns: method  ncalls  act_err  est_err  value
// Methods: fixed, fixed_cc, random_N2, random_N4, random_N8, random_N16,
//          random_cc_N2, random_cc_N4, random_cc_N8, random_cc_N16
void write_comparison(const std::string& filename,
                      const std::string& title,
                      const Function& f,
                      double a, double b, double exact,
                      double acc = 1e-3, double eps = 1e-3)
{
    std::ofstream out(filename);
    out << "# " << title << "  exact=" << exact << "\n";
    out << "# method              ncalls   act_err      est_err      value\n";
    out << std::scientific << std::setprecision(6);

    auto row = [&](const std::string& name, const Result& r) {
        out << std::left << std::setw(20) << name << "  "
            << std::setw(8) << r.ncalls << "  "
            << std::setw(12) << r.act_err << "  "
            << std::setw(12) << r.err_est << "  "
            << r.value << "\n";
        std::cout << std::left << std::setw(20) << name
                  << "  ncalls=" << std::setw(7) << r.ncalls
                  << "  act_err=" << r.act_err << "\n";
    };

    std::cout << "\n--- " << title << " ---\n";
    row("fixed",         run_fixed(f, a, b, acc, eps, exact));
    row("fixed_cc",      run_fixed_cc(f, a, b, acc, eps, exact));
    for (int N : {2, 4, 8, 16, 32, 64, 128, 256}) {
        row("random_N"    + std::to_string(N), run_random(f, a, b, N, acc, eps, exact));
        row("random_cc_N" + std::to_string(N), run_random_cc(f, a, b, N, acc, eps, exact));
    }
}

// Write N-sweep .data file for a single integral (random+CC, averaged over REPS runs).
// Columns: N  ncalls  act_err
void write_sweep(const std::string& filename,
                 const std::string& title,
                 const Function& f,
                 double a, double b, double exact,
                 double acc = 1e-4, double eps = 1e-4,
                 int REPS = 8)
{
    std::ofstream out(filename);
    out << "# N sweep — " << title << "  acc=eps=" << acc << "  averaged over " << REPS << " runs\n";
    out << "# N    ncalls   act_err\n";
    out << std::scientific << std::setprecision(6);

    std::cout << "\n--- N sweep: " << title << " ---\n";
    std::cout << std::left << std::setw(6) << "N"
              << std::setw(10) << "ncalls" << "act_err\n";

    for (int N : {2, 4, 8, 16, 32, 64, 128, 256}) {
        double total_calls = 0, total_err = 0;
        for (int r = 0; r < REPS; ++r) {
            auto res = run_random_cc(f, a, b, N, acc, eps, exact);
            total_calls += res.ncalls;
            total_err   += res.act_err;
        }
        const int    avg_calls = static_cast<int>(total_calls / REPS);
        const double avg_err   = total_err / REPS;
        out << N << "  " << avg_calls << "  " << avg_err << "\n";
        std::cout << std::left << std::setw(6) << N
                  << std::setw(10) << avg_calls << avg_err << "\n";
    }
}



// Print a comparison row
void print_row(const std::string& method, const Result& r, double exact) {
    const double actual_err = std::abs(r.value - exact);
    std::cout << std::left << std::setw(26) << method
              << " val=" << std::setw(10) << std::setprecision(6) << r.value
              << " est_err=" << std::setw(10) << r.err_est
              << " act_err=" << std::setw(10) << actual_err
              << " ncalls=" << r.ncalls << '\n';
}

int main() {

    // Intergrals
    struct Case { std::string name; Function f; double a, b, exact; };


    const std::vector<Case> cases = {
        {"sqrt(x)",        [](double x){ return std::sqrt(x); },            0, 1,  2.0/3.0},
        {"1_sqrt(x)",      [](double x){ return 1.0/std::sqrt(x); },        0, 1,  2.0},
        {"ln(x)_sqrt(x)",  [](double x){ return std::log(x)/std::sqrt(x);}, 0, 1, -4.0},
        {"sqrt(1-x^2)",    [](double x){ return std::sqrt(1-x*x); },        0, 1, M_PI/4.0},
    };

    // Pretty titles for plot labels
    const std::vector<std::string> titles = {
        "integral of sqrt(x) on [0,1] = 2/3",
        "integral of 1/sqrt(x) on [0,1] = 2  (singularity at x=0)",
        "integral of ln(x)/sqrt(x) on [0,1] = -4  (singularity at x=0)",
        "integral of sqrt(1-x^2) on [0,1] = pi/4",
    };


    // Write comparison files: data/<name>.data
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& c = cases[i];
        write_comparison("data/" + c.name + ".data", titles[i],
                         c.f, c.a, c.b, c.exact);
    }

    // Write N-sweep files: data/sweep_<name>.data
    // (only for the three integrals used in the sweep — skip circle)
    for (std::size_t i = 0; i < 3; ++i) {
        const auto& c = cases[i];
        write_sweep("data/sweep_" + c.name + ".data",
                    titles[i], c.f, c.a, c.b, c.exact);
    }

    std::cout << "\nData written to data/\n";

    // std::cout << std::fixed;

    // // 1. Smooth integral: ∫01 √x dx = 2/3
    // {
    //     std::cout << "\n=== ∫₀¹ √x dx  (exact = 2/3) ===\n";
    //     auto f = [](double x) { return std::sqrt(x); };
    //     const double exact = 2.0 / 3.0;

    //     print_row("fixed (open 4-pt)",    run_fixed(f, 0, 1),         exact);
    //     print_row("fixed CC",             run_fixed_cc(f, 0, 1),      exact);
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random N=" + std::to_string(N),
    //                   run_random(f, 0, 1, N), exact);
    //     }
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random+CC N=" + std::to_string(N),
    //                   run_random_cc(f, 0, 1, N), exact);
    //     }
    // }

    // // 2. Endpoint singularity: ∫01 1/√x dx = 2
    // {
    //     std::cout << "\n=== ∫₀¹ 1/√x dx  (exact = 2, singularity at x=0) ===\n";
    //     auto f = [](double x) { return 1.0 / std::sqrt(x); };
    //     const double exact = 2.0;

    //     print_row("fixed (open 4-pt)",    run_fixed(f, 0, 1),         exact);
    //     print_row("fixed CC",             run_fixed_cc(f, 0, 1),      exact);
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random N=" + std::to_string(N),
    //                   run_random(f, 0, 1, N), exact);
    //     }
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random+CC N=" + std::to_string(N),
    //                   run_random_cc(f, 0, 1, N), exact);
    //     }
    // }

    // // 3. Double singularity: ∫01 ln(x)/√x dx = -4
    // {
    //     std::cout << "\n=== ∫₀¹ ln(x)/√x dx  (exact = -4, singularity at x=0) ===\n";
    //     auto f = [](double x) { return std::log(x) / std::sqrt(x); };
    //     const double exact = -4.0;

    //     print_row("fixed (open 4-pt)",    run_fixed(f, 0, 1),         exact);
    //     print_row("fixed CC",             run_fixed_cc(f, 0, 1),      exact);
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random N=" + std::to_string(N),
    //                   run_random(f, 0, 1, N), exact);
    //     }
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random+CC N=" + std::to_string(N),
    //                   run_random_cc(f, 0, 1, N), exact);
    //     }
    // }

    // // 4. Smooth: ∫01 √(1-x²) dx = π/4
    // {
    //     std::cout << "\n=== ∫₀¹ √(1-x²) dx  (exact = π/4) ===\n";
    //     auto f = [](double x) { return std::sqrt(1.0 - x * x); };
    //     const double exact = M_PI / 4.0;

    //     print_row("fixed (open 4-pt)",    run_fixed(f, 0, 1),         exact);
    //     print_row("fixed CC",             run_fixed_cc(f, 0, 1),      exact);
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random N=" + std::to_string(N),
    //                   run_random(f, 0, 1, N), exact);
    //     }
    //     for (int N : {2, 4, 8, 16}) {
    //         print_row("random+CC N=" + std::to_string(N),
    //                   run_random_cc(f, 0, 1, N), exact);
    //     }
    // }

    // // 5. Find optimal N: sweep across integrals and N values
    // {
    //     std::cout << "\n=== Optimal N sweep (acc=eps=1e-4, random+CC) ===\n";
    //     std::cout << std::left
    //               << std::setw(30) << "integral"
    //               << std::setw(6)  << "N"
    //               << std::setw(12) << "ncalls"
    //               << std::setw(14) << "act_err" << '\n';
    //     std::cout << std::string(62, '-') << '\n';

    //     struct Case { std::string name; Function f; double a, b, exact; };
    //     std::vector<Case> cases = {
    //         {"sqrt(x)",        [](double x){ return std::sqrt(x); },            0, 1,  2.0/3.0},
    //         {"1/sqrt(x)",      [](double x){ return 1.0/std::sqrt(x); },        0, 1,  2.0},
    //         {"ln(x)/sqrt(x)",  [](double x){ return std::log(x)/std::sqrt(x);}, 0, 1, -4.0},
    //     };

    //     for (auto& c : cases) {
    //         for (int N : {2, 4, 8, 16, 32}) {
    //             // average over 5 runs (random integrator)
    //             double total_calls = 0;
    //             double total_err   = 0;
    //             const int REPS = 5;
    //             for (int r = 0; r < REPS; ++r) {
    //                 auto res = run_random_cc(c.f, c.a, c.b, N, 1e-4, 1e-4);
    //                 total_calls += res.ncalls;
    //                 total_err   += std::abs(res.value - c.exact);
    //             }
    //             std::cout << std::left
    //                       << std::setw(30) << c.name
    //                       << std::setw(6)  << N
    //                       << std::setw(12) << static_cast<int>(total_calls / REPS)
    //                       << std::setw(14) << total_err / REPS << '\n';
    //         }
    //         std::cout << '\n';
    //     }
    // }

    return 0;
}
