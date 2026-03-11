#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
#include "lsfit.h"

int main() {
    // Rutherford & Soddy (1902) ThX radioactivity data
    pp::vector t  = {1, 2, 3, 4, 6, 9, 10, 13, 15};
    pp::vector y  = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
    pp::vector dy = {6, 5, 4, 4, 4, 3, 3, 2, 2};
    int n = t.size();

    // Fit y(t) = a * exp(-lambda*t) in logarithmic form:
    //   ln(y) = ln(a) - lambda*t   (linear in {ln(a), lambda})
    //
    // Error propagation: if y = f(z), then delta_z ~ |df/dz|^{-1} * delta_f
    //   Equivalently, delta(ln y) = |d ln(y)/dy| * delta_y = delta_y / y
    pp::vector lny(n), dlny(n);
    for (int i = 0; i < n; i++) {
        lny[i]  = std::log(y[i]);
        dlny[i] = dy[i] / y[i];   // delta(ln y) = delta_y / y
    }

    // Basis functions for the linear fit of ln(y):
    //   f0(t) = 1   ->  coefficient c[0] = ln(a)
    //   f1(t) = t   ->  coefficient c[1] = -lambda
    std::vector<std::function<double(double)>> fs {
        [](double /*z*/) { return 1.0; },
        [](double z)     { return z;   }
    };

    auto [c, C] = lsfit(fs, t, lny, dlny);

    double ln_a     = c[0];
    double a        = std::exp(ln_a);
    double lambda   = -c[1];
    double half_life = std::log(2.0) / lambda;

    // Uncertainties from diagonal of covariance matrix
    double sigma_c0 = std::sqrt(C[0, 0]);
    double sigma_c1 = std::sqrt(C[1, 1]);
    double sigma_a  = a * sigma_c0;          // delta(a) = a * delta(ln a)
    double sigma_lam = sigma_c1;             // lambda = -c[1], same magnitude
    // Error propagation: T_{1/2} = ln(2)/lambda
    //   delta(T_{1/2}) = (ln(2)/lambda^2) * delta(lambda) = T_{1/2}/lambda * sigma_lam
    double sigma_T  = half_life / lambda * sigma_lam;

    // Modern accepted half-life of 224Ra (= what was called ThX)
    const double half_life_modern = 3.6319; // days

    std::cout << "=== Least-squares fit: ln(y) = c0 + c1*t ===\n";
    std::cout << "c0 = ln(a) = " << ln_a   << " +/- " << sigma_c0 << "\n";
    std::cout << "c1 = -lambda = " << c[1] << " +/- " << sigma_c1 << "\n";
    std::cout << "\n";
    C.print("Covariance matrix C:");
    std::cout << "\na      = " << a     << " +/- " << sigma_a   << "\n";
    std::cout << "lambda = " << lambda  << " +/- " << sigma_lam << " day^-1\n";
    std::cout << "T_1/2  = " << half_life << " +/- " << sigma_T << " days\n";
    std::cout << "Modern value for 224Ra: " << half_life_modern << " days\n";
    double delta = half_life - half_life_modern;
    double pull = delta / sigma_T;
    std::cout << "Deviation: " << delta << " days"
              << "  (" << pull << " sigma)\n";
    bool agrees_within_uncertainty = pp::approx(half_life, half_life_modern);
    std::cout << "Does it agree with the modern value within the estimated uncertainty? "
              << (agrees_within_uncertainty ? "Yes" : "No") << "\n";

    // Write data file: t  y  dy  (for gnuplot errorbars)
    {
        std::ofstream f("plot.data");
        f << "# t   y   dy\n";
        for (int i = 0; i < n; i++)
            f << t[i] << "  " << y[i] << "  " << dy[i] << "\n";
    }

    // Write model curves for plotting quality of coefficient uncertainties.
    // Columns: x  best  (++), (+-), (-+), (--)
    {
        std::ofstream f("plot.curves.data");
        f << "# x  y_best  y_pp  y_pm  y_mp  y_mm\n";

        double tmin = t[0], tmax = t[0];
        for (int i = 1; i < n; i++) {
            if (t[i] < tmin) tmin = t[i];
            if (t[i] > tmax) tmax = t[i];
        }
        tmin = std::max(0.0, tmin - 0.5);
        tmax += 0.5;

        const int samples = 300;
        for (int i = 0; i < samples; i++) {
            double x = tmin + (tmax - tmin) * i / (samples - 1.0);

            auto model = [x, &fs](const pp::vector& coeffs) {
                double s = 0.0;
                for (int k = 0; k < coeffs.size(); k++) s += coeffs[k] * fs[k](x);
                return std::exp(s); // we fitted ln(y)
            };

            pp::vector c_pp = c, c_pm = c, c_mp = c, c_mm = c;
            c_pp[0] += sigma_c0; c_pp[1] += sigma_c1;
            c_pm[0] += sigma_c0; c_pm[1] -= sigma_c1;
            c_mp[0] -= sigma_c0; c_mp[1] += sigma_c1;
            c_mm[0] -= sigma_c0; c_mm[1] -= sigma_c1;

            f << x << "  "
              << model(c)    << "  "
              << model(c_pp) << "  "
              << model(c_pm) << "  "
              << model(c_mp) << "  "
              << model(c_mm) << "\n";
        }
    }

    return 0;
}
