#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>

using Function = std::function<double(double)>;

double integrate(const Function& f,
				 double a,
				 double b,
				 double acc = 1e-3,
				 double eps = 1e-3,
				 std::optional<double> f2 = std::nullopt,
				 std::optional<double> f3 = std::nullopt) {
	const double h = b - a;

	const double x1 = a + h / 6.0;
	const double x2 = a + 2.0 * h / 6.0;
	const double x3 = a + 4.0 * h / 6.0;
	const double x4 = a + 5.0 * h / 6.0;

	const double v1 = f(x1);
	const double v2 = f2.has_value() ? *f2 : f(x2);
	const double v3 = f3.has_value() ? *f3 : f(x3);
	const double v4 = f(x4);

	const double Q = (2*v1+v2+v3+2*v4)/6*(b-a);
	const double q = (  v1+v2+v3+  v4)/4*(b-a);

	const double err = std::abs(Q - q);
	const double tol = acc + eps * std::abs(Q);

	if (err < tol) {
		return Q;
	}

	const double m = (a + b) / 2.0;
	const double nextAcc = acc / std::sqrt(2.0);

	return integrate(f, a, m, nextAcc, eps, v1, v2) +
		   integrate(f, m, b, nextAcc, eps, v3, v4);
}


double errf_func(double x, double acc = 1e-3, double eps = 1e-3) {
    if (x < 0.0) {
        return -errf_func(-x, acc, eps);
    } else if (x <= 1.0) {
        auto f = [](double t) { return std::exp(-t*t); };
        return 2.0 / std::sqrt(M_PI) * integrate(f, 0.0, x, acc, eps);
    } else {
        auto f = [x](double t) { return std::exp(-(x+(1.0-t)/t)*(x+(1.0-t)/t))/t/t; };
        return 1.0 - 2.0 / std::sqrt(M_PI) * integrate(f, 0.0, 1.0, acc, eps);
    }



}

int main() {

    // ∫01 dx √(x) = 2/3 ,
	int ncalls = 0;
    auto f1 = [&ncalls](double x) { 
		++ncalls;
        return std::sqrt(x); 
    };
    const double result1 = integrate(f1, 0.0, 1.0);
    
    std::cout << "f1 ∫01 dx √(x) = 2/3" << "\n";
    std::cout << "result = " << result1 << '\n';
	std::cout << "ncalls = " << ncalls << '\n';
    // GAI start
    std::cout << "result is within 1e-3 of 2/3: " << (std::abs(result1 - 2.0/3.0) < 1e-3 ? "true" : "false") << '\n';
    // GAI end



    // ∫01 dx 1/√(x) = 2 ,
	ncalls = 0;
	auto f2 = [&ncalls](double x) {
		++ncalls;
		return 1.0 / std::sqrt(x);
	};
    const double result2 = integrate(f2, 0.0, 1.0);
    std::cout << "f2 ∫01 dx 1/√(x) = 2" << "\n";
    std::cout << "result = " << result2 << '\n';
	std::cout << "ncalls = " << ncalls << '\n';
    std::cout << "result is within 1e-3 of 2: " << (std::abs(result2 - 2.0) < 1e-3 ? "true" : "false") << '\n';


    // ∫01 dx √(1-x²) = π/4
    ncalls = 0;
    auto f3 = [&ncalls](double x) {
        ++ncalls;
        return std::sqrt(1.0 - (x * x));
    };
    const double result3 = integrate(f3, 0.0, 1.0);
    std::cout << "f3 ∫01 dx √(1-x²) = π/4" << "\n";
    std::cout << "result = " << result3 << '\n';
    std::cout << "ncalls = " << ncalls << '\n';
    std::cout << "result is within 1e-3 of π/4: " << (std::abs(result3 - M_PI/4.0) < 1e-3 ? "true" : "false") << '\n';


    // ∫01 dx ln(x)/√(x) = -4
	ncalls = 0;
	auto f4 = [&ncalls](double x) {
		++ncalls;
		return std::log(x) / std::sqrt(x);
	};
    const double result4 = integrate(f4, 0.0, 1.0);
    std::cout << "f4 ∫01 dx ln(x)/√(x) = -4" << "\n";
    std::cout << "result = " << result4 << '\n';
	std::cout << "ncalls = " << ncalls << '\n';
    std::cout << "result is within 1e-3 of -4: " << (std::abs(result4 + 4.0) < 1e-3 ? "true" : "false") << '\n';

    ncalls = 0;
    std::ofstream errf_out("erf_out.data");
    errf_out << "# acc result\n";
    for (double x = 0.0; x <= 3.0; x += 0.1) {
        const double result = errf_func(x, 1e-3, 1e-3);
        errf_out << x << ' ' << result << '\n';
    }
    errf_out.close();


    std::ofstream acc_out("erf_acc.data");
    acc_out << "# acc result - exact result of erf(1)\n";
    double exact_result = errf_func(1.0, 1e-3, 1e-3);
    for (double acc = 0.1; acc >= 1e-6; acc /= 10.0) {
        ncalls = 0;
        const double result = errf_func(1.0, acc, 0.0);
        acc_out << acc << ' ' << std::abs(result - exact_result) << '\n';
    }

	acc_out.close();
	return 0;
}
