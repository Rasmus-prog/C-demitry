#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include "sfuns.h"


int main() {

    // ∫01 dx √(x) = 2/3 ,
	int ncalls1 = 0;
    auto f1 = [&ncalls1](double x) { 
		++ncalls1;
        return std::sqrt(x); 
    };
    const auto [result1, error1] = integrate(f1, 0.0, 1.0);
    
    std::cout << "f1 ∫01 dx √(x) = 2/3" << "\n";
    std::cout << "result = " << result1 << '\n';
	std::cout << "error estimate = " << error1 << '\n';
	std::cout << "ncalls = " << ncalls1 << '\n';
    // GAI start
    std::cout << "result is within 1e-3 of 2/3: " << (std::abs(result1 - 2.0/3.0) < 1e-3 ? "true" : "false") << '\n';
    // GAI end



    // ∫01 dx 1/√(x) = 2 ,
	int ncalls2 = 0;
	auto f2 = [&ncalls2](double x) {
		++ncalls2;
		return 1.0 / std::sqrt(x);
	};
    const auto [result2, error2] = integrate(f2, 0.0, 1.0);
    std::cout << "f2 ∫01 dx 1/√(x) = 2" << "\n";
    std::cout << "result = " << result2 << '\n';
	std::cout << "error estimate = " << error2 << '\n';
	std::cout << "ncalls = " << ncalls2 << '\n';
    std::cout << "result is within 1e-3 of 2: " << (std::abs(result2 - 2.0) < 1e-3 ? "true" : "false") << '\n';

    // ∫01 dx √(1-x²) = π/4
    int ncalls3 = 0;
    auto f3 = [&ncalls3](double x) {
        ++ncalls3;
        return std::sqrt(1.0 - (x * x));
    };
    const auto [result3, error3] = integrate(f3, 0.0, 1.0);
    std::cout << "f3 ∫01 dx √(1-x²) = π/4" << "\n";
    std::cout << "result = " << result3 << '\n';
    std::cout << "error estimate = " << error3 << '\n';
    std::cout << "ncalls = " << ncalls3 << '\n';
    std::cout << "result is within 1e-3 of π/4: " << (std::abs(result3 - M_PI/4.0) < 1e-3 ? "true" : "false") << '\n';

    // ∫01 dx ln(x)/√(x) = -4
	int ncalls4 = 0;
	auto f4 = [&ncalls4](double x) {
		++ncalls4;
		return std::log(x) / std::sqrt(x);
	};
    const auto [result4, error4] = integrate(f4, 0.0, 1.0);
    std::cout << "f4 ∫01 dx ln(x)/√(x) = -4" << "\n";
    std::cout << "result = " << result4 << '\n';
	std::cout << "error estimate = " << error4 << '\n';
	std::cout << "ncalls = " << ncalls4 << '\n';
    std::cout << "result is within 1e-3 of -4: " << (std::abs(result4 + 4.0) < 1e-3 ? "true" : "false") << '\n';

    std::ofstream errf_out("erf_out.data");
    errf_out << "# x my_erf tabulated_erf abs_diff\n";
    for (double x = 0.0; x <= 3.0; x += 0.1) {
        const double my_result = errf_func(x, 1e-3, 1e-3);
        const double tabulated = std::erf(x);
        errf_out << x << ' ' << my_result << ' ' << tabulated << ' ' << std::abs(my_result - tabulated) << '\n';
    }
    errf_out.close();


    std::ofstream acc_out("erf_acc.data");
    acc_out << "# acc abs(my_erf(1)-exact)\n";
    const double exact_result = 0.84270079294971486934;
    for (double acc = 0.1; acc >= 1e-6; acc /= 10.0) {
        const double result = errf_func(1.0, acc, 0.0);
        acc_out << acc << ' ' << std::abs(result - exact_result) << '\n';
    }

	acc_out.close();

    std::cout << "erf(1) with acc=1e-6, eps=0: " << errf_func(1.0, 1e-6, 0.0) << '\n';
    std::cout << "reference erf(1): " << exact_result << '\n';

    std::cout << "transformation results " << "\n";
    //  ∫01 dx 1/√(x) = 2 , but with the Clenshaw–Curtis variable transformation
    int ncalls2_pre = ncalls2;
    ncalls2 = 0;
    const auto [result5, error5] = integrate_clenshaw_curtis(f2, 0.0, 1.0);
    std::cout << "f2 ∫01 dx 1/√(x) = 2 with the Clenshaw–Curtis variable transformation" << "\n";
    std::cout << "result = " << result5 << '\n';
	std::cout << "error estimate = " << error5 << '\n';
    std::cout << "ncalls = " << ncalls2 << '\n';
    std::cout << "ncalls with transformation is " << ncalls2 << " vs " << ncalls2_pre << " without transformation" << "\n";
    std::cout << "result is within 1e-3 of 2: " << (std::abs(result5 - 2.0) < 1e-3 ? "true" : "false") << '\n';

    // ∫01 dx ln(x)/√(x) = -4 . but with the Clenshaw–Curtis variable transformation
    int ncalls4_pre = ncalls4;
    ncalls4 = 0;
    const auto [result6, error6] = integrate_clenshaw_curtis(f4, 0.0, 1.0);
    std::cout << "f4 ∫01 dx ln(x)/√(x) = -4 with the Clenshaw–Curtis variable transformation" << "\n";
    std::cout << "result = " << result6 << '\n';
	std::cout << "error estimate = " << error6 << '\n';
    std::cout << "ncalls = " << ncalls4 << '\n';
    std::cout << "ncalls with transformation is " << ncalls4 << " vs " << ncalls4_pre << " without transformation" << "\n";
    std::cout << "result is within 1e-3 of -4: " << (std::abs(result6 + 4.0) < 1e-3 ? "true" : "false") << '\n';

    // Test your implementation on some (converging)
    // ∫0∞ dx e^(-x) = 1
    int ncalls9 = 0;
    auto f5 = [&ncalls9](double x) {
        ++ncalls9;
        return std::exp(-x);
    };
    const auto [result7, error7] = integrate(f5, 0.0, std::numeric_limits<double>::infinity(), 1e-4, 1e-4);
    std::cout << "f5 ∫0∞ dx e^(-x) = 1" << "\n";
    std::cout << "result = " << result7 << '\n';
	std::cout << "error estimate = " << error7 << '\n';
    std::cout << "ncalls = " << ncalls9 << '\n';
    std::cout << "result is within 1e-3 of 1: " << (std::abs(result7 - 1.0) < 1e-3 ? "true" : "false") << '\n';
	return 0;
}
