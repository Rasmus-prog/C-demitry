#include<iostream>
#include<complex>
#include<cmath>
#include"sfuns.h"
int main(){
    using complex = std::complex<double>;
    double ePI =std::exp(std::numbers::pi);
    double PIe =std::pow(std::numbers::pi, std::numbers::e);
    double sqrt2=std::sqrt(2.0);
    double fifthrt2=std::pow(2.0, 1.0/5.0);
    constexpr double  π = 3.14159265358979324;
    constexpr double  E = 2.71828182845904523;
    constexpr complex I = complex(0,1);

    std::cout << "sqrt(2) = " << sqrt2 << std::endl;
    std::cout << "fifth root of 2 = " << fifthrt2 << std::endl;
    std::cout << "e^pi = " << ePI << std::endl;
    std::cout << "pi^e = " << PIe << std::endl;
    std::cout << "log(I)=" << std::log(I)   <<"\n";
    std::cout << "   I^I=" << std::pow(I,I) <<"\n";
    std::cout << "   pi^I=" << std::pow(π,I) <<"\n";
    std::cout << "   e^I=" << std::pow(E,I) <<"\n";
   
    
    std::cout << "My gamma function" << std::endl;
    for (double x = 1; x <= 10; ++x) {
    std::cout << "fgamma(" << x << ") = " << sfuns::fgamma(x) << '\n';
    }
    std::cout << "Standard gamma function" << std::endl;
    for (double x = 1; x <= 10; ++x) {
    std::cout << "std::tgamma(" << x << ") = " << std::tgamma(x) << '\n';
    }

    std::cout << "My log gamma function" << std::endl;
    for (double x = 1; x <= 10; ++x) {
    std::cout << "lngamma(" << x << ") = " << sfuns::lngamma(x) << '\n';
    }
    std::cout << "Standard log gamma function" << std::endl;
    for (double x = 1; x <= 10; ++x) {
    std::cout << "std::lgamma(" << x << ") = " << std::lgamma(x) << '\n';
    }   


    return 0;
}