#include<iostream>
#include"sfuns.h"
int main(){
    double sqrt2=std::sqrt(2.0);
    std::cout << "sqrt(2) = " << sqrt2 << std::endl;
    double fifthrt2=std::pow(2.0, 1.0/5.0);
    std::cout << "fifth root of 2 = " << fifthrt2 << std::endl;
    double ePI =std::exp(std::numbers::pi);
    std::cout << "e^pi = " << ePI << std::endl;
    double PIe =std::pow(std::numbers::pi, std::numbers::e);
    std::cout << "pi^e = " << PIe << std::endl;
   
    

    for (double x = 1; x <= 10; ++x) {
    std::cout << "fgamma(" << x << ") = " << sfuns::fgamma(x) << '\n';
    }
    return 0;
}