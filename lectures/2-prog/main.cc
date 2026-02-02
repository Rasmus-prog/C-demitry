#include<iostream>
#include"hello.h"
#include"sfuns.h"
int main(){
    hello();
    double x = 1;
    double y = sfuns::fgamma(x);
    std::cout << "fgamma(" << x << ") = " << y << "\n";
    return 0;
}