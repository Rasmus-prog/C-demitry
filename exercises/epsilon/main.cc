#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>

bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9) {
    double diff = std::abs(a - b);
    if (diff <= acc) return true;  // absolute tolerance
    double max_ab = std::max(std::abs(a), std::abs(b));
    return diff <= eps * max_ab;   // relative tolerance
}

int main()
{
    std::cout << "Machine epsilon" << "\n";

    float       f=1.0f; while((float)      (1.0f+f) != 1.0f){f/=2.0f;} f*=2.0f;
    double      d=1.0d; while((double)     (1.0d+d) != 1.0d){d/=2.0d;} d*=2.0d;
    long double l=1.0L; while((long double)(1.0L+l) != 1.0L){l/=2.0L;} l*=2.0L;
    std::cout << "      float eps=" << f << "\n";
    std::cout << "     double eps=" << d << "\n";
    std::cout << "long double eps=" << l << "\n";


    std::cout << std::numeric_limits<float>::epsilon() << "\n";
    std::cout << std::numeric_limits<double>::epsilon() << "\n";
    std::cout << std::numeric_limits<long double>::epsilon() << "\n";
    
    std::cout << std::pow(2,-52) << "\n";
    std::cout << std::pow(2,-23) << "\n";

    std::cout << "Non-commutativity of addition" << "\n";

    double epsilon=std::pow(2,-52);
    double tiny=epsilon/2;
    double a=1+tiny+tiny;
    double b=tiny+tiny+1;

    std::cout << "a=" << a << "\n";
    std::cout << "b=" << b << "\n";

    std::cout << "a==b ? " << (a==b ? "true":"false") << "\n"; // the order of the additions matters, something big plus something small removes the small, but something small plus something big does not remove the small
    std::cout << "a>1  ? " << (a>1  ? "true":"false") << "\n"; // something big plus something small removes the small
    std::cout << "b>1  ? " << (b>1  ? "true":"false") << "\n"; // something small plus something big does not remove the small

    std::cout << std::fixed << std::setprecision(17);
    std::cout << "       tiny=" << tiny << "\n";
    std::cout << "1+tiny+tiny=" << a << "\n";
    std::cout << "tiny+tiny+1=" << b << "\n";
    // we see what stated earlier


    std::cout << "Comparing doubles: introduction" << "\n";

    double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
    double d2 = 8*0.1;
    std::cout << "d1==d2? " << (d1==d2 ? "true":"false") << "\n"; 


    std::cout << "d1=" << d1 << "\n";
    std::cout << "d2=" << d2 << "\n";
    

    std::cout << "Comparing doubles: the task" << "\n";
    std::cout << "approx(d1, d2)? " << (approx(d1, d2) ? "true":"false") << "\n";


    return 0;
}

