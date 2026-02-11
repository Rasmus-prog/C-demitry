#include <iostream>
#include <cmath>
#include <random>
#include <complex>
#include "vec.h"

// Generic test function for any vec<T> type
template<typename T>
void run_tests(const std::string& type_name, vec<T> u, vec<T> v) {
    std::cout << "\n========== TESTING WITH " << type_name << " ==========\n\n";
    std::cout << "u=" << u << "\n";
    std::cout << "v=" << v << "\n\n";

    vec<T> t;

    // Test unary negation
    std::cout << "Test: unary negation operator(-)" << "\n";
    t = vec<T>(-u.x, -u.y, -u.z);
    std::cout << "-u = " << (-u) << "\n";
    std::cout << "t  = " << t << "\n";
    if(approx(t, -u)) {
        std::cout << "test 'unary negation' passed" << "\n\n";
    } else {
        std::cout << "test 'unary negation' FAILED" << "\n\n";
    }

    // Test subtraction
    std::cout << "Test: binary subtraction operator(-)" << "\n";
    t = vec<T>(u.x - v.x, u.y - v.y, u.z - v.z);
    std::cout << "u-v = " << (u-v) << "\n";
    std::cout << "t   = " << t << "\n";
    if(approx(t, u-v)) {
        std::cout << "test 'subtraction' passed" << "\n\n";
    } else {
        std::cout << "test 'subtraction' FAILED" << "\n\n";
    }

    // Test addition
    std::cout << "Test: binary addition operator(+)" << "\n";
    t = vec<T>(u.x + v.x, u.y + v.y, u.z + v.z);
    std::cout << "u+v = " << (u+v) << "\n";
    std::cout << "t   = " << t << "\n";
    if(approx(t, u+v)) {
        std::cout << "test 'addition' passed" << "\n\n";
    } else {
        std::cout << "test 'addition' FAILED" << "\n\n";
    }

    // Dot product
    std::cout << "Test: dot product operator(%)" << "\n";
    auto dot_result = u % v;
    std::cout << "u % v = " << dot_result << "\n";
    std::cout << "test 'dot product' passed" << "\n\n";

    // Test compound assignment operators
    std::cout << "Test: compound assignment operators" << "\n";
    vec<T> a = u;
    a += v;
    std::cout << "u += v: " << a << " (expected " << (u+v) << ")" << "\n";
    if(approx(a, u+v)) {
        std::cout << "test '+=' passed" << "\n";
    } else {
        std::cout << "test '+=' FAILED" << "\n";
    }
    
    a = u;
    a -= v;
    std::cout << "u -= v: " << a << " (expected " << (u-v) << ")" << "\n";
    if(approx(a, u-v)) {
        std::cout << "test '-=' passed" << "\n\n";
    } else {
        std::cout << "test '-=' FAILED" << "\n\n";
    }

    // Test print method
    std::cout << "Test: print method" << "\n";
    u.print("Vector u: ");
    v.print("Vector v: ");
}

int main(){
    std::random_device rd;
    std::mt19937 gen(rd());

    // Test with DOUBLE
    std::uniform_real_distribution<> dis_double(0.0, 1.0);
    vec<double> ud(dis_double(gen), dis_double(gen), dis_double(gen));
    vec<double> vd(dis_double(gen), dis_double(gen), dis_double(gen));
    run_tests("DOUBLE", ud, vd);

    // Test with FLOAT
    std::uniform_real_distribution<float> dis_float(0.0f, 1.0f);
    vec<float> uf(dis_float(gen), dis_float(gen), dis_float(gen));
    vec<float> vf(dis_float(gen), dis_float(gen), dis_float(gen));
    run_tests("FLOAT", uf, vf);

    // Test with COMPLEX<DOUBLE>
    std::uniform_real_distribution<> dis_complex(0.0, 1.0);
    vec<std::complex<double>> uc(
        std::complex<double>(dis_complex(gen), dis_complex(gen)),
        std::complex<double>(dis_complex(gen), dis_complex(gen)),
        std::complex<double>(dis_complex(gen), dis_complex(gen))
    );
    vec<std::complex<double>> vc(
        std::complex<double>(dis_complex(gen), dis_complex(gen)),
        std::complex<double>(dis_complex(gen), dis_complex(gen)),
        std::complex<double>(dis_complex(gen), dis_complex(gen))
    );
    run_tests("COMPLEX<DOUBLE>", uc, vc);

    // Test with INT
    std::uniform_int_distribution<int> dis_int(1, 10);
    vec<int> ui(dis_int(gen), dis_int(gen), dis_int(gen));
    vec<int> vi(dis_int(gen), dis_int(gen), dis_int(gen));
    run_tests("INT", ui, vi);

    // Test with LONG (demonstrates any scalar type that supports arithmetic)
    std::uniform_int_distribution<long> dis_long(1, 10);
    vec<long> ul(dis_long(gen), dis_long(gen), dis_long(gen));
    vec<long> vl(dis_long(gen), dis_long(gen), dis_long(gen));
    run_tests("LONG", ul, vl);

    // Test with LONG DOUBLE (demonstrates any scalar type that supports arithmetic)
    std::uniform_real_distribution<long double> dis_longdbl(0.0L, 1.0L);
    vec<long double> uld(dis_longdbl(gen), dis_longdbl(gen), dis_longdbl(gen));
    vec<long double> vld(dis_longdbl(gen), dis_longdbl(gen), dis_longdbl(gen));
    run_tests("LONG DOUBLE", uld, vld);

    return 0;
}
