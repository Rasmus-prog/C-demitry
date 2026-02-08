#include <iostream>
#include <cmath>
#include <random>
#include "vec.h"

int main(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    vec u(dis(gen), dis(gen), dis(gen));
    vec v(dis(gen), dis(gen), dis(gen));
    
    std::cout << "u=" << u << "\n";
    std::cout << "v=" << v << "\n\n";

    vec t;

    // Test unary negation
    std::cout << "Test: unary negation operator(-)" << "\n";
    t = vec(-u.x, -u.y, -u.z);
    std::cout << "-u = " << (-u) << "\n";
    std::cout << "t  = " << t << "\n";
    if(approx(t, -u)) {
        std::cout << "test 'unary negation' passed" << "\n\n";
    } else {
        std::cout << "test 'unary negation' FAILED" << "\n\n";
    }

    // Test subtraction
    std::cout << "Test: binary subtraction operator(-)" << "\n";
    t = vec(u.x - v.x, u.y - v.y, u.z - v.z);
    std::cout << "u-v = " << (u-v) << "\n";
    std::cout << "t   = " << t << "\n";
    if(approx(t, u-v)) {
        std::cout << "test 'subtraction' passed" << "\n\n";
    } else {
        std::cout << "test 'subtraction' FAILED" << "\n\n";
    }

    // Test addition
    std::cout << "Test: binary addition operator(+)" << "\n";
    t = vec(u.x + v.x, u.y + v.y, u.z + v.z);
    std::cout << "u+v = " << (u+v) << "\n";
    std::cout << "t   = " << t << "\n";
    if(approx(t, u+v)) {
        std::cout << "test 'addition' passed" << "\n\n";
    } else {
        std::cout << "test 'addition' FAILED" << "\n\n";
    }

    // Test scalar multiplication
    std::cout << "Test: scalar multiplication operator(*)" << "\n";
    double c = dis(gen);
    t = vec(u.x*c, u.y*c, u.z*c);
    vec tmp1 = u*c;
    vec tmp2 = c*u;
    std::cout << "u*c = " << tmp1 << "\n";
    std::cout << "t   = " << t << "\n";
    if(approx(t, u*c)) {
        std::cout << "test 'multiplication (u*c)' passed" << "\n";
    } else {
        std::cout << "test 'multiplication (u*c)' FAILED" << "\n";
    }
    std::cout << "c*u = " << tmp2 << "\n";
    if(approx(t, c*u)) {
        std::cout << "test 'multiplication (c*u)' passed" << "\n\n";
    } else {
        std::cout << "test 'multiplication (c*u)' FAILED" << "\n\n";
    }

    // Test scalar division
    std::cout << "Test: scalar division operator(/)" << "\n";
    double d = dis(gen) + 0.5; // ensure non-zero
    t = vec(u.x/d, u.y/d, u.z/d);
    std::cout << "u/" << d << " = " << (u/d) << "\n";
    std::cout << "t           = " << t << "\n";
    if(approx(t, u/d)) {
        std::cout << "test 'division' passed" << "\n\n";
    } else {
        std::cout << "test 'division' FAILED" << "\n\n";
    }

    // Test dot product
    std::cout << "Test: dot product operator(%)" << "\n";
    double dot_product = u.x*v.x + u.y*v.y + u.z*v.z;
    double result = u % v;
    std::cout << "u % v = " << result << "\n";
    std::cout << "dot   = " << dot_product << "\n";
    if(approx(dot_product, result)) {
        std::cout << "test 'dot product' passed" << "\n\n";
    } else {
        std::cout << "test 'dot product' FAILED" << "\n\n";
    }

    // Test compound assignment operators
    std::cout << "Test: compound assignment operators" << "\n";
    vec a = u;
    vec b = v;
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
        std::cout << "test '-=' passed" << "\n";
    } else {
        std::cout << "test '-=' FAILED" << "\n";
    }
    
    a = u;
    a *= c;
    std::cout << "u *= " << c << ": " << a << " (expected " << (u*c) << ")" << "\n";
    if(approx(a, u*c)) {
        std::cout << "test '*=' passed" << "\n";
    } else {
        std::cout << "test '*=' FAILED" << "\n";
    }
    
    a = u;
    a /= d;
    std::cout << "u /= " << d << ": " << a << " (expected " << (u/d) << ")" << "\n";
    if(approx(a, u/d)) {
        std::cout << "test '/=' passed" << "\n\n";
    } else {
        std::cout << "test '/=' FAILED" << "\n\n";
    }

    // Test print method
    std::cout << "Test: print method" << "\n";
    u.print("Vector u: ");
    v.print("Vector v: ");

    return 0;
}
