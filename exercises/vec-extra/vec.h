#pragma once
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>

template<typename T>
struct vec {
    T x, y, z;

    // constructors
    vec(T x, T y, T z) : x(x), y(y), z(z) {}  // parameterized
    vec() : vec(T(0), T(0), T(0)) {}           // default
    vec(const vec&) = default;                 // copy
    vec(vec&&) = default;                      // move
    ~vec() = default;                          // destructor

    // assignment
    vec& operator=(const vec&) = default;      // copy assignment
    vec& operator=(vec&&) = default;           // move assignment

    // arithmetic
    vec& operator+=(const vec&);
    vec& operator-=(const vec&);
    vec& operator*=(T);
    vec& operator/=(T);

    // utility
    void set(T a, T b, T c) { x = a; y = b; z = c; }
    void print(const std::string& s = "") const;              // for debugging

    // stream output
    template<typename U>
    friend std::ostream& operator<<(std::ostream&, const vec<U>&);
};

// non-member operators
template<typename T>
vec<T> operator-(const vec<T>&);

template<typename T>
vec<T> operator+(const vec<T>&, const vec<T>&);

template<typename T>
vec<T> operator-(const vec<T>&, const vec<T>&);

template<typename T>
vec<T> operator*(const vec<T>&, T);

template<typename T>
vec<T> operator*(T, const vec<T>&);

template<typename T>
vec<T> operator/(const vec<T>&, T);

template<typename T>
T operator%(const vec<T>&, const vec<T>&);  // dot product

// approximate equality
template<typename T>
bool approx(const vec<T>&, const vec<T>&, double acc = 1e-6, double eps = 1e-6);

bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9);
bool approx(float a, float b, float acc = 1e-6, float eps = 1e-6);

template<typename T>
bool approx(std::complex<T> a, std::complex<T> b, double acc = 1e-6, double eps = 1e-6);


// Compound assignment operators
template<typename T>
vec<T>& vec<T>::operator+=(const vec<T>& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

template<typename T>
vec<T>& vec<T>::operator-=(const vec<T>& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

template<typename T>
vec<T>& vec<T>::operator*=(T s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

template<typename T>
vec<T>& vec<T>::operator/=(T s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

// Non-member operators
template<typename T>
vec<T> operator-(const vec<T>& v) {
    return vec<T>(-v.x, -v.y, -v.z);
}

template<typename T>
vec<T> operator+(const vec<T>& u, const vec<T>& v) {
    return vec<T>(u.x + v.x, u.y + v.y, u.z + v.z);
}

template<typename T>
vec<T> operator-(const vec<T>& u, const vec<T>& v) {
    return vec<T>(u.x - v.x, u.y - v.y, u.z - v.z);
}

template<typename T>
vec<T> operator*(const vec<T>& v, T s) {
    return vec<T>(v.x * s, v.y * s, v.z * s);
}

template<typename T>
vec<T> operator*(T s, const vec<T>& v) {
    return vec<T>(v.x * s, v.y * s, v.z * s);
}

template<typename T>
vec<T> operator/(const vec<T>& v, T s) {
    return vec<T>(v.x / s, v.y / s, v.z / s);
}

// Dot product
template<typename T>
T operator%(const vec<T>& u, const vec<T>& v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// Utility
template<typename T>
void vec<T>::print(const std::string& s) const {
    std::cout << s << "(" << x << ", " << y << ", " << z << ")" << std::endl;
}

// Stream output
template<typename T>
std::ostream& operator<<(std::ostream& os, const vec<T>& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

// Approximate equality for vec
template<typename T>
bool approx(const vec<T>& u, const vec<T>& v, double acc, double eps) {
    auto diff_x = u.x - v.x;
    auto diff_y = u.y - v.y;
    auto diff_z = u.z - v.z;
    
    auto mag_diff_sq = std::norm(diff_x) + std::norm(diff_y) + std::norm(diff_z);
    auto diff = std::sqrt(mag_diff_sq);
    
    if (diff <= acc) return true;
    
    auto mag_u_sq = std::norm(u.x) + std::norm(u.y) + std::norm(u.z);
    auto mag_v_sq = std::norm(v.x) + std::norm(v.y) + std::norm(v.z);
    auto mag_u = std::sqrt(mag_u_sq);
    auto mag_v = std::sqrt(mag_v_sq);
    auto max_mag = std::max(mag_u, mag_v);
    
    return diff <= eps * max_mag;
}

// Approximate equality for std::complex
template<typename T>
bool approx(std::complex<T> a, std::complex<T> b, double acc, double eps) {
    auto diff = std::abs(a - b);
    if (diff <= acc) return true;
    auto max_ab = std::max(std::abs(a), std::abs(b));
    return diff <= eps * max_ab;
}

// Approximate equality for doubles
inline bool approx(double a, double b, double acc, double eps) {
    double diff = std::abs(a - b);
    if (diff <= acc) return true;
    double max_ab = std::max(std::abs(a), std::abs(b));
    return diff <= eps * max_ab;
}

// Approximate equality for floats
inline bool approx(float a, float b, float acc, float eps) {
    float diff = std::abs(a - b);
    if (diff <= acc) return true;
    float max_ab = std::max(std::abs(a), std::abs(b));
    return diff <= eps * max_ab;
}

