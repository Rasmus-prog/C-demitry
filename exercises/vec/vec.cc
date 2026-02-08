#include "vec.h"
#include <cmath>
#include <iomanip>

// Compound assignment operators
vec& vec::operator+=(const vec& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

vec& vec::operator-=(const vec& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

vec& vec::operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

vec& vec::operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

// Non-member operators
vec operator-(const vec& v) {
    return vec(-v.x, -v.y, -v.z);
}

vec operator+(const vec& u, const vec& v) {
    return vec(u.x + v.x, u.y + v.y, u.z + v.z);
}

vec operator-(const vec& u, const vec& v) {
    return vec(u.x - v.x, u.y - v.y, u.z - v.z);
}

vec operator*(const vec& v, double s) {
    return vec(v.x * s, v.y * s, v.z * s);
}

vec operator*(double s, const vec& v) {
    return vec(v.x * s, v.y * s, v.z * s);
}

vec operator/(const vec& v, double s) {
    return vec(v.x / s, v.y / s, v.z / s);
}

// Dot product
double operator%(const vec& u, const vec& v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// Utility
void vec::print(const std::string& s) const {
    std::cout << s << "(" << x << ", " << y << ", " << z << ")" << std::endl;
}

// Stream output
std::ostream& operator<<(std::ostream& os, const vec& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

// Approximate equality
bool approx(const vec& u, const vec& v, double acc, double eps) {
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;
    double diff = std::sqrt(dx * dx + dy * dy + dz * dz);
    
    if (diff <= acc) return true;  // absolute tolerance
    
    double mag_u = std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    double mag_v = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    double max_mag = std::max(mag_u, mag_v);
    
    return diff <= eps * max_mag;  // relative tolerance
}

// Approximate equality for doubles
bool approx(double a, double b, double acc, double eps) {
    double diff = std::abs(a - b);
    if (diff <= acc) return true;  // absolute tolerance
    double max_ab = std::max(std::abs(a), std::abs(b));
    return diff <= eps * max_ab;   // relative tolerance
}
