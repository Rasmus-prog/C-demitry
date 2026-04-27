#pragma once
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>
class lcg {
public:
    using result_type = std::uint32_t;

    lcg(result_type seed, result_type a = 1664525u, result_type c = 1013904223u);
    static constexpr result_type min() { return 0u; }
    static constexpr result_type max() { return 0xffffffffu; }
    result_type operator()();

private:
    result_type state, a, c;
};