#include "lcg.h"

lcg::lcg(result_type seed, result_type a, result_type c)
    : state(seed), a(a), c(c) {}

lcg::result_type lcg::operator()() {
    state = a * state + c;
    return state;
}