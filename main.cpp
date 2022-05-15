#include <array>
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <iostream>

// #include "cuPC-S.h"

const int NUMBER_OF_LEVELS = 50;

auto std_normal_qnorm(const double p) -> double
{
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, p);
}

auto threshold_array(const int n, const double alpha) -> std::array<double, NUMBER_OF_LEVELS>
{
    std::array<double, NUMBER_OF_LEVELS> thr{0.0};
    // my loop range is exclusive of the last i, so I don't subtract one here
    const int n_thr = (NUMBER_OF_LEVELS < (n - 3)) ? NUMBER_OF_LEVELS : n;
    const double half = 0.5;
    for (size_t i = 0; i < n_thr; i++) {
        thr[i] = abs(std_normal_qnorm(half * alpha) / sqrt(n - i - 3));
    }
    return thr;
}

void call_skeleton()
{
    std::cout << "starting `call skeleton`" << std::endl;
    const int n = 10;
    int p = 2;
    const double alpha = 0.05;
    int max_level = 14;
    const size_t sepset_size = p * p * 14;
    std::array<double, 4> pmax{0.0};
    std::array<int, 4> G = {0, 1, 1, 0};
    std::array<double, 4> C = {1.0, 0.4, 0.4, 1.0};
    std::array<double, NUMBER_OF_LEVELS> Th = threshold_array(n, alpha);
    int sepset[sepset_size];
    memset(sepset, 0, sizeof sepset);  // but using memset is easy
    std::cout << "calling Skeleton" << std::endl;
    // Skeleton(C.data(), &p, G.data(), Th.data(), 0, &max_level, pmax.data(), sepset);
    std::cout << "done with call" << std::endl;
}

auto main(int argc, char const *argv[]) -> int
{
    std::cout << "entering `main`" << std::endl;
    call_skeleton();
    return 0;
}
