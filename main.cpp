#include <array>
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <iostream>

#include "bed_marker_test_set.h"
#include "corr_test_set.h"
#include "cuPC-S.h"
#include "kendall.h"

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
    memset(sepset, 0, sizeof(sepset));
    int l = 0;
    Skeleton(C.data(), &p, G.data(), Th.data(), &l, &max_level, pmax.data(), sepset);
}

auto corr_matrix_size(size_t num_markers, size_t num_individuals) -> size_t
{
    return num_markers * (num_markers - 1) / 2;
}

void call_cu_corr()
{
    const size_t num_markers = TEST_NUM_MARKERS;
    const size_t num_individuals = TEST_NUM_INDIVIDUALS;
    const size_t cm_size = corr_matrix_size(num_markers, num_individuals);
    float marker_corr[cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));
    cu_corr_npn(test_a, num_markers, num_individuals, marker_corr);

    for (size_t i = 0; i < cm_size; i++) {
        std::cout << marker_corr[i] << std::endl;
    }
}

void call_cu_bed_corr()
{
    const size_t num_markers = BMT_NUM_MARKERS;
    const size_t num_individuals = BMT_NUM_INDIVIDUALS;
    const size_t cm_size = corr_matrix_size(num_markers, num_individuals);
    float marker_corr[cm_size];
    memset(marker_corr, 0.0, sizeof(marker_corr));
    cu_bed_corr_npn(bmt_a, num_markers, num_individuals, marker_corr);

    for (size_t i = 0; i < cm_size; i++) {
        std::cout << marker_corr[i] << std::endl;
    }
}

auto main(int argc, char const *argv[]) -> int
{
    call_cu_bed_corr();
    std::cout << "All done!" << std::endl;
    return 0;
}
