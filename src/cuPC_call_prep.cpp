#include <boost/math/distributions/normal.hpp>

const int NUMBER_OF_LEVELS = 14;

double std_normal_qnorm(const double p)
{
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, p);
}

std::array<double, NUMBER_OF_LEVELS> threshold_array(const int n, const double alpha)
{
    std::array<double, NUMBER_OF_LEVELS> thr{0.0};
    // my loop range is exclusive of the last i, so I don't subtract one here
    const int n_thr = (NUMBER_OF_LEVELS < (n - 3)) ? NUMBER_OF_LEVELS : n;
    const double half = 0.5;
    for (size_t i = 0; i < n_thr; i++)
    {
        thr[i] = abs(std_normal_qnorm(half * alpha) / sqrt(n - i - 3));
    }
    return thr;
}