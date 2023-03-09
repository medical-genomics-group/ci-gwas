#include <mps/cuPC_call_prep.h>

#include <boost/math/distributions/normal.hpp>

double std_normal_qnorm(const double p)
{
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, p);
}

std::vector<double> threshold_array(const int n, const double alpha)
{
    std::vector<double> thr;
    // my loop range is exclusive of the last i, so I don't subtract one here
    const double half = 0.5;
    for (size_t i = 0; i < NUMBER_OF_LEVELS; i++)
    {
        thr.push_back(abs(std_normal_qnorm(half * alpha) / sqrt(n - i - 3)));
    }
    return thr;
}
