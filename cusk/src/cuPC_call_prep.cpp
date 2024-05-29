#include <mps/cuPC_call_prep.h>

#include <boost/math/distributions/normal.hpp>

const int MAX_LEVEL = 14;

float std_normal_qnorm(const float p)
{
    boost::math::normal dist(0.0, 1.0);
    return quantile(dist, p);
}

std::vector<float> threshold_array(const int n, const float alpha)
{
    std::vector<float> thr;
    // my loop range is exclusive of the last i, so I don't subtract one here
    const float half = 0.5;
    for (size_t i = 0; i < MAX_LEVEL + 1; i++)
    {
        thr.push_back(abs(std_normal_qnorm(half * alpha)) / sqrt(n - i - 3));
    }
    return thr;
}

float hetcor_threshold(const float alpha)
{
    return abs(std_normal_qnorm(half * alpha));
}