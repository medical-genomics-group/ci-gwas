#pragma once

#include <vector>

const int NUMBER_OF_LEVELS = 14;

double std_normal_qnorm(const double p);

std::vector<double> threshold_array(const int n, const double alpha);
