#pragma once

#include <vector>

const int NUMBER_OF_LEVELS = 14;

float std_normal_qnorm(const float p);

std::vector<float> threshold_array(const int n, const float alpha);
