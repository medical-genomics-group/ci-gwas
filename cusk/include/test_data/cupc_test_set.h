#pragma once

#include <vector>

const float ALPHA_N10 = 0.00001;
const int N_N10 = 10;
const int SAMPLE_SIZE_N10 = 10000;

const std::vector<float> A_N10{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0,
                               0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                               0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
                               0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0,
                               0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0};

std::vector<float> C_N10 = {
    1.00000000e+00,  -1.90712503e-02, -2.83009343e-02, -3.85512912e-02, -2.87644192e-02,
    6.51970726e-02,  1.92298400e-02,  2.06432873e-02,  -2.06871270e-02, 3.15843695e-03,
    -1.90712503e-02, 1.00000000e+00,  4.91443578e-02,  1.48816325e-02,  2.23954368e-02,
    6.48272344e-02,  5.73356420e-02,  6.73839902e-02,  -1.01657510e-02, 4.09517380e-02,
    -2.83009343e-02, 4.91443578e-02,  1.00000000e+00,  6.57071811e-03,  1.01099579e-02,
    2.91768748e-02,  -8.96051233e-03, 3.85163844e-02,  -8.24654502e-04, -1.77101023e-02,
    -3.85512912e-02, 1.48816325e-02,  6.57071811e-03,  1.00000000e+00,  -2.30972490e-02,
    1.44911049e-02,  1.41387533e-02,  -4.33176707e-04, 6.46368709e-02,  8.67405203e-03,
    -2.87644192e-02, 2.23954368e-02,  1.01099579e-02,  -2.30972490e-02, 1.00000000e+00,
    -1.14501869e-02, -4.99909018e-02, 4.87836476e-02,  -3.98139887e-02, 5.78182720e-02,
    6.51970726e-02,  6.48272344e-02,  2.91768748e-02,  1.44911049e-02,  -1.14501869e-02,
    1.00000000e+00,  1.46948556e-02,  5.81364085e-02,  3.46636236e-02,  -4.00708409e-02,
    1.92298400e-02,  5.73356420e-02,  -8.96051233e-03, 1.41387533e-02,  -4.99909018e-02,
    1.46948556e-02,  1.00000000e+00,  1.43198822e-02,  -1.85019453e-02, 5.23717711e-02,
    2.06432873e-02,  6.73839902e-02,  3.85163844e-02,  -4.33176707e-04, 4.87836476e-02,
    5.81364085e-02,  1.43198822e-02,  1.00000000e+00,  -3.20720779e-02, 3.82462770e-03,
    -2.06871270e-02, -1.01657510e-02, -8.24654502e-04, 6.46368709e-02,  -3.98139887e-02,
    3.46636236e-02,  -1.85019453e-02, -3.20720779e-02, 1.00000000e+00,  -3.38017912e-02,
    3.15843695e-03,  4.09517380e-02,  -1.77101023e-02, 8.67405203e-03,  5.78182720e-02,
    -4.00708409e-02, 5.23717711e-02,  3.82462770e-03,  -3.38017912e-02, 1.00000000e+00};