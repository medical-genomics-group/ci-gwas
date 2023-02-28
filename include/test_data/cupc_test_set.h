#pragma once

#define CUPCT1_M 5
#define CUPCT1_P 3
#define CUPCT1_W 2

/**
 * This is a row-major matrix with the following structure:
 * m + p rows, w + p cols,
 * {
 * 0.47, 0.39, 0.74, 0.53, 0.62,
 * 0.64, 0.51, 0.35, 0.4,  0.15,
 * 0.80, 0.15, 0.98, 0.21, 0.94,
 * 0.50, 0.00, 0.33, 0.49, 0.25, // here only one forward entry
 * 0.00, 0.00, 0.28, 0.38, 0.24, // here no forward entries
 * 0.00, 0.00, 0.00, 0.21, 0.06, // the last three lines are the phenotype corrs in the lower right corner
 * 0.00, 0.00, 0.21, 0.00, 0.99,
 * 0.00, 0.00, 0.06, 0.99, 0.00
 * }
 */
const float cupct1_c[40] = {0.47, 0.39, 0.74, 0.53, 0.62, 0.64, 0.51, 0.35, 0.4, 0.15, 0.8,
                            0.15, 0.98, 0.21, 0.94, 0.5, 0., 0.33, 0.49, 0.25, 0., 0.,
                            0.28, 0.38, 0.24, 0., 0., 0., 0.21, 0.06, 0., 0., 0.21,
                            0., 0.99, 0., 0., 0.06, 0.99, 0.};

/**
 * This is a row-major compound matrix that consists of two parts.
 * First part:
 * m rows, 2 * w + p cols
 * {
 * 0, 0, 0, 0, 1, 0, 0,
 * 0, 0, 1, 0, 0, 0, 0,
 * 0, 1, 1, 0, 1, 0, 1,
 * 0, 1, 0, 0, 0, 0, 0,
 * 0, 0, 0, 0, 0, 0, 0,
 * }
 * Second part:
 * p rows, m + p cols
 * {
 * 1, 0, 1, 0, 0, 0, 0, 0,
 * 0, 0, 0, 0, 0, 0, 0, 1,
 * 0, 0, 1, 0, 0, 0, 1, 0
 * }
 */
const int cupct1_g[59] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
                          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};
