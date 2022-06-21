#pragma once

#define BMT_NUM_INDIVIDUALS 10
#define BMT_NUM_MARKERS 3
#define BMT_NUM_PHEN 2

// That should be this matrix in padded .bed format, column-major,
// without magic numbers
// [
//     [1, 2, 1],
//     [1, 1, 1],
//     [1, 0, 1],
//     [0, 2, 0],
//     [0, 1, 1],
//     [0, 0, 1],
//     [2, 2, 2],
//     [1, 1, 1],
//     [0, 0, 0],
//     [0, 2, 0],
// ]
// Column-major
const unsigned char bmt_marker_vals[9] = {0x2a, 0xb0, 0x00, 0xcb, 0xb2, 0x0c, 0x2a, 0xba, 0x00};
const float bmt_marker_mean[3] = {0.6, 1.1, 0.8};
const float bmt_marker_std[3] = {0.66332496, 0.83066239, 0.6};

// Column-major
const float bmt_phen_vals[20] = {1.18972481,  1.08972481, 0.98972481, 0.08972481,  -0.31027519,
                                 -0.41027519, 2.28972481, 1.08972481, -0.11027519, 0.08972481,
                                 1.89522803,  1.59522803, 1.29522803, 0.79522803,  1.49522803,
                                 1.19522803,  2.99522803, 1.59522803, 0.19522803,  0.79522803};
