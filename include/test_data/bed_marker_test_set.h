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
const float bmt_marker_mean[3] = {1.4, 0.9, 1.2};
const float bmt_marker_std[3] = {0.66332496, 0.83066239, 0.6};

// Column-major
const float bmt_phen_vals[20] = {0.72031609,  0.59822862, 0.47614114,  -0.62264611, -1.110996,
                                 -1.23308348, 2.06327829, 0.59822862,  -0.86682106, -0.62264611,
                                 0.71351355,  0.2937997,  -0.12591416, -0.82543724, 0.15389508,
                                 -0.26581877, 2.25246434, 0.2937997,   -1.66486495, -0.82543724};
