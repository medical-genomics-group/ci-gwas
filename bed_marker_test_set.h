#pragma once

#define BMT_NUM_INDIVIDUALS 10
#define BMT_NUM_MARKERS 3

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
const unsigned char bmt_a[9] = {0x2a, 0xb0, 0x00, 0xcb, 0xb2, 0x0c, 0x2a, 0xba, 0x00};