#include <math.h>
#include <mps/bed_lut_gpu.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/gpuerrors.h>
#include <stdint.h>

// this is from https://forums.developer.nvidia.com/t/integer-square-root/198642/2
__device__ unsigned long long int umul_wide(unsigned int a, unsigned int b)
{
    unsigned long long int r;
    asm("mul.wide.u32 %0,%1,%2;\n\t" : "=l"(r) : "r"(a), "r"(b));
    return r;
}

// this is from https://forums.developer.nvidia.com/t/integer-square-root/198642/2
__device__ uint32_t isqrtll(uint64_t a)
{
    uint64_t rem, arg;
    uint32_t b, r, s, scal;

    arg = a;
    /* Normalize argument */
    scal = __clzll(a) & ~1;
    a = a << scal;
    b = a >> 32;
    /* Approximate rsqrt accurately. Make sure it's an underestimate! */
    float fb, fr;
    fb = (float)b;
    asm("rsqrt.approx.ftz.f32 %0,%1; \n\t" : "=f"(fr) : "f"(fb));
    r = (uint32_t)fmaf(1.407374884e14f, fr, -438.0f);
    /* Compute sqrt(a) as a * rsqrt(a) */
    s = __umulhi(r, b);
    /* NR iteration combined with back multiply */
    s = s * 2;
    rem = a - umul_wide(s, s);
    r = __umulhi((uint32_t)(rem >> 32) + 1, r);
    s = s + r;
    /* Denormalize result */
    s = s >> (scal >> 1);
    /* Make sure we get the floor correct; can be off by one to either side */
    rem = arg - umul_wide(s, s);
    if ((int64_t)rem < 0)
        s--;
    else if (rem >= ((uint64_t)s * 2 + 1))
        s++;
    return (arg == 0) ? 0 : s;
}

// Compute row and column indices from the linear index into
// upper triangular matrix.
// __device__ void row_col_ix_from_linear_ix(
//     const size_t lin_ix, const size_t num_rows, size_t *row_ix, size_t *col_ix
// )
// {
//     float l = num_rows - 1;
//     float b = 2 * l - 1;
//     float c = 2 * (l - (float)lin_ix);
//     float row = std::floor((-b + sqrt(b * b + 4 * c)) / -2.0) + 1.0;
//     float h = (-(row * row) + row * (2.0 * l + 1.0)) / 2.0;
//     // offset of 1 because we do not include the diagonal
//     float col = (float)lin_ix - h + row + 1;
//     *row_ix = (size_t)row;
//     *col_ix = (size_t)col;
// }

// Compute row and column indices from the linear index into
// upper triangular matrix.
__device__ void row_col_ix_from_linear_ix(
    const size_t lin_ix, const size_t num_rows, size_t *row_ix, size_t *col_ix
)
{
    uint32_t M = num_rows - 1;
    uint64_t num_ix = (uint64_t)M * ((uint64_t)M + 1) / 2;
    uint64_t ii = num_ix - 1 - lin_ix;
    uint32_t K = (isqrtll(8 * ii + 1) - 1) / 2;
    uint32_t rix = M - 1 - K;
    *row_ix = (size_t)rix;
    *col_ix = (size_t)(rix + 1 + lin_ix - num_ix + (K + 1) * (K + 2) / 2);
}

// Compute row and column indices from a linear index into
// a upper triangular matrix.
// This kernel exists only for testing purposes of the device function.
__global__ void ix_from_linear(
    const size_t lin_ix, const size_t num_rows, size_t *row_ix, size_t *col_ix
)
{
    row_col_ix_from_linear_ix(lin_ix, num_rows, row_ix, col_ix);
}

__global__ void bed_marker_phen_corr_pearson(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
)
{
    size_t tix = threadIdx.x;
    size_t lin_ix = blockIdx.x;
    size_t mv_ix = lin_ix / num_phen;
    size_t phen_ix = lin_ix - (num_phen * mv_ix);
    size_t mv_start_ix = mv_ix * col_len_bytes;
    size_t phen_start_ix = phen_ix * num_individuals;

    float thread_sum_mv_phen = 0.0;
    float thread_sum_phen = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums_mv_phen[NUMTHREADS];
    __shared__ float sums_phen[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t curr_mv_byte_ix = 4 * (size_t)(marker_vals[mv_start_ix + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_valid = gpu_bed_lut_b[(curr_mv_byte_ix + j)];
            float mv_val = gpu_bed_lut_a[(curr_mv_byte_ix + j)];
            float phen_val = phen_vals[(phen_start_ix + (4 * i) + j)];
            if (!isnan(phen_val))
            {
                thread_sum_mv_phen += mv_valid * mv_val * phen_val;
                thread_sum_phen += mv_valid * phen_val;
                thread_num_valid += mv_valid;
            }
        }
    }

    sums_mv_phen[tix] = thread_sum_mv_phen;
    sums_phen[tix] = thread_sum_phen;
    nums_valid[tix] = thread_num_valid;

    __syncthreads();
    if (tix == 0)
    {
        float s_mv_phen = 0.0;
        float s_phen = 0.0;
        float s_valid = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mv_phen += sums_mv_phen[i];
            s_phen += sums_phen[i];
            s_valid += nums_valid[i];
        }

        results[lin_ix] = (s_mv_phen - marker_mean[mv_ix] * s_phen) / (s_valid * marker_std[mv_ix]);
    }
}

__global__ void bed_marker_phen_corr_pearson_scan(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t num_phen,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
)
{
    size_t lin_ix = blockIdx.x;
    size_t mv_ix = lin_ix / num_phen;
    size_t phen_ix = lin_ix - (num_phen * mv_ix);
    size_t mv_start_ix = mv_ix * col_len_bytes;
    size_t phen_start_ix = phen_ix * num_individuals;

    float thread_sum_mv_phen = 0.0;
    float thread_sum_phen = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums_mv_phen[NUMTHREADS];
    __shared__ float sums_phen[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tx; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t curr_mv_byte_ix = 4 * (size_t)(marker_vals[mv_start_ix + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_valid = gpu_bed_lut_b[(curr_mv_byte_ix + j)];
            float mv_val = gpu_bed_lut_a[(curr_mv_byte_ix + j)];
            float phen_val = phen_vals[(phen_start_ix + (4 * i) + j)];
            if (!isnan(phen_val))
            {
                thread_sum_mv_phen += mv_valid * mv_val * phen_val;
                thread_sum_phen += mv_valid * phen_val;
                thread_num_valid += mv_valid;
            }
        }
    }

    sums_mv_phen[tx] = thread_sum_mv_phen;
    sums_phen[tx] = thread_sum_phen;
    nums_valid[tx] = thread_num_valid;

    __syncthreads();
    // consolidate thread_sums
    float tmp_mv_phen = 0.0;
    float tmp_phen = 0.0;
    float tmp_n = 0.0;
    for (int step = 1; step < NUMTHREADS; step = step * 2)
    {
        if (tx < step)
        {
            tmp_mv_phen = sums_mv_phen[tx];
            tmp_phen = sums_phen[tx];
            tmp_n = nums_valid[tx];
        }
        else
        {
            tmp_mv_phen = sums_mv_phen[tx] + sums_mv_phen[tx - step];
            tmp_phen = sums_phen[tx] + sums_phen[tx - step];
            tmp_n = nums_valid[tx] + nums_valid[tx - step];
        }
        __syncthreads();
        sums_mv_phen[tx] = tmp_mv_phen;
        sums_phen[tx] = tmp_phen;
        nums_valid[tx] = tmp_n;
        __syncthreads();
    }

    if (tx == 0)
    {
        float s_mv_phen = sums_mv_phen[NUMTHREADS - 1];
        float s_phen = sums_phen[NUMTHREADS - 1];
        float n = nums_valid[NUMTHREADS - 1];

        results[lin_ix] = (s_mv_phen - marker_mean[mv_ix] * s_phen) / (n * marker_std[mv_ix]);
    }
}

// Compute Pearson's r between a pair of standardized phenotype vectors.
__global__ void phen_corr_pearson(
    const float *phen_vals, const size_t num_individuals, const size_t num_phen, float *results
)
{
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(bx, num_phen, &row, &col);
    size_t col_start_a = row * num_individuals;
    size_t col_start_b = col * num_individuals;

    float thread_sum = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];
    for (size_t i = tx; i < num_individuals; i += NUMTHREADS)
    {
        float val_a = phen_vals[(col_start_a + i)];
        float val_b = phen_vals[(col_start_b + i)];
        if (!(isnan(val_a) || isnan(val_b)))
        {
            thread_sum += val_a * val_b;
            thread_num_valid++;
        }
    }

    sums[tx] = thread_sum;
    nums_valid[tx] = thread_num_valid;

    __syncthreads();
    if (tx == 0)
    {
        float s = 0.0;
        float n = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s += sums[i];
            n += nums_valid[i];
        }
        results[bx] = s / n;
    }
}

// Compute Pearson's r between a pair of standardized phenotype vectors.
__global__ void phen_corr_pearson_scan(
    const float *phen_vals, const size_t num_individuals, const size_t num_phen, float *results
)
{
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(bx, num_phen, &row, &col);
    size_t col_start_a = row * num_individuals;
    size_t col_start_b = col * num_individuals;

    float thread_sum = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tx; i < num_individuals; i += NUMTHREADS)
    {
        float val_a = phen_vals[(col_start_a + i)];
        float val_b = phen_vals[(col_start_b + i)];
        if (!(isnan(val_a) || isnan(val_b)))
        {
            thread_sum += val_a * val_b;
            thread_num_valid++;
        }
    }

    sums[tx] = thread_sum;
    nums_valid[tx] = thread_num_valid;

    // consolidate thread sums
    float tmp_s = 0.0;
    float tmp_n = 0.0;
    __syncthreads();
    for (int step = 1; step < NUMTHREADS; step = step * 2)
    {
        if (tx < step)
        {
            tmp_s = sums[tx];
            tmp_n = nums_valid[tx];
        }
        else
        {
            tmp_s = sums[tx] + sums[tx - step];
            tmp_n = nums_valid[tx] + nums_valid[tx - step];
        }
        __syncthreads();
        sums[tx] = tmp_s;
        nums_valid[tx] = tmp_n;
        __syncthreads();
    }

    if (tx == 0)
    {
        float s = sums[NUMTHREADS - 1];
        float n = nums_valid[NUMTHREADS - 1];
        results[bx] = s / n;
    }
}

// Compute Pearson's r between a pair of marker vectors.
__global__ void bed_marker_corr_pearson(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_markers, &row, &col);
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum_mvp = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums_mvp[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t mv_byte_ix_a = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t mv_byte_ix_b = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_a = gpu_bed_lut_a[(mv_byte_ix_a + j)];
            float mv_b = gpu_bed_lut_a[(mv_byte_ix_b + j)];
            float valid_op = gpu_bed_lut_b[(mv_byte_ix_a + j)] * gpu_bed_lut_b[(mv_byte_ix_b + j)];
            thread_sum_mvp += valid_op * mv_a * mv_b;
            thread_num_valid += valid_op;
        }
    }

    sums_mvp[tix] = thread_sum_mvp;
    nums_valid[tix] = thread_num_valid;

    __syncthreads();
    // printf("[tix: %lu]: syncing... \n", tix);
    if (tix == 0)
    {
        float s_mvps = 0.0;
        float n = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mvps += sums_mvp[i];
            n += nums_valid[i];
        }

        float num = (s_mvps / n) - (marker_mean[row] * marker_mean[col]);
        float denom = marker_std[row] * marker_std[col];

        results[blockIdx.x] = num / denom;
    }
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
__global__ void bed_marker_corr_pearson_npn(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_markers, &row, &col);
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum[9] = {0.0};

    __shared__ float thread_sums[NUMTHREADS][9];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t bix = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += valid_op * 1.f;
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        // produce single sum
        float s[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            for (size_t j = 0; j < 9; j++)
            {
                s[j] += thread_sums[i][j];
            }
        }
        float p =
            ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
             (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q =
            ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
             (s[5] * (s[6] + s[7])));
        float t =
            ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
             (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u =
            ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
             (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        results[blockIdx.x] = sin(M_PI / 2 * kendall_corr);
    }
}

__global__ void bed_marker_corr_pearson_npn_scan(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_markers, &row, &col);
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum[9] = {0.0};

    __shared__ float thread_sums[NUMTHREADS][9];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t bix = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += valid_op * 1.f;
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    // consolidate thread_sums
    float tmp[9] = {0.0};
    __syncthreads();
    for (int step = 1; step < NUMTHREADS; step = step * 2)
    {
        if (tx < step)
        {
            for (size_t i = 0; i < 9; ++i)
            {
                tmp[i] = thread_sums[tx][i];
            }
        }
        else
        {
            for (size_t i = 0; i < 9; ++i)
            {
                tmp[i] = thread_sums[tx][i] + thread_sums[tx - step][i];
            }
        }
        __syncthreads();
        for (size_t i = 0; i < 9; ++i)
        {
            thread_sums[tx][i] = tmp[i];
        }
        __syncthreads();
    }

    if (tix == 0)
    {
        // produce single sum
        float *s = thread_sums[NUMTHREADS - 1];
        float p =
            ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
             (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q =
            ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
             (s[5] * (s[6] + s[7])));
        float t =
            ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
             (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u =
            ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
             (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        results[blockIdx.x] = sin(M_PI / 2 * kendall_corr);
    }
}

// A O(n) runtime pearson implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
__global__ void bed_marker_corr_pearson_batched(
    const unsigned char *row_marker_vals,
    const unsigned char *col_marker_vals,
    const size_t num_col_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const float *marker_mean_row,
    const float *marker_std_row,
    const float *marker_mean_col,
    const float *marker_std_col,
    float *results
)
{
    // figure out which corr we are computing
    // we have a 2d grid
    size_t tix = threadIdx.x;
    size_t col = blockIdx.x;
    size_t row = blockIdx.y;
    size_t data_start_b = col * col_len_bytes;
    size_t data_start_a = row * col_len_bytes;

    float thread_sum_mvp = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums_mvp[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        // look up table indices
        size_t aix = 4 * (size_t)(row_marker_vals[data_start_a + i]);
        size_t bix = 4 * (size_t)(col_marker_vals[data_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            thread_sum_mvp += valid_op * val_a * val_b;
            thread_num_valid += valid_op;
        }
    }

    sums_mvp[tix] = thread_sum_mvp;
    nums_valid[tix] = thread_num_valid;

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        float s_mvps = 0.0;
        float n = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mvps += sums_mvp[i];
            n += nums_valid[i];
        }

        float num = (s_mvps / n) - (marker_mean_row[row] * marker_mean_col[col]);
        float denom = marker_std_row[row] * marker_std_col[col];

        results[num_col_markers * row + col] = num / denom;
    }
}

// A O(n) runtime pearson implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
// These are pearson correlations between the markers that are held constant
// in a batch row.
__global__ void bed_marker_corr_pearson_batched_row(
    const unsigned char *row_marker_vals,
    const size_t num_rows,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    float *results
)
{
    // 1d grid
    size_t tix = threadIdx.x;
    size_t row, col;
    row_col_ix_from_linear_ix(blockIdx.x, num_rows, &row, &col);
    size_t data_start_a = col * col_len_bytes;
    size_t data_start_b = row * col_len_bytes;

    float thread_sum_mvp = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float nums_valid[NUMTHREADS];
    __shared__ float sums_mvp[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(row_marker_vals[data_start_a + i]);
        size_t bix = 4 * (size_t)(row_marker_vals[data_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            thread_sum_mvp += valid_op * val_a * val_b;
            thread_num_valid += valid_op;
        }
    }

    sums_mvp[tix] = thread_sum_mvp;
    nums_valid[tix] = thread_num_valid;

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        float s_mvps = 0.0;
        float n = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mvps += sums_mvp[i];
            n += nums_valid[i];
        }

        float num = s_mvps / n - marker_mean[row] * marker_mean[col];
        float denom = marker_std[row] * marker_std[col];

        results[blockIdx.x] = num / denom;
    }
}

__global__ void bed_marker_corr_pearson_npn_batched(
    const unsigned char *row_marker_vals,
    const unsigned char *col_marker_vals,
    const size_t num_col_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
)
{
    // figure out which corr we are computing
    // we have a 2d grid
    size_t tix = threadIdx.x;
    size_t col = blockIdx.x;
    size_t row = blockIdx.y;
    size_t data_start_b = col * col_len_bytes;
    size_t data_start_a = row * col_len_bytes;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        // look up table indices
        size_t aix = 4 * (size_t)(row_marker_vals[data_start_a + i]);
        size_t bix = 4 * (size_t)(col_marker_vals[data_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += valid_op * 1.f;
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        // produce single sum
        float s[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            for (size_t j = 0; j < 9; j++)
            {
                s[j] += thread_sums[i][j];
            }
        }
        float p =
            ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
             (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q =
            ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
             (s[5] * (s[6] + s[7])));
        float t =
            ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
             (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u =
            ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
             (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        // TODO:: are these indices correct?
        results[num_col_markers * row + col] = sin(M_PI / 2 * kendall_corr);
    }
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
// These are the correlations between the markers that are held constant
// in a batch row.
__global__ void bed_marker_corr_pearson_npn_batched_row(
    const unsigned char *row_marker_vals,
    const size_t num_rows,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results
)
{
    // 1d grid
    size_t tix = threadIdx.x;
    size_t row, col;
    row_col_ix_from_linear_ix(blockIdx.x, num_rows, &row, &col);
    size_t data_start_a = col * col_len_bytes;
    size_t data_start_b = row * col_len_bytes;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(row_marker_vals[data_start_a + i]);
        size_t bix = 4 * (size_t)(row_marker_vals[data_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += valid_op * 1.f;
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        thread_sums[tix][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        // produce single sum
        float s[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            for (size_t j = 0; j < 9; j++)
            {
                s[j] += thread_sums[i][j];
            }
        }
        float p =
            ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
             (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q =
            ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
             (s[5] * (s[6] + s[7])));
        float t =
            ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
             (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u =
            ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
             (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        // TODO:: are these indices correct?
        results[blockIdx.x] = sin(M_PI / 2 * kendall_corr);
    }
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
// TODO: use scan sum to compute this in O(log n) time.
__global__ void bed_marker_corr_pearson_npn_sparse(
    const unsigned char *marker_vals,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const size_t row_in,
    float *results
)
{
    size_t row = row_in;
    size_t col = row + bx + 1;
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    for (size_t i = tx; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t bix = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += valid_op * 1.f;
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        thread_sums[tx][i] = thread_sum[i];
    }

    // consolidate thread_sums
    __syncthreads();
    if (tx == 0)
    {
        // produce single sum
        float s[9] = {0.0};
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            for (size_t j = 0; j < 9; j++)
            {
                s[j] += thread_sums[i][j];
            }
        }
        float p =
            ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
             (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q =
            ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
             (s[5] * (s[6] + s[7])));
        float t =
            ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
             (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u =
            ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
             (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        results[bx] = sin(M_PI / 2 * kendall_corr);
    }
}

// A O((0.24n / t) + log t) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
__global__ void bed_marker_corr_pearson_npn_sparse_scan(
    const unsigned char *marker_vals,
    const size_t num_individuals,
    const size_t col_len_bytes,
    const size_t row_in,
    float *results
)
{
    size_t row = row_in;
    size_t col = row + bx + 1;
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum[9] = {0.0};
    __shared__ float thread_sums[NUMTHREADS][9];

    for (size_t i = tx; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t bix = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            float valid_op = gpu_bed_lut_b[(aix + j)] * gpu_bed_lut_b[(bix + j)];
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += valid_op * 1.f;
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        thread_sums[tx][i] = thread_sum[i];
    }

    // consolidate thread_sums
    float tmp[9] = {0.0};
    __syncthreads();
    for (int step = 1; step < NUMTHREADS; step = step * 2)
    {
        if (tx < step)
        {
            for (size_t i = 0; i < 9; ++i)
            {
                tmp[i] = thread_sums[tx][i];
            }
        }
        else
        {
            for (size_t i = 0; i < 9; ++i)
            {
                tmp[i] = thread_sums[tx][i] + thread_sums[tx - step][i];
            }
        }
        __syncthreads();
        for (size_t i = 0; i < 9; ++i)
        {
            thread_sums[tx][i] = tmp[i];
        }
        __syncthreads();
    }

    if (tx == 0)
    {
        // produce single sum
        float *s = thread_sums[NUMTHREADS - 1];
        float p =
            ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
             (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q =
            ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
             (s[5] * (s[6] + s[7])));
        float t =
            ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
             (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u =
            ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
             (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        results[bx] = sin(M_PI / 2 * kendall_corr);
    }
}

__global__ void bed_marker_phen_corr_pearson_sparse(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_individuals,
    const size_t num_phen,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    const size_t corr_width,
    const size_t row_in,
    float *results
)
{
    size_t row = row_in;
    size_t phen_ix = bx;
    size_t mv_start_ix = row * col_len_bytes;
    size_t phen_start_ix = phen_ix * num_individuals;

    float thread_sum_mv_phen = 0.0;
    float thread_sum_phen = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums_mv_phen[NUMTHREADS];
    __shared__ float sums_phen[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tx; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t curr_mv_byte_ix = 4 * (size_t)(marker_vals[mv_start_ix + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_val = gpu_bed_lut_a[(curr_mv_byte_ix + j)];
            float phen_val = phen_vals[(phen_start_ix + (4 * i) + j)];
            if (!isnan(phen_val))
            {
                float valid_op = gpu_bed_lut_b[(curr_mv_byte_ix + j)];
                thread_sum_mv_phen += valid_op * (mv_val * phen_val);
                thread_sum_phen += valid_op * phen_val;
                thread_num_valid += valid_op;
            }
        }
    }

    sums_mv_phen[tx] = thread_sum_mv_phen;
    sums_phen[tx] = thread_sum_phen;
    nums_valid[tx] = thread_num_valid;

    __syncthreads();
    if (tx == 0)
    {
        float s_mv_phen = 0.0;
        float s_phen = 0.0;
        float n = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mv_phen += sums_mv_phen[i];
            s_phen += sums_phen[i];
            n += nums_valid[i];
        }

        results[corr_width + phen_ix] = (s_mv_phen - *marker_mean * s_phen) / (n * *marker_std);
    }
}

__global__ void bed_marker_phen_corr_pearson_sparse_scan(
    const unsigned char *marker_vals,
    const float *phen_vals,
    const size_t num_individuals,
    const size_t num_phen,
    const size_t col_len_bytes,
    const float *marker_mean,
    const float *marker_std,
    const size_t corr_width,
    const size_t row_in,
    float *results
)
{
    size_t row = row_in;
    size_t phen_ix = bx;
    size_t mv_start_ix = row * col_len_bytes;
    size_t phen_start_ix = phen_ix * num_individuals;

    float thread_sum_mv_phen = 0.0;
    float thread_sum_phen = 0.0;
    float thread_num_valid = 0.0;

    __shared__ float sums_mv_phen[NUMTHREADS];
    __shared__ float sums_phen[NUMTHREADS];
    __shared__ float nums_valid[NUMTHREADS];

    for (size_t i = tx; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t curr_mv_byte_ix = 4 * (size_t)(marker_vals[mv_start_ix + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_val = gpu_bed_lut_a[(curr_mv_byte_ix + j)];
            float phen_val = phen_vals[(phen_start_ix + (4 * i) + j)];
            if (!isnan(phen_val))
            {
                float valid_op = gpu_bed_lut_b[(curr_mv_byte_ix + j)];
                thread_sum_mv_phen += valid_op * (mv_val * phen_val);
                thread_sum_phen += valid_op * phen_val;
                thread_num_valid += valid_op;
            }
        }
    }

    sums_mv_phen[tx] = thread_sum_mv_phen;
    sums_phen[tx] = thread_sum_phen;
    nums_valid[tx] = thread_num_valid;

    __syncthreads();
    // consolidate thread_sums
    float tmp_mv_phen = 0.0;
    float tmp_phen = 0.0;
    float tmp_n = 0.0;
    for (int step = 1; step < NUMTHREADS; step = step * 2)
    {
        if (tx < step)
        {
            tmp_mv_phen = sums_mv_phen[tx];
            tmp_phen = sums_phen[tx];
            tmp_n = nums_valid[tx];
        }
        else
        {
            tmp_mv_phen = sums_mv_phen[tx] + sums_mv_phen[tx - step];
            tmp_phen = sums_phen[tx] + sums_phen[tx - step];
            tmp_n = nums_valid[tx] + nums_valid[tx - step];
        }
        __syncthreads();
        sums_mv_phen[tx] = tmp_mv_phen;
        sums_phen[tx] = tmp_phen;
        nums_valid[tx] = tmp_n;
        __syncthreads();
    }

    if (tx == 0)
    {
        float s_mv_phen = sums_mv_phen[NUMTHREADS - 1];
        float s_phen = sums_phen[NUMTHREADS - 1];
        float n = nums_valid[NUMTHREADS - 1];

        results[corr_width + phen_ix] = (s_mv_phen - *marker_mean * s_phen) / (n * *marker_std);
    }
}
