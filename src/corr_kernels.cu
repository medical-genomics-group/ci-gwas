#include <math.h>
#include <mps/bed_lut.h>
#include <mps/corr_host.h>
#include <mps/corr_kernels.h>
#include <mps/gpuerrors.h>

// Compute row and column indices from the linear index into
// upper triangular matrix.
__device__ void row_col_ix_from_linear_ix(const size_t lin_ix,
                                          const size_t num_rows,
                                          size_t *row_ix,
                                          size_t *col_ix)
{
    float l = num_rows - 1;
    float b = 2 * l - 1;
    float c = 2 * (l - (float)lin_ix);
    float row = std::floor((-b + sqrt(b * b + 4 * c)) / -2.0) + 1.0;
    float h = (-(row * row) + row * (2.0 * l + 1.0)) / 2.0;
    // offset of 1 because we do not include the diagonal
    float col = (float)lin_ix - h + row + 1;
    *row_ix = (size_t)row;
    *col_ix = (size_t)col;
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
    float *results)
{
    size_t tix = threadIdx.x;
    size_t lin_ix = blockIdx.x;
    size_t mv_ix = lin_ix / num_phen;
    size_t phen_ix = lin_ix - (num_phen * mv_ix);
    size_t mv_start_ix = mv_ix * col_len_bytes;
    size_t phen_start_ix = phen_ix * num_individuals;

    float thread_sum_mv_phen = 0.0;
    float thread_sum_phen = 0.0;
    __shared__ float thread_sums_mv_phen[NUMTHREADS];
    __shared__ float thread_sums_phen[NUMTHREADS];
    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t curr_mv_byte_ix = 4 * (size_t)(marker_vals[mv_start_ix + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_val = gpu_bed_lut_a[(curr_mv_byte_ix + j)];
            float phen_val = phen_vals[(phen_start_ix + (4 * i) + j)];

            thread_sum_mv_phen += (mv_val * phen_val);
            thread_sum_phen += phen_val;
        }
    }

    thread_sums_mv_phen[tix] = thread_sum_mv_phen;
    thread_sums_phen[tix] = thread_sum_phen;

    __syncthreads();
    if (tix == 0)
    {
        float s_mv_phen = 0.0;
        float s_phen = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mv_phen += thread_sums_mv_phen[i];
            s_phen += thread_sums_phen[i];
        }

        results[lin_ix] = (s_mv_phen - marker_mean[mv_ix] * s_phen) /
                          ((float)(num_individuals)*marker_std[mv_ix]);
    }
}

// Compute Pearson's r between a pair of standardized phenotype vectors.
__global__ void phen_corr_pearson(
    const float *phen_vals,
    const size_t num_individuals,
    const size_t num_phen,
    float *results)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_phen, &row, &col);
    size_t col_start_a = row * num_individuals;
    size_t col_start_b = col * num_individuals;

    float thread_sum = 0.0;
    __shared__ float thread_sums[NUMTHREADS];
    for (size_t i = tix; i < num_individuals; i += NUMTHREADS)
    {
        float val_a = phen_vals[(col_start_a + i)];
        float val_b = phen_vals[(col_start_b + i)];
        thread_sum += val_a * val_b;
    }

    thread_sums[tix] = thread_sum;

    __syncthreads();
    if (tix == 0)
    {
        float s = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s += thread_sums[i];
        }
        results[blockIdx.x] = s / (float)(num_individuals);
    }
}

// Compute Pearson's r between a pair of marker vectors.
__global__ void bed_marker_corr_pearson(const unsigned char *marker_vals,
                                        const size_t num_markers,
                                        const size_t num_individuals,
                                        const size_t col_len_bytes,
                                        const float *marker_mean,
                                        const float *marker_std,
                                        float *results)
{
    size_t tix = threadIdx.x;
    size_t row;
    size_t col;
    row_col_ix_from_linear_ix(blockIdx.x, num_markers, &row, &col);
    size_t col_start_a = row * col_len_bytes;
    size_t col_start_b = col * col_len_bytes;

    float thread_sum_mvp = 0.0;
    __shared__ float thread_sums_mvp[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t mv_byte_ix_a = 4 * (size_t)(marker_vals[col_start_a + i]);
        size_t mv_byte_ix_b = 4 * (size_t)(marker_vals[col_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float mv_a = gpu_bed_lut_a[(mv_byte_ix_a + j)];
            float mv_b = gpu_bed_lut_a[(mv_byte_ix_b + j)];
            thread_sum_mvp += mv_a * mv_b;
        }
    }

    thread_sums_mvp[tix] = thread_sum_mvp;

    __syncthreads();
    // printf("[tix: %lu]: syncing... \n", tix);
    if (tix == 0)
    {
        float s_mvps = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mvps += thread_sums_mvp[i];
        }

        float num = (s_mvps / (float)num_individuals) - (marker_mean[row] * marker_mean[col]);
        float denom = marker_std[row] * marker_std[col];

        results[blockIdx.x] = num / denom;
    }
}

// A O(n) runtime Kendall implementation for compressed genomic marker data.
// The compression format is expected to be col-major .bed without NaN
// and without leading magic numbers.
__global__ void bed_marker_corr_kendall_npn(
    const unsigned char *marker_vals,
    const size_t num_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results)
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
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += 1.f;
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
        float p = ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
                   (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q = ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
                   (s[5] * (s[6] + s[7])));
        float t = ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
                   (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u = ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
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
    float *results)
{
    // figure out which corr we are computing
    // we have a 2d grid
    size_t tix = threadIdx.x;
    size_t col = blockIdx.x;
    size_t row = blockIdx.y;
    size_t data_start_b = col * col_len_bytes;
    size_t data_start_a = row * col_len_bytes;

    float thread_sum_mvp = 0.0;
    __shared__ float thread_sums_mvp[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        // look up table indices
        size_t aix = 4 * (size_t)(row_marker_vals[data_start_a + i]);
        size_t bix = 4 * (size_t)(col_marker_vals[data_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            thread_sum_mvp += mv_a * mv_b;
        }
    }

    thread_sums_mvp[tix] = thread_sum_mvp;

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        float s_mvps = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mvps += thread_sums_mvp[i];
        }

        float num = (s_mvps / (float)num_individuals) - (marker_mean_row[row] * marker_mean_col[col]);
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
    float *results)
{
    // 1d grid
    size_t tix = threadIdx.x;
    size_t row, col;
    row_col_ix_from_linear_ix(blockIdx.x, num_rows, &row, &col);
    size_t data_start_a = col * col_len_bytes;
    size_t data_start_b = row * col_len_bytes;

    float thread_sum_mvp = 0.0;
    __shared__ float thread_sums_mvp[NUMTHREADS];

    for (size_t i = tix; i < col_len_bytes; i += NUMTHREADS)
    {
        size_t aix = 4 * (size_t)(row_marker_vals[data_start_a + i]);
        size_t bix = 4 * (size_t)(row_marker_vals[data_start_b + i]);
        for (size_t j = 0; (j < 4) && (i * 4 + j < num_individuals); j++)
        {
            float val_a = gpu_bed_lut_a[(aix + j)];
            float val_b = gpu_bed_lut_a[(bix + j)];
            thread_sum_mvp += mv_a * mv_b;
        }
    }

    thread_sums_mvp[tix] = thread_sum_mvp;

    // consolidate thread_sums
    __syncthreads();
    if (tix == 0)
    {
        float s_mvps = 0.0;
        for (size_t i = 0; i < NUMTHREADS; i++)
        {
            s_mvps += thread_sums_mvp[i];
        }

        float num = (s_mvps / (float)num_individuals) - (marker_mean[row] * marker_mean[col]);
        float denom = marker_std[row] * marker_std[col];

        results[blockIdx.x] = num / denom;
    }
}

__global__ void bed_marker_corr_kendall_npn_batched(
    const unsigned char *row_marker_vals,
    const unsigned char *col_marker_vals,
    const size_t num_col_markers,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results)
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
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += 1.f;
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
        float p = ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
                   (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q = ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
                   (s[5] * (s[6] + s[7])));
        float t = ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
                   (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u = ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
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
__global__ void bed_marker_corr_kendall_npn_batched_row(
    const unsigned char *row_marker_vals,
    const size_t num_rows,
    const size_t num_individuals,
    const size_t col_len_bytes,
    float *results)
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
            size_t comp_ix = (size_t)((3 * val_a + val_b));
            thread_sum[comp_ix] += 1.f;
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
        float p = ((s[0] * (s[4] + s[5] + s[7] + s[8])) + (s[1] * (s[5] + s[8])) +
                   (s[3] * (s[7] + s[8])) + (s[4] * s[8]));
        float q = ((s[1] * (s[3] + s[6])) + (s[2] * (s[3] + s[4] + s[6] + s[7])) + (s[4] * s[6]) +
                   (s[5] * (s[6] + s[7])));
        float t = ((s[0] * (s[1] + s[2])) + (s[1] * s[2]) + (s[3] * (s[4] + s[5])) + (s[4] * s[5]) +
                   (s[6] * (s[7] + s[8])) + (s[7] * s[8]));
        float u = ((s[0] * (s[3] + s[6])) + (s[1] * (s[4] + s[7])) + (s[2] * (s[5] + s[8])) +
                   (s[3] * s[6]) + (s[4] * s[7]) + (s[5] * s[8]));

        float kendall_corr = (p - q) / sqrt((p + q + t) * (p + q + u));

        // TODO:: are these indices correct?
        results[blockIdx.x] = sin(M_PI / 2 * kendall_corr);
    }
}