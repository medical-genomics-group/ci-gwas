#include <math.h>
#include <mps/cuPC-S.h>
#include <mps/hetcor-cuPC-S.h>
#include <mps/gpuerrors.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>

__device__ void print_sepset(int *var_ixs, int *time_index, int l)
{
    int a = var_ixs[0];
    int b = var_ixs[1];
    if ((time_index[a] > 0) || (time_index[b] > 0)) {
        if (l == 1) {
            printf("(%d, %d) || removed at ti (%d, %d): %d \n", a, b,
            time_index[a], time_index[b], var_ixs[2]);
        } else if (l == 2 ) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d \n", a, b,
            time_index[a], time_index[b], var_ixs[2], var_ixs[3]);
        } else if (l == 3) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d \n", a, b,
            time_index[a], time_index[b], var_ixs[2], var_ixs[3], var_ixs[4]);
        } else if (l == 4) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d \n", a, b,
            time_index[a], time_index[b], var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5]);
        } else if (l == 5) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d \n", a, b,
            time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6]);
        } else if (l == 6) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d \n", a, b,
            time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7]);
        } else if (l == 7) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8]);
        } else if (l == 8) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9]);
        } else if (l == 9) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9], var_ixs[10]);
        } else if (l == 10) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9], var_ixs[10], var_ixs[11]);
        } else if (l == 11) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9], var_ixs[10], var_ixs[11], var_ixs[12]);
        } else if (l == 12) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9], var_ixs[10], var_ixs[11], var_ixs[12], var_ixs[13]);
        } else if (l == 13) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9], var_ixs[10], var_ixs[11], var_ixs[12], var_ixs[13], var_ixs[14]);
        } else if (l == 14) {
            printf("(%d, %d) || removed at ti (%d, %d): %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n", a, b, time_index[a], time_index[b],
            var_ixs[2], var_ixs[3], var_ixs[4], var_ixs[5], var_ixs[6], var_ixs[7], var_ixs[8], var_ixs[9], var_ixs[10], var_ixs[11], var_ixs[12], var_ixs[13], var_ixs[14], var_ixs[15]);
        }
    }
}

/** @brief Computes the skeleton using hetcor correlations and estimates of effective sample sizes.
 *
 * @param[in]  C  Pointer to full, square, correlation matrix
 * @param[in]  P  Pointer to number of variables
 * @param[in]  G  Pointer to adjacency matrix
 * @param[in]  N  Pointer to effective sample size matrix
 * @param[in]  Th  Pointer to alpha / 2 percentile
 * @param[in]  l  Pointer to current level
 * @param[in]  maxlevel  Pointer to maximal level
 * @param[in]  time_index  Pointer to time indices of all variables
 * @return return_name return description
 */
void hetcor_skeleton(
    float *C, int *P, int *G, float *N, float *Th, int *l, const int *maxlevel, const int *time_index
)
{
    float *C_cuda;  // Copy of C array in GPU
    int *G_cuda;  // Copy of G Array in GPU
    float *N_cuda; // Copy of N Array in GPU
    int *nprime_cuda;
    int *GPrime_cuda;
    int *mutex_cuda;
    int *time_index_cuda;

    int n = *P;
    float th = *Th;
    int nprime = 0;
    dim3 BLOCKS_PER_GRID;
    dim3 THREADS_PER_BLOCK;

    bool FinishFlag = false;

    HANDLE_ERROR(cudaMalloc((void **)&mutex_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&nprime_cuda, 1 * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&GPrime_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&C_cuda, n * n * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void **)&G_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&N_cuda, n * n * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void **)&time_index_cuda, n * sizeof(int)));
    // copy correlation matrix from CPU to GPU
    HANDLE_ERROR(cudaMemcpy(C_cuda, C, n * n * sizeof(float), cudaMemcpyHostToDevice));
    // copy effective sample size matrix from CPU to GPU
    HANDLE_ERROR(cudaMemcpy(N_cuda, N, n * n * sizeof(float), cudaMemcpyHostToDevice));
    // copy time indices from host to device
    HANDLE_ERROR(cudaMemcpy(time_index_cuda, time_index, n * sizeof(int), cudaMemcpyHostToDevice));
    // initialize a 0 matrix
    HANDLE_ERROR(cudaMemset(mutex_cuda, 0, n * n * sizeof(int)));

    CudaCheckError();
    //----------------------------------------------------------
    for (*l = 0; *l <= ML && !FinishFlag && *l <= *maxlevel; *l = *l + 1)
    {
        CudaCheckError();
        if (*l == 0)
        {
            printf("Starting lvl 0\n");
            fflush(stdout);

            if ((n * n) < 1024)
            {
                BLOCKS_PER_GRID = dim3(1, 1, 1);
                THREADS_PER_BLOCK = dim3(32, 32, 1);
                cal_Indepl0_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                    C_cuda, G_cuda, N_cuda, n, th
                );
                CudaCheckError();
            }
            else
            {
                BLOCKS_PER_GRID = dim3(ceil(((float)(n)) / 32.0), ceil(((float)(n)) / 32.0), 1);
                THREADS_PER_BLOCK = dim3(32, 32, 1);
                cal_Indepl0_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                    C_cuda, G_cuda, N_cuda, n, th
                );
                CudaCheckError();
            }
            BLOCKS_PER_GRID = dim3(n * n, 1, 1);
            THREADS_PER_BLOCK = dim3(ML, 1, 1);
            CudaCheckError();
        }
        else
        {
            //================================> Start Scan Process <===============================
            HANDLE_ERROR(cudaMemset(nprime_cuda, 0, 1 * sizeof(int)));
            BLOCKS_PER_GRID = dim3(1, n, 1);
            THREADS_PER_BLOCK = dim3(1024, 1, 1);
            scan_compact<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, n * sizeof(int)>>>(
                GPrime_cuda, G_cuda, n, nprime_cuda
            );
            CudaCheckError();
            HANDLE_ERROR(cudaMemcpy(&nprime, nprime_cuda, 1 * sizeof(int), cudaMemcpyDeviceToHost));

            printf("nprime: %i \n", nprime);
            fflush(stdout);

            //================================> Begin The Gaussian CI Test
            //<==============================
            // CHeck whether a CI test is possible
            if (nprime - 1 < *l)
            {  // if not:
                *l = *l - 1;
                FinishFlag = true;
                break;
            }

            if (*l == 1)
            {
                printf("Starting lvl 1\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL1, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL1, 1, 1);
                cal_Indepl1_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 2)
            {
                printf("Starting lvl 2\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL2, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL2, 1, 1);
                cal_Indepl2_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 3)
            {
                printf("Starting lvl 3\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL3, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL3, 1, 1);
                cal_Indepl3_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 4)
            {
                printf("Starting lvl 4\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL4, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL4, 1, 1);
                cal_Indepl4_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 5)
            {
                printf("Starting lvl 5\n");
                fflush(stdout);

                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL5, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL5, 1, 1);
                cal_Indepl5_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 6)
            {
                printf("Starting lvl 6\n");
                fflush(stdout);

                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL6, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL6, 1, 1);
                cal_Indepl6_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
                CudaCheckError();
            }
            else if (*l == 7)
            {
                printf("Starting lvl 7\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL7, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL7, 1, 1);
                cal_Indepl7_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 8)
            {
                printf("Starting lvl 8\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL8, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL8, 1, 1);
                cal_Indepl8_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 9)
            {
                printf("Starting lvl 9\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL9, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL9, 1, 1);
                cal_Indepl9_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 10)
            {
                printf("Starting lvl 10\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL10, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL10, 1, 1);
                cal_Indepl10_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 11)
            {
                printf("Starting lvl 11\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL11, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL11, 1, 1);
                cal_Indepl11_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 12)
            {
                printf("Starting lvl 12\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL12, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL12, 1, 1);
                cal_Indepl12_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 13)
            {
                printf("Starting lvl 13\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL13, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL13, 1, 1);
                cal_Indepl13_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else if (*l == 14)
            {
                printf("Starting lvl 14\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL14, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL14, 1, 1);
                cal_Indepl14_ess<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, N_cuda, GPrime_cuda, mutex_cuda, n, th, time_index_cuda
                );
                CudaCheckError();
            }
            else
            {
                // TODO: add PC serial
            }
        }
    }  // if l > 0

    // Copy Graph G from GPU to CPU
    HANDLE_ERROR(cudaMemcpy(G, G_cuda, n * n * sizeof(int), cudaMemcpyDeviceToHost));
    // Free allocated space
    HANDLE_ERROR(cudaFree(C_cuda));
    HANDLE_ERROR(cudaFree(GPrime_cuda));
    HANDLE_ERROR(cudaFree(G_cuda));
    HANDLE_ERROR(cudaFree(mutex_cuda));
}  // Skeleton

__global__ void cal_Indepl0_ess(float *C, int *G, float *N, int n, float th)
{
    int row = blockDim.x * bx + tx;
    int col = blockDim.y * by + ty;
    if (row < col && col < n)
    {
        float res = C[row * n + col];
        res = abs(0.5 * log(abs((1 + res) / (1 - res))));
        float loc_th = th / sqrt(N[row * n + col] - 3.0);
        if (res < loc_th)
        {
            G[row * n + col] = 0;
            G[col * n + row] = 0;
        }
        else
        {
            G[row * n + col] = 1;
            G[col * n + row] = 1;
        }
    }
    if (row == col && col < n)
    {
        G[row * n + col] = 0;
        G[col * n + row] = 0;
    }
}

__global__ void cal_Indepl1_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 1;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx;
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    __shared__ int NoEdgeFlag;
    float M0;
    float H[2][2];
    float M1[2];
    float rho, Z;
    extern __shared__ int G_Chunk[];
    int var_ixs[3];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if ((SizeOfArr % ParGivenL1) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL1;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL1 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL1) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL1] = GPrime[XIdx * n + tx + cnt * ParGivenL1];
        }
        __syncthreads();
    }

    if ((SizeOfArr % (ParGivenL1 * NumOfBlockForEachNodeL1)) == 0)
    {
        NumOfGivenJump = SizeOfArr / (ParGivenL1 * NumOfBlockForEachNodeL1);
    }
    else
    {
        NumOfGivenJump = SizeOfArr / (ParGivenL1 * NumOfBlockForEachNodeL1) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        __syncthreads();
        NbrIdxPointer = tx + bx * ParGivenL1 + d1 * ParGivenL1 * NumOfBlockForEachNodeL1;
        NoEdgeFlag = 1;
        __syncthreads();
        if (NbrIdxPointer < SizeOfArr)
        {
            NbrIdx = G_Chunk[NbrIdxPointer];
            M1[0] = C[XIdx * n + NbrIdx];
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if (d2 == NbrIdxPointer)
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (time_index[NbrIdx] > max(time_index[XIdx], time_index[YIdx]))
                {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    M1[1] = C[YIdx * n + NbrIdx];
                    
                    H[0][0] = 1 - (M1[0] * M1[0]);
                    H[0][1] = M0 - (M1[0] * M1[1]);
                    H[1][1] = 1 - (M1[1] * M1[1]);

                    rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                    Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx;
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 3, n) - 1.0 - 3.0);

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl2_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 2;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[2];
    int NbrIdx[2];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    float M0;
    float M1[2][2];
    float M2[2][2];
    float M2Inv[2][2];
    float M1MulM2Inv[2][2];
    float H[2][2];
    float rho;
    float Z;
    int var_ixs[4];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 2)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL2) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL2;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL2 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL2) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL2] = GPrime[XIdx * n + tx + cnt * ParGivenL2];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 2, &NumOfComb);
    if ((NumOfComb % (ParGivenL2 * NumOfBlockForEachNodeL2)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL2 * NumOfBlockForEachNodeL2);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL2 * NumOfBlockForEachNodeL2) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if ((tx + bx * ParGivenL2 + d1 * ParGivenL2 * NumOfBlockForEachNodeL2) < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                2,
                tx + bx * ParGivenL2 + d1 * ParGivenL2 * NumOfBlockForEachNodeL2 + 1
            );
            NbrIdx[0] = G_Chunk[NbrIdxPointer[0] - 1];
            NbrIdx[1] = G_Chunk[NbrIdxPointer[1] - 1];
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[0][0] = 1;

            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            pseudoinversel2(M2, M2Inv);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level))
                {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 2; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }
                    H[0][0] = 1 - (M1MulM2Inv[0][0] * M1[0][0] + M1MulM2Inv[0][1] * M1[0][1]);
                    H[0][1] = M0 - (M1MulM2Inv[0][0] * M1[1][0] + M1MulM2Inv[0][1] * M1[1][1]);
                    H[1][1] = 1 - (M1MulM2Inv[1][0] * M1[1][0] + M1MulM2Inv[1][1] * M1[1][1]);

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 4, n) - 2.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = 0.5 * abs(log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl3_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 3;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[3];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[3];
    int var_ixs[5];
    float M0;
    float M1[2][3];
    float M2[3][3];
    float M2Inv[3][3];
    float M1MulM2Inv[2][3];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 3)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL3) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL3;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL3 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL3) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL3] = GPrime[XIdx * n + tx + cnt * ParGivenL3];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 3, &NumOfComb);
    if ((NumOfComb % (ParGivenL3 * NumOfBlockForEachNodeL3)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL3 * NumOfBlockForEachNodeL3);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL3 * NumOfBlockForEachNodeL3) + 1;
    }
    __syncthreads();
    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL3 + d1 * ParGivenL3 * NumOfBlockForEachNodeL3 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                3,
                tx + bx * ParGivenL3 + d1 * ParGivenL3 * NumOfBlockForEachNodeL3 + 1
            );
            NbrIdx[0] = G_Chunk[NbrIdxPointer[0] - 1];
            NbrIdx[1] = G_Chunk[NbrIdxPointer[1] - 1];
            NbrIdx[2] = G_Chunk[NbrIdxPointer[2] - 1];
            M2[0][0] = 1;
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[0][2] = C[NbrIdx[0] * n + NbrIdx[2]];
            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[1][2] = C[NbrIdx[1] * n + NbrIdx[2]];
            M2[2][0] = M2[0][2];
            M2[2][1] = M2[1][2];
            M2[2][2] = 1;

            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][2] = C[XIdx * n + NbrIdx[2]];

            pseudoinversel3(M2, M2Inv);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];

                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    M1[1][2] = C[YIdx * n + NbrIdx[2]];
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 3; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 3; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 3; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 5, n) - 3.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl4_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 4;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[4];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[4];

    float M0;
    float M1[2][4];
    float M2[4][4];
    float M2Inv[4][4];
    float M1MulM2Inv[2][4];
    float H[2][2];
    float rho;
    float Z;
    int var_ixs[6];

    float v[4][4];
    float w[4], rv1[4];
    float res1[4][4];
    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 4)
    {
        return;
    }
    if ((SizeOfArr % ParGivenL4) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL4;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL4 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL4) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL4] = GPrime[XIdx * n + tx + cnt * ParGivenL4];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 4, &NumOfComb);
    if ((NumOfComb % (ParGivenL4 * NumOfBlockForEachNodeL4)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL4 * NumOfBlockForEachNodeL4);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL4 * NumOfBlockForEachNodeL4) + 1;
    }
    __syncthreads();
    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL4 + d1 * ParGivenL4 * NumOfBlockForEachNodeL4 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                4,
                tx + bx * ParGivenL4 + d1 * ParGivenL4 * NumOfBlockForEachNodeL4 + 1
            );
            NbrIdx[0] = G_Chunk[NbrIdxPointer[0] - 1];
            NbrIdx[1] = G_Chunk[NbrIdxPointer[1] - 1];
            NbrIdx[2] = G_Chunk[NbrIdxPointer[2] - 1];
            NbrIdx[3] = G_Chunk[NbrIdxPointer[3] - 1];

            M2[0][0] = 1;
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[0][2] = C[NbrIdx[0] * n + NbrIdx[2]];
            M2[0][3] = C[NbrIdx[0] * n + NbrIdx[3]];

            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[1][2] = C[NbrIdx[1] * n + NbrIdx[2]];
            M2[1][3] = C[NbrIdx[1] * n + NbrIdx[3]];

            M2[2][0] = M2[0][2];
            M2[2][1] = M2[1][2];
            M2[2][2] = 1;
            M2[2][3] = C[NbrIdx[2] * n + NbrIdx[3]];

            M2[3][0] = M2[0][3];
            M2[3][1] = M2[1][3];
            M2[3][2] = M2[2][3];
            M2[3][3] = 1;

            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][2] = C[XIdx * n + NbrIdx[2]];
            M1[0][3] = C[XIdx * n + NbrIdx[3]];
            pseudoinversel4(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];

                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    M1[1][2] = C[YIdx * n + NbrIdx[2]];
                    M1[1][3] = C[YIdx * n + NbrIdx[3]];
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 4; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 4; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 4; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 6, n) - 4.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl5_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 5;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[5];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[5];
    int var_ixs[7];

    float M0;
    float M1[2][5];
    float M2[5][5];
    float M2Inv[5][5];
    float M1MulM2Inv[2][5];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[5][5];
    float w[5], rv1[5];
    float res1[5][5];
    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 5)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL5) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL5;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL5 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL5) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL5] = GPrime[XIdx * n + tx + cnt * ParGivenL5];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 5, &NumOfComb);
    if ((NumOfComb % (ParGivenL5 * NumOfBlockForEachNodeL5)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL5 * NumOfBlockForEachNodeL5);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL5 * NumOfBlockForEachNodeL5) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL5 + d1 * ParGivenL5 * NumOfBlockForEachNodeL5 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                5,
                tx + bx * ParGivenL5 + d1 * ParGivenL5 * NumOfBlockForEachNodeL5 + 1
            );
            for (int tmp = 0; tmp < 5; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            M2[0][0] = 1;
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[0][2] = C[NbrIdx[0] * n + NbrIdx[2]];
            M2[0][3] = C[NbrIdx[0] * n + NbrIdx[3]];
            M2[0][4] = C[NbrIdx[0] * n + NbrIdx[4]];

            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[1][2] = C[NbrIdx[1] * n + NbrIdx[2]];
            M2[1][3] = C[NbrIdx[1] * n + NbrIdx[3]];
            M2[1][4] = C[NbrIdx[1] * n + NbrIdx[4]];

            M2[2][0] = M2[0][2];
            M2[2][1] = M2[1][2];
            M2[2][2] = 1;
            M2[2][3] = C[NbrIdx[2] * n + NbrIdx[3]];
            M2[2][4] = C[NbrIdx[2] * n + NbrIdx[4]];

            M2[3][0] = M2[0][3];
            M2[3][1] = M2[1][3];
            M2[3][2] = M2[2][3];
            M2[3][3] = 1;
            M2[3][4] = C[NbrIdx[3] * n + NbrIdx[4]];

            M2[4][0] = M2[0][4];
            M2[4][1] = M2[1][4];
            M2[4][2] = M2[2][4];
            M2[4][3] = M2[3][4];
            M2[4][4] = 1;

            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][2] = C[XIdx * n + NbrIdx[2]];
            M1[0][3] = C[XIdx * n + NbrIdx[3]];
            M1[0][4] = C[XIdx * n + NbrIdx[4]];
            pseudoinversel5(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation

                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    M1[1][2] = C[YIdx * n + NbrIdx[2]];
                    M1[1][3] = C[YIdx * n + NbrIdx[3]];
                    M1[1][4] = C[YIdx * n + NbrIdx[4]];
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 5; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 5; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 5; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 7, n) - 5.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl6_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 6;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[6];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[6];
    int var_ixs[8];

    float M0;
    float M1[2][6];
    float M2[6][6];
    float M2Inv[6][6];
    float M1MulM2Inv[2][6];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[6][6];
    float w[6], rv1[6];
    float res1[6][6];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 6)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL6) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL6;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL6 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL6) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL6] = GPrime[XIdx * n + tx + cnt * ParGivenL6];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 6, &NumOfComb);
    if ((NumOfComb % (ParGivenL6 * NumOfBlockForEachNodeL6)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL6 * NumOfBlockForEachNodeL6);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL6 * NumOfBlockForEachNodeL6) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL6 + d1 * ParGivenL6 * NumOfBlockForEachNodeL6 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                6,
                tx + bx * ParGivenL6 + d1 * ParGivenL6 * NumOfBlockForEachNodeL6 + 1
            );
            for (int tmp = 0; tmp < 6; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            M2[0][0] = 1;
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[0][2] = C[NbrIdx[0] * n + NbrIdx[2]];
            M2[0][3] = C[NbrIdx[0] * n + NbrIdx[3]];
            M2[0][4] = C[NbrIdx[0] * n + NbrIdx[4]];
            M2[0][5] = C[NbrIdx[0] * n + NbrIdx[5]];

            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[1][2] = C[NbrIdx[1] * n + NbrIdx[2]];
            M2[1][3] = C[NbrIdx[1] * n + NbrIdx[3]];
            M2[1][4] = C[NbrIdx[1] * n + NbrIdx[4]];
            M2[1][5] = C[NbrIdx[1] * n + NbrIdx[5]];

            M2[2][0] = M2[0][2];
            M2[2][1] = M2[1][2];
            M2[2][2] = 1;
            M2[2][3] = C[NbrIdx[2] * n + NbrIdx[3]];
            M2[2][4] = C[NbrIdx[2] * n + NbrIdx[4]];
            M2[2][5] = C[NbrIdx[2] * n + NbrIdx[5]];

            M2[3][0] = M2[0][3];
            M2[3][1] = M2[1][3];
            M2[3][2] = M2[2][3];
            M2[3][3] = 1;
            M2[3][4] = C[NbrIdx[3] * n + NbrIdx[4]];
            M2[3][5] = C[NbrIdx[3] * n + NbrIdx[5]];

            M2[4][0] = M2[0][4];
            M2[4][1] = M2[1][4];
            M2[4][2] = M2[2][4];
            M2[4][3] = M2[3][4];
            M2[4][4] = 1;
            M2[4][5] = C[NbrIdx[4] * n + NbrIdx[5]];

            M2[5][0] = M2[0][5];
            M2[5][1] = M2[1][5];
            M2[5][2] = M2[2][5];
            M2[5][3] = M2[3][5];
            M2[5][4] = M2[4][5];
            M2[5][5] = 1;

            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][2] = C[XIdx * n + NbrIdx[2]];
            M1[0][3] = C[XIdx * n + NbrIdx[3]];
            M1[0][4] = C[XIdx * n + NbrIdx[4]];
            M1[0][5] = C[XIdx * n + NbrIdx[5]];

            pseudoinversel6(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    M1[1][2] = C[YIdx * n + NbrIdx[2]];
                    M1[1][3] = C[YIdx * n + NbrIdx[3]];
                    M1[1][4] = C[YIdx * n + NbrIdx[4]];
                    M1[1][5] = C[YIdx * n + NbrIdx[5]];
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 6; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 6; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 6; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 8, n) - 6.0 - 3.0);
                    
                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl7_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 7;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[7];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[7];
    int var_ixs[9];

    float M0;
    float M1[2][7];
    float M2[7][7];
    float M2Inv[7][7];
    float M1MulM2Inv[2][7];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[7][7];
    float w[7], rv1[7];
    float res1[7][7];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 7)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL7) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL7;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL7 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL7) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL7] = GPrime[XIdx * n + tx + cnt * ParGivenL7];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 7, &NumOfComb);
    if ((NumOfComb % (ParGivenL7 * NumOfBlockForEachNodeL7)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL7 * NumOfBlockForEachNodeL7);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL7 * NumOfBlockForEachNodeL7) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL7 + d1 * ParGivenL7 * NumOfBlockForEachNodeL7 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                7,
                tx + bx * ParGivenL7 + d1 * ParGivenL7 * NumOfBlockForEachNodeL7 + 1
            );
            for (int tmp = 0; tmp < 7; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            M2[0][0] = 1;
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[0][2] = C[NbrIdx[0] * n + NbrIdx[2]];
            M2[0][3] = C[NbrIdx[0] * n + NbrIdx[3]];
            M2[0][4] = C[NbrIdx[0] * n + NbrIdx[4]];
            M2[0][5] = C[NbrIdx[0] * n + NbrIdx[5]];
            M2[0][6] = C[NbrIdx[0] * n + NbrIdx[6]];

            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[1][2] = C[NbrIdx[1] * n + NbrIdx[2]];
            M2[1][3] = C[NbrIdx[1] * n + NbrIdx[3]];
            M2[1][4] = C[NbrIdx[1] * n + NbrIdx[4]];
            M2[1][5] = C[NbrIdx[1] * n + NbrIdx[5]];
            M2[1][6] = C[NbrIdx[1] * n + NbrIdx[6]];

            M2[2][0] = M2[0][2];
            M2[2][1] = M2[1][2];
            M2[2][2] = 1;
            M2[2][3] = C[NbrIdx[2] * n + NbrIdx[3]];
            M2[2][4] = C[NbrIdx[2] * n + NbrIdx[4]];
            M2[2][5] = C[NbrIdx[2] * n + NbrIdx[5]];
            M2[2][6] = C[NbrIdx[2] * n + NbrIdx[6]];

            M2[3][0] = M2[0][3];
            M2[3][1] = M2[1][3];
            M2[3][2] = M2[2][3];
            M2[3][3] = 1;
            M2[3][4] = C[NbrIdx[3] * n + NbrIdx[4]];
            M2[3][5] = C[NbrIdx[3] * n + NbrIdx[5]];
            M2[3][6] = C[NbrIdx[3] * n + NbrIdx[6]];

            M2[4][0] = M2[0][4];
            M2[4][1] = M2[1][4];
            M2[4][2] = M2[2][4];
            M2[4][3] = M2[3][4];
            M2[4][4] = 1;
            M2[4][5] = C[NbrIdx[4] * n + NbrIdx[5]];
            M2[4][6] = C[NbrIdx[4] * n + NbrIdx[6]];

            M2[5][0] = M2[0][5];
            M2[5][1] = M2[1][5];
            M2[5][2] = M2[2][5];
            M2[5][3] = M2[3][5];
            M2[5][4] = M2[4][5];
            M2[5][5] = 1;
            M2[5][6] = C[NbrIdx[5] * n + NbrIdx[6]];

            M2[6][0] = M2[0][6];
            M2[6][1] = M2[1][6];
            M2[6][2] = M2[2][6];
            M2[6][3] = M2[3][6];
            M2[6][4] = M2[4][6];
            M2[6][5] = M2[5][6];
            M2[6][6] = 1;
            pseudoinversel7(M2, M2Inv, v, rv1, w, res1);

            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][2] = C[XIdx * n + NbrIdx[2]];
            M1[0][3] = C[XIdx * n + NbrIdx[3]];
            M1[0][4] = C[XIdx * n + NbrIdx[4]];
            M1[0][5] = C[XIdx * n + NbrIdx[5]];
            M1[0][6] = C[XIdx * n + NbrIdx[6]];

            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation

                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    M1[1][2] = C[YIdx * n + NbrIdx[2]];
                    M1[1][3] = C[YIdx * n + NbrIdx[3]];
                    M1[1][4] = C[YIdx * n + NbrIdx[4]];
                    M1[1][5] = C[YIdx * n + NbrIdx[5]];
                    M1[1][6] = C[YIdx * n + NbrIdx[6]];
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 7; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 7; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 7; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 9, n) - 7.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl8_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 8;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[8];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[8];
    int var_ixs[10];

    float M0;
    float M1[2][8];
    float M2[8][8];
    float M2Inv[8][8];
    float M1MulM2Inv[2][8];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[8][8];
    float w[8], rv1[8];
    float res1[8][8];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 8)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL8) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL8;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL8 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL8) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL8] = GPrime[XIdx * n + tx + cnt * ParGivenL8];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 8, &NumOfComb);
    if ((NumOfComb % (ParGivenL8 * NumOfBlockForEachNodeL8)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL8 * NumOfBlockForEachNodeL8);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL8 * NumOfBlockForEachNodeL8) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL8 + d1 * ParGivenL8 * NumOfBlockForEachNodeL8 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                8,
                tx + bx * ParGivenL8 + d1 * ParGivenL8 * NumOfBlockForEachNodeL8 + 1
            );
            for (int tmp = 0; tmp < 8; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            M2[0][0] = 1;
            M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
            M2[0][2] = C[NbrIdx[0] * n + NbrIdx[2]];
            M2[0][3] = C[NbrIdx[0] * n + NbrIdx[3]];
            M2[0][4] = C[NbrIdx[0] * n + NbrIdx[4]];
            M2[0][5] = C[NbrIdx[0] * n + NbrIdx[5]];
            M2[0][6] = C[NbrIdx[0] * n + NbrIdx[6]];
            M2[0][7] = C[NbrIdx[0] * n + NbrIdx[7]];

            M2[1][0] = M2[0][1];
            M2[1][1] = 1;
            M2[1][2] = C[NbrIdx[1] * n + NbrIdx[2]];
            M2[1][3] = C[NbrIdx[1] * n + NbrIdx[3]];
            M2[1][4] = C[NbrIdx[1] * n + NbrIdx[4]];
            M2[1][5] = C[NbrIdx[1] * n + NbrIdx[5]];
            M2[1][6] = C[NbrIdx[1] * n + NbrIdx[6]];
            M2[1][7] = C[NbrIdx[1] * n + NbrIdx[7]];

            M2[2][0] = M2[0][2];
            M2[2][1] = M2[1][2];
            M2[2][2] = 1;
            M2[2][3] = C[NbrIdx[2] * n + NbrIdx[3]];
            M2[2][4] = C[NbrIdx[2] * n + NbrIdx[4]];
            M2[2][5] = C[NbrIdx[2] * n + NbrIdx[5]];
            M2[2][6] = C[NbrIdx[2] * n + NbrIdx[6]];
            M2[2][7] = C[NbrIdx[2] * n + NbrIdx[7]];

            M2[3][0] = M2[0][3];
            M2[3][1] = M2[1][3];
            M2[3][2] = M2[2][3];
            M2[3][3] = 1;
            M2[3][4] = C[NbrIdx[3] * n + NbrIdx[4]];
            M2[3][5] = C[NbrIdx[3] * n + NbrIdx[5]];
            M2[3][6] = C[NbrIdx[3] * n + NbrIdx[6]];
            M2[3][7] = C[NbrIdx[3] * n + NbrIdx[7]];

            M2[4][0] = M2[0][4];
            M2[4][1] = M2[1][4];
            M2[4][2] = M2[2][4];
            M2[4][3] = M2[3][4];
            M2[4][4] = 1;
            M2[4][5] = C[NbrIdx[4] * n + NbrIdx[5]];
            M2[4][6] = C[NbrIdx[4] * n + NbrIdx[6]];
            M2[4][7] = C[NbrIdx[4] * n + NbrIdx[7]];

            M2[5][0] = M2[0][5];
            M2[5][1] = M2[1][5];
            M2[5][2] = M2[2][5];
            M2[5][3] = M2[3][5];
            M2[5][4] = M2[4][5];
            M2[5][5] = 1;
            M2[5][6] = C[NbrIdx[5] * n + NbrIdx[6]];
            M2[5][7] = C[NbrIdx[5] * n + NbrIdx[7]];

            M2[6][0] = M2[0][6];
            M2[6][1] = M2[1][6];
            M2[6][2] = M2[2][6];
            M2[6][3] = M2[3][6];
            M2[6][4] = M2[4][6];
            M2[6][5] = M2[5][6];
            M2[6][6] = 1;
            M2[6][7] = C[NbrIdx[6] * n + NbrIdx[7]];

            M2[7][0] = M2[0][7];
            M2[7][1] = M2[1][7];
            M2[7][2] = M2[2][7];
            M2[7][3] = M2[3][7];
            M2[7][4] = M2[4][7];
            M2[7][5] = M2[5][7];
            M2[7][6] = M2[6][7];
            M2[7][7] = 1;

            M1[0][0] = C[XIdx * n + NbrIdx[0]];
            M1[0][1] = C[XIdx * n + NbrIdx[1]];
            M1[0][2] = C[XIdx * n + NbrIdx[2]];
            M1[0][3] = C[XIdx * n + NbrIdx[3]];
            M1[0][4] = C[XIdx * n + NbrIdx[4]];
            M1[0][5] = C[XIdx * n + NbrIdx[5]];
            M1[0][6] = C[XIdx * n + NbrIdx[6]];
            M1[0][7] = C[XIdx * n + NbrIdx[7]];
            pseudoinversel8(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    M1[1][0] = C[YIdx * n + NbrIdx[0]];
                    M1[1][1] = C[YIdx * n + NbrIdx[1]];
                    M1[1][2] = C[YIdx * n + NbrIdx[2]];
                    M1[1][3] = C[YIdx * n + NbrIdx[3]];
                    M1[1][4] = C[YIdx * n + NbrIdx[4]];
                    M1[1][5] = C[YIdx * n + NbrIdx[5]];
                    M1[1][6] = C[YIdx * n + NbrIdx[6]];
                    M1[1][7] = C[YIdx * n + NbrIdx[7]];
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 8; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 8; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 8; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 10, n) - 8.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl9_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 9;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[9];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[9];
    int var_ixs[11];

    float M0;
    float M1[2][9];
    float M2[9][9];
    float M2Inv[9][9];
    float M1MulM2Inv[2][9];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[9][9];
    float w[9], rv1[9];
    float res1[9][9];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 9)
    {
        return;
    }
    if ((SizeOfArr % ParGivenL9) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL9;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL9 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL9) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL9] = GPrime[XIdx * n + tx + cnt * ParGivenL9];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 9, &NumOfComb);
    if ((NumOfComb % (ParGivenL9 * NumOfBlockForEachNodeL9)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL9 * NumOfBlockForEachNodeL9);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL9 * NumOfBlockForEachNodeL9) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL9 + d1 * ParGivenL9 * NumOfBlockForEachNodeL9 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                9,
                tx + bx * ParGivenL9 + d1 * ParGivenL9 * NumOfBlockForEachNodeL9 + 1
            );
            for (int tmp = 0; tmp < 9; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            for (int c1 = 0; c1 < 9; c1++)
            {
                for (int c2 = 0; c2 < 9; c2++)
                {
                    if (c1 > c2)
                    {
                        M2[c1][c2] = M2[c2][c1];
                    }
                    else if (c1 == c2)
                    {
                        M2[c1][c1] = 1;
                    }
                    else
                    {
                        M2[c1][c2] = C[NbrIdx[c1] * n + NbrIdx[c2]];
                    }
                }
            }

            for (int c1 = 0; c1 < 9; c1++)
            {
                M1[0][c1] = C[XIdx * n + NbrIdx[c1]];
            }

            pseudoinversel9(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)) ||
                    (d2 == (NbrIdxPointer[8] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    for (int c1 = 0; c1 < 9; c1++)
                    {
                        M1[1][c1] = C[YIdx * n + NbrIdx[c1]];
                    }
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 9; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 9; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 9; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    var_ixs[10] = NbrIdx[8];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 11, n) - 9.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl10_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 10;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[10];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[10];
    int var_ixs[12];

    float M0;
    float M1[2][10];
    float M2[10][10];
    float M2Inv[10][10];
    float M1MulM2Inv[2][10];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[10][10];
    float w[10], rv1[10];
    float res1[10][10];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 10)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL10) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL10;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL10 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL10) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL10] = GPrime[XIdx * n + tx + cnt * ParGivenL10];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 10, &NumOfComb);
    if ((NumOfComb % (ParGivenL10 * NumOfBlockForEachNodeL10)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL10 * NumOfBlockForEachNodeL10);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL10 * NumOfBlockForEachNodeL10) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL10 + d1 * ParGivenL10 * NumOfBlockForEachNodeL10 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                10,
                tx + bx * ParGivenL10 + d1 * ParGivenL10 * NumOfBlockForEachNodeL10 + 1
            );
            for (int tmp = 0; tmp < 10; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            for (int c1 = 0; c1 < 10; c1++)
            {
                for (int c2 = 0; c2 < 10; c2++)
                {
                    if (c1 > c2)
                    {
                        M2[c1][c2] = M2[c2][c1];
                    }
                    else if (c1 == c2)
                    {
                        M2[c1][c1] = 1;
                    }
                    else
                    {
                        M2[c1][c2] = C[NbrIdx[c1] * n + NbrIdx[c2]];
                    }
                }
            }

            for (int c1 = 0; c1 < 10; c1++)
            {
                M1[0][c1] = C[XIdx * n + NbrIdx[c1]];
            }

            pseudoinversel10(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)) ||
                    (d2 == (NbrIdxPointer[8] - 1)) || (d2 == (NbrIdxPointer[9] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    for (int c1 = 0; c1 < 10; c1++)
                    {
                        M1[1][c1] = C[YIdx * n + NbrIdx[c1]];
                    }
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 10; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 10; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 10; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    var_ixs[10] = NbrIdx[8];
                    var_ixs[11] = NbrIdx[9];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 12, n) - 10.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl11_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 11;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[11];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[11];
    int var_ixs[13];

    float M0;
    float M1[2][11];
    float M2[11][11];
    float M2Inv[11][11];
    float M1MulM2Inv[2][11];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[11][11];
    float w[11], rv1[11];
    float res1[11][11];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 11)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL11) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL11;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL11 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL11) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL11] = GPrime[XIdx * n + tx + cnt * ParGivenL11];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 11, &NumOfComb);
    if ((NumOfComb % (ParGivenL11 * NumOfBlockForEachNodeL11)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL11 * NumOfBlockForEachNodeL11);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL11 * NumOfBlockForEachNodeL11) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL11 + d1 * ParGivenL11 * NumOfBlockForEachNodeL11 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                11,
                tx + bx * ParGivenL11 + d1 * ParGivenL11 * NumOfBlockForEachNodeL11 + 1
            );
            for (int tmp = 0; tmp < 11; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            for (int c1 = 0; c1 < 11; c1++)
            {
                for (int c2 = 0; c2 < 11; c2++)
                {
                    if (c1 > c2)
                    {
                        M2[c1][c2] = M2[c2][c1];
                    }
                    else if (c1 == c2)
                    {
                        M2[c1][c1] = 1;
                    }
                    else
                    {
                        M2[c1][c2] = C[NbrIdx[c1] * n + NbrIdx[c2]];
                    }
                }
            }

            for (int c1 = 0; c1 < 11; c1++)
            {
                M1[0][c1] = C[XIdx * n + NbrIdx[c1]];
            }

            pseudoinversel11(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)) ||
                    (d2 == (NbrIdxPointer[8] - 1)) || (d2 == (NbrIdxPointer[9] - 1)) ||
                    (d2 == (NbrIdxPointer[10] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    for (int c1 = 0; c1 < 11; c1++)
                    {
                        M1[1][c1] = C[YIdx * n + NbrIdx[c1]];
                    }
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 11; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 11; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 11; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    var_ixs[10] = NbrIdx[8];
                    var_ixs[11] = NbrIdx[9];
                    var_ixs[12] = NbrIdx[10];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 13, n) - 11.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl12_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 12;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[12];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[12];
    int var_ixs[14];

    float M0;
    float M1[2][12];
    float M2[12][12];
    float M2Inv[12][12];
    float M1MulM2Inv[2][12];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[12][12];
    float w[12], rv1[12];
    float res1[12][12];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 12)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL12) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL12;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL12 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL12) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL12] = GPrime[XIdx * n + tx + cnt * ParGivenL12];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 12, &NumOfComb);
    if ((NumOfComb % (ParGivenL12 * NumOfBlockForEachNodeL12)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL12 * NumOfBlockForEachNodeL12);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL12 * NumOfBlockForEachNodeL12) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL12 + d1 * ParGivenL12 * NumOfBlockForEachNodeL12 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                12,
                tx + bx * ParGivenL12 + d1 * ParGivenL12 * NumOfBlockForEachNodeL12 + 1
            );
            for (int tmp = 0; tmp < 12; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            for (int c1 = 0; c1 < 12; c1++)
            {
                for (int c2 = 0; c2 < 12; c2++)
                {
                    if (c1 > c2)
                    {
                        M2[c1][c2] = M2[c2][c1];
                    }
                    else if (c1 == c2)
                    {
                        M2[c1][c1] = 1;
                    }
                    else
                    {
                        M2[c1][c2] = C[NbrIdx[c1] * n + NbrIdx[c2]];
                    }
                }
            }

            for (int c1 = 0; c1 < 12; c1++)
            {
                M1[0][c1] = C[XIdx * n + NbrIdx[c1]];
            }

            pseudoinversel12(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)) ||
                    (d2 == (NbrIdxPointer[8] - 1)) || (d2 == (NbrIdxPointer[9] - 1)) ||
                    (d2 == (NbrIdxPointer[10] - 1)) || (d2 == (NbrIdxPointer[11] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    for (int c1 = 0; c1 < 12; c1++)
                    {
                        M1[1][c1] = C[YIdx * n + NbrIdx[c1]];
                    }
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 12; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 12; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 12; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    var_ixs[10] = NbrIdx[8];
                    var_ixs[11] = NbrIdx[9];
                    var_ixs[12] = NbrIdx[10];
                    var_ixs[13] = NbrIdx[11];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 14, n) - 12.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl13_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 13;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[13];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[13];
    int var_ixs[15];

    float M0;
    float M1[2][13];
    float M2[13][13];
    float M2Inv[13][13];
    float M1MulM2Inv[2][13];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[13][13];
    float w[13], rv1[13];
    float res1[13][13];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 13)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL13) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL13;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL13 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL13) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL13] = GPrime[XIdx * n + tx + cnt * ParGivenL13];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 13, &NumOfComb);
    if ((NumOfComb % (ParGivenL13 * NumOfBlockForEachNodeL13)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL13 * NumOfBlockForEachNodeL13);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL13 * NumOfBlockForEachNodeL13) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL13 + d1 * ParGivenL13 * NumOfBlockForEachNodeL13 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                13,
                tx + bx * ParGivenL13 + d1 * ParGivenL13 * NumOfBlockForEachNodeL13 + 1
            );
            for (int tmp = 0; tmp < 13; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            for (int c1 = 0; c1 < 13; c1++)
            {
                for (int c2 = 0; c2 < 13; c2++)
                {
                    if (c1 > c2)
                    {
                        M2[c1][c2] = M2[c2][c1];
                    }
                    else if (c1 == c2)
                    {
                        M2[c1][c1] = 1;
                    }
                    else
                    {
                        M2[c1][c2] = C[NbrIdx[c1] * n + NbrIdx[c2]];
                    }
                }
            }

            for (int c1 = 0; c1 < 13; c1++)
            {
                M1[0][c1] = C[XIdx * n + NbrIdx[c1]];
            }

            pseudoinversel13(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)) ||
                    (d2 == (NbrIdxPointer[8] - 1)) || (d2 == (NbrIdxPointer[9] - 1)) ||
                    (d2 == (NbrIdxPointer[10] - 1)) || (d2 == (NbrIdxPointer[11] - 1)) ||
                    (d2 == (NbrIdxPointer[12] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    for (int c1 = 0; c1 < 13; c1++)
                    {
                        M1[1][c1] = C[YIdx * n + NbrIdx[c1]];
                    }
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 13; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 13; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 13; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    var_ixs[10] = NbrIdx[8];
                    var_ixs[11] = NbrIdx[9];
                    var_ixs[12] = NbrIdx[10];
                    var_ixs[13] = NbrIdx[11];
                    var_ixs[14] = NbrIdx[12];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 15, n) - 13.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl14_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
)
{
    int level = 14;
    float loc_th;
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer[14];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    int NbrIdx[14];
    int var_ixs[16];

    float M0;
    float M1[2][14];
    float M2[14][14];
    float M2Inv[14][14];
    float M1MulM2Inv[2][14];
    float H[2][2];
    float rho;
    float Z;
    // pseudo-inverse parameter
    float v[14][14];
    float w[14], rv1[14];
    float res1[14][14];

    extern __shared__ int G_Chunk[];

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];
    if (SizeOfArr <= 14)
    {
        return;
    }

    if ((SizeOfArr % ParGivenL14) == 0)
    {
        NumberOfJump = SizeOfArr / ParGivenL14;
    }
    else
    {
        NumberOfJump = SizeOfArr / ParGivenL14 + 1;
    }
    // Copy Row Xid from GPrime to G_chunck
    for (int cnt = 0; cnt < NumberOfJump; cnt++)
    {
        if ((tx + cnt * ParGivenL14) < SizeOfArr)
        {
            G_Chunk[tx + cnt * ParGivenL14] = GPrime[XIdx * n + tx + cnt * ParGivenL14];
        }
        __syncthreads();
    }

    BINOM(SizeOfArr, 14, &NumOfComb);
    if ((NumOfComb % (ParGivenL14 * NumOfBlockForEachNodeL14)) == 0)
    {
        NumOfGivenJump = NumOfComb / (ParGivenL14 * NumOfBlockForEachNodeL14);
    }
    else
    {
        NumOfGivenJump = NumOfComb / (ParGivenL14 * NumOfBlockForEachNodeL14) + 1;
    }

    for (int d1 = 0; d1 < NumOfGivenJump; d1++)
    {
        __syncthreads();
        if (NoEdgeFlag == 1)
        {
            return;
        }
        if (tx + bx * ParGivenL14 + d1 * ParGivenL14 * NumOfBlockForEachNodeL14 < NumOfComb)
        {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            IthCombination(
                NbrIdxPointer,
                SizeOfArr,
                14,
                tx + bx * ParGivenL14 + d1 * ParGivenL14 * NumOfBlockForEachNodeL14 + 1
            );
            for (int tmp = 0; tmp < 14; tmp++)
            {
                NbrIdx[tmp] = G_Chunk[NbrIdxPointer[tmp] - 1];
            }

            for (int c1 = 0; c1 < 14; c1++)
            {
                for (int c2 = 0; c2 < 14; c2++)
                {
                    if (c1 > c2)
                    {
                        M2[c1][c2] = M2[c2][c1];
                    }
                    else if (c1 == c2)
                    {
                        M2[c1][c1] = 1;
                    }
                    else
                    {
                        M2[c1][c2] = C[NbrIdx[c1] * n + NbrIdx[c2]];
                    }
                }
            }

            for (int c1 = 0; c1 < 14; c1++)
            {
                M1[0][c1] = C[XIdx * n + NbrIdx[c1]];
            }

            pseudoinversel14(M2, M2Inv, v, rv1, w, res1);
            for (int d2 = 0; d2 < SizeOfArr; d2++)
            {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1)) ||
                    (d2 == (NbrIdxPointer[2] - 1)) || (d2 == (NbrIdxPointer[3] - 1)) ||
                    (d2 == (NbrIdxPointer[4] - 1)) || (d2 == (NbrIdxPointer[5] - 1)) ||
                    (d2 == (NbrIdxPointer[6] - 1)) || (d2 == (NbrIdxPointer[7] - 1)) ||
                    (d2 == (NbrIdxPointer[8] - 1)) || (d2 == (NbrIdxPointer[9] - 1)) ||
                    (d2 == (NbrIdxPointer[10] - 1)) || (d2 == (NbrIdxPointer[11] - 1)) ||
                    (d2 == (NbrIdxPointer[12] - 1)) || (d2 == (NbrIdxPointer[13] - 1)))
                {
                    continue;
                }
                YIdx = G_Chunk[d2];
                if (!valid_time_conditioning(XIdx, YIdx, NbrIdx, time_index, level)) {
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1)
                {
                    NoEdgeFlag = 0;
                    M0 = C[XIdx * n + YIdx];
                    // Beginning Of the Indep Test Calculation
                    for (int c1 = 0; c1 < 14; c1++)
                    {
                        M1[1][c1] = C[YIdx * n + NbrIdx[c1]];
                    }
                    // Begin to calculate I2Inv Using pseudo-inverse
                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 14; c2++)
                        {
                            M1MulM2Inv[c1][c2] = 0;
                            for (int c3 = 0; c3 < 14; c3++)
                                M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                        }
                    }

                    for (int c1 = 0; c1 < 2; c1++)
                    {
                        for (int c2 = 0; c2 < 2; c2++)
                        {
                            H[c1][c2] = 0;
                            for (int c3 = 0; c3 < 14; c3++)
                                H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                        }
                    }
                    H[0][0] = 1 - H[0][0];
                    H[0][1] = M0 - H[0][1];
                    H[1][1] = 1 - H[1][1];

                    var_ixs[0] = XIdx;
                    var_ixs[1] = YIdx;
                    var_ixs[2] = NbrIdx[0];
                    var_ixs[3] = NbrIdx[1];
                    var_ixs[4] = NbrIdx[2];
                    var_ixs[5] = NbrIdx[3];
                    var_ixs[6] = NbrIdx[4];
                    var_ixs[7] = NbrIdx[5];
                    var_ixs[8] = NbrIdx[6];
                    var_ixs[9] = NbrIdx[7];
                    var_ixs[10] = NbrIdx[8];
                    var_ixs[11] = NbrIdx[9];
                    var_ixs[12] = NbrIdx[10];
                    var_ixs[13] = NbrIdx[11];
                    var_ixs[14] = NbrIdx[12];
                    var_ixs[15] = NbrIdx[13];
                    loc_th = th / sqrt(mean_ess(N, var_ixs, 16, n) - 14.0 - 3.0);

                    rho = H[0][1] / (sqrt(abs(H[0][0] * H[1][1])));
                    Z = abs(0.5 * log(abs((1 + rho) / (1 - rho))));

                    if (Z < loc_th)
                    {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0)
                        {  // lock
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            print_sepset(var_ixs, time_index, level);
                        }
                    }
                }
            }
        }
    }
}


__device__ bool valid_time_conditioning(int a, int b, int *S, int *time_index, int l)
{
    int maxtime = max(time_index[a], time_index[b]);
    for (int nix = 0; nix < l; nix++) {
        if (time_index[S[nix]] > maxtime)
        {
            return false;
        }    
    }
    return true;
}


__device__ float mean_ess(float *N, int var_ixs[], int l, int n)
{
    float s = 0.0;
    int ix_a;
    int ix_b;
    int loc_val;
    int num_s = 0;
    for (int i = 0; i < l; i++)
    {
        ix_a = var_ixs[i];
        for (int j = 0; j < i; j++) {
            ix_b = var_ixs[j];
            loc_val = N[ix_a * n + ix_b];
            if (!isnan(loc_val)) {
                s += loc_val;
                num_s += 1;
            }
        }
    }
    float res = s / (float)num_s;
    return res;
}