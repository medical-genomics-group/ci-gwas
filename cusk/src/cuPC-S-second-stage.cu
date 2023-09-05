#include <math.h>
#include <mps/cuPC-S.h>
#include <mps/gpuerrors.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>

/*
NOTES REGARDING MODIFICATIONS

PURPOSE:
The idea of the second stage is to re-run cuPC with a reduced set of variables, with the aim to
learn better separation sets. Single-stage cuPC has the problems that i) the separation sets are
order-dependent ii) when run at an alpha that controls the FDR for adjacencies, many links are
removed due to small marginal correlations. Therefore, we run stage 1 at a small alpha to find the
skeleton, and stage 2 with a large alpha to find more conditional independencies (as opposed to
marginal ones).
//

DESCRIPTION:
At each level, we compare the partial correlations (pcorr) found with all separation sets, i.e. we
don't stop at the first one that is sufficient for independence. We then pick the one that reduces
the pcorr the most.
If some separation set is found, the link is removed from the graph.
However, it is still checked in the next level, but in a modified routine, equivalent to l1:
we add one variable at a time to the found separation set, recording the pcorrs, and again choosing
the set with the lowest pcorr. If no new set has a lower pcorr than the one from the previous level,
this link has found its final separation set and does not need to be updated anymore.
As long as not separation set was found, a link remains in the graph and the search for a sepset
follows the usual routine of checking all combinations at a level.

So we have one new routine which handles the l1-like sepset extensions for all links that have a
(non-empty) sepset already. This is called before the normal l2, l3, etc.

The normal l1, l2, l3... routines have to be modified:
the local minimum for each thread should be traced, then compared among all threads in a block,
and then among all blocks.
The easiest way to do this seems to be to have a datastructure
pcorr_cuda of dims n * n * max-degree-after-l0 (float);
I hope that max-degree-after-l0 isn't too large, otherwise this might get large.
The mutex section is completely removed from the routines, they simply write the pcorr
obtained to pcorr_cuda, which is then reduced in a final step.
The reduction can probably be done using scan/reduce, but I would probably in my first iteration
of the impl just send one thread per edge.
//
*/

//========================> Main Function Parameter <========================
// Description : this function just calculate one Stage of PC stable algorithm
//@param C          = Correlation matrix
//@param VarSize    = Number of Nodes in Dataset
//@param Stage      = Number of Neighbor in each dimension of Neighbor Matrix
//@param G          = Is the Graph array
//@param TH         = The Th for deleting each edge
//@param Nbr        = Neighbor Matrix with format of:
//[i , j , |Neighbor idx 1|,|Neighbor idx 2| , ...]
//@param Nrow       = Number Of row in Nbr matrix
//@param Ncol       = Number of Col in Nbr matrix
//============================================================================

__device__ int is_in_arr(int x, int *arr, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (arr[i] == x)
        {
            return 1;
        }
    }
    return 0;
}

__global__ void print_matrix(float *M, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", M[i * n + j]);
        }
        printf("\n");
    }
}

__global__ void print_matrix(int *M, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%d ", M[i * n + j]);
        }
        printf("\n");
    }
}

__global__ void unfinished_initialize(int *uf, int n)
{
    if (bx == by)
    {
        uf[bx * n + by] = 0;
    }
    else
    {
        uf[bx * n + by] = 1;
    }
}

__global__ void pMax_initialize(float *pMax, int n) { pMax[bx * n + by] = 1.0; }

__global__ void pcorr_initialize(float *pcorr, int n)
{
    pcorr[(bx * n + by) * PCORR_MAX_DEGREE + tx] = 1.0;
}

void cusk_second_stage(
    float *C, int *P, int *G, float *Th, int *l, const int *maxlevel, float *pMax, int *SepSet
)
{
    float *C_cuda;  // Copy of C array in GPU
    float *pMax_cuda;
    int *G_cuda;  // Copy of G Array in GPU
    int *nprime_cuda;
    int *SepSet_cuda;
    int *GPrime_cuda;
    int *mutex_cuda;
    float *pcorr_cuda;
    int *unfinished_cuda;
    int *unfinished_prime_cuda;

    int n = *P;
    int nprime = 0;
    dim3 BLOCKS_PER_GRID;
    dim3 THREADS_PER_BLOCK;

    bool FinishFlag = false;

    *l = 0;

    // initialize sepset element selection matrix
    HANDLE_ERROR(cudaMalloc((void **)&pcorr_cuda, n * n * PCORR_MAX_DEGREE * sizeof(float)));
    // marks which sepsets should be updated
    HANDLE_ERROR(cudaMalloc((void **)&unfinished_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&unfinished_prime_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&mutex_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&mutex_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&nprime_cuda, 1 * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&SepSet_cuda, n * n * ML * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&GPrime_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&C_cuda, n * n * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void **)&G_cuda, n * n * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&pMax_cuda, n * n * sizeof(float)));
    // copy skeleton from CPU to GPU
    HANDLE_ERROR(cudaMemcpy(G_cuda, G, n * n * sizeof(int), cudaMemcpyHostToDevice));
    // copy correlation matrix from CPU to GPU
    HANDLE_ERROR(cudaMemcpy(C_cuda, C, n * n * sizeof(float), cudaMemcpyHostToDevice));
    // initialize a 0 matrix
    HANDLE_ERROR(cudaMemset(mutex_cuda, 0, n * n * sizeof(int)));
    pMax_initialize<<<dim3(n, n, 1), dim3(1, 1, 1)>>>(pMax_cuda, n);
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
                marginal_pMax<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                    C_cuda, G_cuda, Th[0], pMax_cuda, n
                );
                CudaCheckError();
            }
            else
            {
                BLOCKS_PER_GRID = dim3(ceil(((float)(n)) / 32.0), ceil(((float)(n)) / 32.0), 1);
                THREADS_PER_BLOCK = dim3(32, 32, 1);
                marginal_pMax<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                    C_cuda, G_cuda, Th[0], pMax_cuda, n
                );
                CudaCheckError();
            }
            BLOCKS_PER_GRID = dim3(n * n, 1, 1);
            THREADS_PER_BLOCK = dim3(ML, 1, 1);
            SepSet_initialize<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(SepSet_cuda, n);
            CudaCheckError();
            pcorr_initialize<<<dim3(n, n, 1), dim3(PCORR_MAX_DEGREE, 1, 1)>>>(pcorr_cuda, n);
            CudaCheckError();
            unfinished_initialize<<<dim3(n, n, 1), dim3(1, 1, 1)>>>(unfinished_cuda, n);
            CudaCheckError();
        }
        else
        {
            //================================> Start Scan Process <===============================
            printf("Compacting adjacencies\n");
            fflush(stdout);
            HANDLE_ERROR(cudaMemset(nprime_cuda, 0, 1 * sizeof(int)));
            BLOCKS_PER_GRID = dim3(1, n, 1);
            THREADS_PER_BLOCK = dim3(1024, 1, 1);
            scan_compact<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, n * sizeof(int)>>>(
                GPrime_cuda, G_cuda, n, nprime_cuda
            );
            CudaCheckError();
            HANDLE_ERROR(cudaMemcpy(&nprime, nprime_cuda, 1 * sizeof(int), cudaMemcpyDeviceToHost));

            // printf("pMax_cuda: \n");
            // print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(pMax_cuda, n);
            // cudaDeviceSynchronize();

            // printf("G_cuda: \n");
            // print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(G_cuda, n);
            // cudaDeviceSynchronize();

            // printf("GPrime_cuda: \n");
            // print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(GPrime_cuda, n);
            // cudaDeviceSynchronize();

            // Check if the max degree is too large
            if (nprime > PCORR_MAX_DEGREE)
            {
                printf("max degree exceeds allowed value\n");
                fflush(stdout);
                break;
            }

            printf("Compacting lists of unfinished sepsets\n");
            fflush(stdout);
            // compact the list of adjacencies for which sepsets should be updated
            scan_compact<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, n * sizeof(int)>>>(
                unfinished_prime_cuda, unfinished_cuda, n, nprime_cuda
            );
            CudaCheckError();

            // printf("unfinished_prime_cuda: \n");
            // print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(unfinished_prime_cuda, n);
            // cudaDeviceSynchronize();

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
                printf("Calculating l1\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL1, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL1, 1, 1);
                check_sepsets_l1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }

            // find best sepset candidate, remove edges with sepsets that
            // yield Z below threshold
            printf("Selecting likely non-colliders\n");
            fflush(stdout);
            BLOCKS_PER_GRID = dim3(n, n, 1);
            THREADS_PER_BLOCK = dim3(1, 1, 1);
            select_non_colliders<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                G_cuda, GPrime_cuda, SepSet_cuda, pMax_cuda, pcorr_cuda, Th[*l], *l, n
            );
            CudaCheckError();
        }
    }  // if l > 0

    // Copy Graph G from GPU to CPU
    HANDLE_ERROR(cudaMemcpy(G, G_cuda, n * n * sizeof(int), cudaMemcpyDeviceToHost));
    // Copy separation set from GPU to CPU
    HANDLE_ERROR(cudaMemcpy(SepSet, SepSet_cuda, n * n * ML * sizeof(int), cudaMemcpyDeviceToHost));
    // Copy  Pmax from GPU to CPU
    HANDLE_ERROR(cudaMemcpy(pMax, pMax_cuda, n * n * sizeof(float), cudaMemcpyDeviceToHost));
    // Preprocess pMax
    float temp = 0;
    for (int i = 0; i < n; i++)
    {
        pMax[i * n + i] = 1;
        for (int j = (i + 1); j < n; j++)
        {
            if (G[i * n + j] == 0)
            {
                temp = fmax(pMax[j * n + i], pMax[i * n + j]);
                pMax[j * n + i] = temp;
                pMax[i * n + j] = temp;
            }
            else
            {
                pMax[j * n + i] = -100000;
                pMax[i * n + j] = -100000;
            }
        }
    }
    // Free allocated space
    HANDLE_ERROR(cudaFree(SepSet_cuda));
    HANDLE_ERROR(cudaFree(C_cuda));
    HANDLE_ERROR(cudaFree(GPrime_cuda));
    HANDLE_ERROR(cudaFree(G_cuda));
    HANDLE_ERROR(cudaFree(mutex_cuda));
    HANDLE_ERROR(cudaFree(pMax_cuda));
    HANDLE_ERROR(cudaFree(unfinished_cuda));
    HANDLE_ERROR(cudaFree(unfinished_prime_cuda));
    HANDLE_ERROR(cudaFree(pcorr_cuda));
}  // Skeleton

__global__ void select_non_colliders(
    int *G, int *GPrime, int *Sepset, float *pMax, float *pcorrs, float th, int l, int n
)
{
    int XIdx = bx;
    int YIdx = by;
    bool new_min = false;
    int SizeOfArr = GPrime[XIdx * n + n - 1];
    // set to previous pMax (the marginal correlation)
    float marginal_corr = pMax[XIdx * n + YIdx];
    int sepset_pos = 0;

    for (int NbrIdxPointer = 0; NbrIdxPointer < SizeOfArr; NbrIdxPointer++)
    {
        float cond_corr = pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer];
        if (cond_corr < marginal_corr)
        {
            Sepset[(XIdx * n + YIdx) * ML + sepset_pos] = GPrime[XIdx * n + NbrIdxPointer];
            sepset_pos++;
        }
    }
}

__global__ void find_min_pcorr(
    int *G,
    int *GPrime,
    int *Sepset,
    float *pMax,
    float *pcorrs,
    int *unfinished,
    float th,
    int l,
    int n
)
{
    int min_nbr_idx;
    int XIdx = bx;
    int YIdx = by;
    if (unfinished[XIdx * n + YIdx] == 0)
    {
        return;
    }
    bool new_min = false;
    int SizeOfArr = GPrime[XIdx * n + n - 1];
    // set to previous pMax
    float min_z = pMax[XIdx * n + YIdx];

    for (int NbrIdxPointer = 0; NbrIdxPointer < SizeOfArr; NbrIdxPointer++)
    {
        float curr_z = pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer];
        if (curr_z < min_z)
        {
            new_min = true;
            min_z = curr_z;
            min_nbr_idx = GPrime[XIdx * n + NbrIdxPointer];
        }
    }

    if (new_min)
    {
        Sepset[(XIdx * n + YIdx) * ML + (l - 1)] = min_nbr_idx;
        pMax[XIdx * n + YIdx] = min_z;
        if (min_z < th)
        {
            G[XIdx * n + YIdx] = 0;
            G[YIdx * n + XIdx] = 0;
        }
    }
    else
    {
        // mark as finished
        unfinished[XIdx * n + YIdx] = 0;
    }
}

/*
Like cal_Indepl0, but without addition of edges for marginally-dependent pairs.
*/
__global__ void marginal_pMax(float *C, int *G, float th, float *pMax, int n)
{
    int row = blockDim.x * bx + tx;
    int col = blockDim.y * by + ty;
    if (row < col && col < n)
    {
        float res = C[row * n + col];
        res = abs(0.5 * log(abs((1 + res) / (1 - res))));
        if (res < th)
        {
            pMax[row * n + col] = res;
            pMax[col * n + row] = res;
            G[row * n + col] = 0;
            G[col * n + row] = 0;
        }
    }
    if (row == col && col < n)
    {
        G[row * n + col] = 0;
        G[col * n + row] = 0;
    }
}

__global__ void check_sepsets_l1(
    float *C, int *G, int *GPrime, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx;
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float H[2][2];
    float M1[2];
    float rho, Z;
    extern __shared__ int G_Chunk[];

    SizeOfArr = GPrime[XIdx * n + n - 1];
    num_unfinished_sepsets = unfinished_prime[XIdx * n + n - 1];
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
        NbrIdxPointer = tx + bx * ParGivenL1 + d1 * ParGivenL1 * NumOfBlockForEachNodeL1;
        if (NbrIdxPointer < SizeOfArr)
        {
            NbrIdx = G_Chunk[NbrIdxPointer];
            M1[0] = C[XIdx * n + NbrIdx];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[XIdx * n + d2];

                if (YIdx == NbrIdx)
                {
                    continue;
                }

                // printf("updating pcorr: x: %d, y: %d, z: %d \n", XIdx, YIdx, NbrIdx);

                M0 = C[XIdx * n + YIdx];
                M1[1] = C[YIdx * n + NbrIdx];

                H[0][0] = 1 - (M1[0] * M1[0]);
                H[0][1] = M0 - (M1[0] * M1[1]);
                H[1][1] = 1 - (M1[1] * M1[1]);

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}