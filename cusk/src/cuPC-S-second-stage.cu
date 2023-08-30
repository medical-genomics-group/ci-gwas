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
                cal_Indepl0<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                    C_cuda, G_cuda, Th[0], pMax_cuda, n
                );
                CudaCheckError();
            }
            else
            {
                BLOCKS_PER_GRID = dim3(ceil(((float)(n)) / 32.0), ceil(((float)(n)) / 32.0), 1);
                THREADS_PER_BLOCK = dim3(32, 32, 1);
                cal_Indepl0<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                    C_cuda, G_cuda, Th[0], pMax_cuda, n
                );
                CudaCheckError();
            }
            BLOCKS_PER_GRID = dim3(n * n, 1, 1);
            THREADS_PER_BLOCK = dim3(ML, 1, 1);
            SepSet_initialize<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(SepSet_cuda, n);
            CudaCheckError();
            pcorr_initialize<<<BLOCKS_PER_GRID, dim3(PCORR_MAX_DEGREE, 1, 1)>>>(pcorr_cuda, n);
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

            printf("pMax_cuda: \n");
            print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(pMax_cuda, n);
            cudaDeviceSynchronize();

            printf("G_cuda: \n");
            print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(G_cuda, n);
            cudaDeviceSynchronize();

            printf("GPrime_cuda: \n");
            print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(GPrime_cuda, n);
            cudaDeviceSynchronize();

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

            printf("unfinished_prime_cuda: \n");
            print_matrix<<<dim3(1, 1, 1), dim3(1, 1, 1)>>>(unfinished_prime_cuda, n);
            cudaDeviceSynchronize();

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
                check_sepsets_l1<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 2)
            {
                printf("Starting lvl 2\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL2, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL2, 1, 1);
                check_sepsets_l2<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 3)
            {
                printf("Starting lvl 3\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL3, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL3, 1, 1);
                check_sepsets_l3<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 4)
            {
                printf("Starting lvl 4\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL4, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL4, 1, 1);
                check_sepsets_l4<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 5)
            {
                printf("Starting lvl 5\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL5, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL5, 1, 1);
                check_sepsets_l5<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 6)
            {
                printf("Starting lvl 6\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL6, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL6, 1, 1);
                check_sepsets_l6<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 7)
            {
                printf("Starting lvl 7\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL7, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL7, 1, 1);
                check_sepsets_l7<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 8)
            {
                printf("Starting lvl 8\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL8, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL8, 1, 1);
                check_sepsets_l8<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 9)
            {
                printf("Starting lvl 9\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL9, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL9, 1, 1);
                check_sepsets_l9<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 10)
            {
                printf("Starting lvl 10\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL10, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL10, 1, 1);
                check_sepsets_l10<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 11)
            {
                printf("Starting lvl 11\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL11, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL11, 1, 1);
                check_sepsets_l11<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 12)
            {
                printf("Starting lvl 12\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL12, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL12, 1, 1);
                check_sepsets_l12<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 13)
            {
                printf("Starting lvl 13\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL13, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL13, 1, 1);
                check_sepsets_l13<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else if (*l == 14)
            {
                printf("Starting lvl 14\n");
                fflush(stdout);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL14, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL14, 1, 1);
                check_sepsets_l14<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int)>>>(
                    C_cuda, G_cuda, GPrime_cuda, SepSet_cuda, pcorr_cuda, unfinished_prime_cuda, n
                );
                CudaCheckError();
            }
            else
            {
                // TODO: add PC serial
            }

            // find best sepset candidate, remove edges with sepsets that
            // yield Z below threshold
            printf("Selecting sepsets with minimal pcorr\n");
            fflush(stdout);
            BLOCKS_PER_GRID = dim3(n, n, 1);
            THREADS_PER_BLOCK = dim3(1, 1, 1);
            find_min_pcorr<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(
                G_cuda,
                GPrime_cuda,
                SepSet_cuda,
                pMax_cuda,
                pcorr_cuda,
                unfinished_cuda,
                Th[*l],
                *l,
                n
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
            printf(
                "New min pcorr at x: %d, y: %d, nbrptr: %d, z: %f \n",
                XIdx,
                YIdx,
                NbrIdxPointer,
                curr_z
            );
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
                YIdx = unfinished_prime[d2];

                if (YIdx == NbrIdx)
                {
                    continue;
                }

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

__global__ void check_sepsets_l2(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[2];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][2];
    float M2[2][2];
    float M2Inv[2][2];
    float M1MulM2Inv[2][2];
    float H[2][2];
    float rho, Z;
    extern __shared__ int G_Chunk[];

    int pos_of_new_sepset_elem = 1;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
                }

                M2[0][1] = C[NbrIdx[0] * n + NbrIdx[1]];
                M2[1][0] = M2[0][1];
                M2[1][1] = 1;
                M2[0][0] = 1;

                M1[0][1] = C[XIdx * n + NbrIdx[1]];
                M1[0][0] = C[XIdx * n + NbrIdx[0]];
                pseudoinversel2(M2, M2Inv);

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

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l3(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[3];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][3];
    float M2[3][3];
    float M2Inv[3][3];
    float M1MulM2Inv[2][3];
    float H[2][2];
    float rho, Z;
    extern __shared__ int G_Chunk[];

    int pos_of_new_sepset_elem = 2;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
                }

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
                        for (int c3 = 0; c3 < 3; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l4(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[4];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][4];
    float M2[4][4];
    float M2Inv[4][4];
    float M1MulM2Inv[2][4];
    float H[2][2];
    float rho;
    float Z;

    float v[4][4];
    float w[4], rv1[4];
    float res1[4][4];
    extern __shared__ int G_Chunk[];

    int pos_of_new_sepset_elem = 3;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
                }

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
                        for (int c3 = 0; c3 < 4; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l5(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[5];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
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

    int pos_of_new_sepset_elem = 4;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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
                        for (int c3 = 0; c3 < 5; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l6(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[6];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
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

    int pos_of_new_sepset_elem = 5;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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
                        for (int c3 = 0; c3 < 6; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l7(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[7];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][7];
    float M2[7][7];
    float M2Inv[7][7];
    float M1MulM2Inv[2][7];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[7][7];
    float w[7], rv1[7];
    float res1[7][7];

    int pos_of_new_sepset_elem = 6;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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
                        for (int c3 = 0; c3 < 7; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l8(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[8];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][8];
    float M2[8][8];
    float M2Inv[8][8];
    float M1MulM2Inv[2][8];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[8][8];
    float w[8], rv1[8];
    float res1[8][8];

    int pos_of_new_sepset_elem = 7;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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
                        for (int c3 = 0; c3 < 8; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}
__global__ void check_sepsets_l9(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[9];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][9];
    float M2[9][9];
    float M2Inv[9][9];
    float M1MulM2Inv[2][9];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[9][9];
    float w[9], rv1[9];
    float res1[9][9];

    int pos_of_new_sepset_elem = 8;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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
                        for (int c3 = 0; c3 < 9; c3++) H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                    }
                }
                H[0][0] = 1 - H[0][0];
                H[0][1] = M0 - H[0][1];
                H[1][1] = 1 - H[1][1];

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l10(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[10];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][10];
    float M2[10][10];
    float M2Inv[10][10];
    float M1MulM2Inv[2][10];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[10][10];
    float w[10], rv1[10];
    float res1[10][10];

    int pos_of_new_sepset_elem = 9;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}
__global__ void check_sepsets_l11(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[11];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][11];
    float M2[11][11];
    float M2Inv[11][11];
    float M1MulM2Inv[2][11];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[11][11];
    float w[11], rv1[11];
    float res1[11][11];

    int pos_of_new_sepset_elem = 10;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
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

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l12(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[12];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][12];
    float M2[12][12];
    float M2Inv[12][12];
    float M1MulM2Inv[2][12];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[12][12];
    float w[12], rv1[12];
    float res1[12][12];

    int pos_of_new_sepset_elem = 11;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l13(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[13];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][13];
    float M2[13][13];
    float M2Inv[13][13];
    float M1MulM2Inv[2][13];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[13][13];
    float w[13], rv1[13];
    float res1[13][13];

    int pos_of_new_sepset_elem = 12;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}

__global__ void check_sepsets_l14(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
)
{
    int YIdx;
    int XIdx = by;
    int NbrIdxPointer;
    int NbrIdx[14];
    int SizeOfArr;
    int num_unfinished_sepsets;
    int NumberOfJump;
    int NumOfGivenJump;
    float M0;
    float M1[2][14];
    float M2[14][14];
    float M2Inv[14][14];
    float M1MulM2Inv[2][14];
    float H[2][2];
    float rho;
    float Z;
    extern __shared__ int G_Chunk[];
    // pseudo-inverse parameter
    float v[14][14];
    float w[14], rv1[14];
    float res1[14][14];

    int pos_of_new_sepset_elem = 13;
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
            NbrIdx[pos_of_new_sepset_elem] = G_Chunk[NbrIdxPointer];
            for (int d2 = 0; d2 < num_unfinished_sepsets; d2++)
            {
                YIdx = unfinished_prime[d2];
                if (YIdx == NbrIdx[pos_of_new_sepset_elem])
                {
                    continue;
                }

                for (int pos = 0; pos < pos_of_new_sepset_elem; pos++)
                {
                    NbrIdx[pos] = SepSet[(XIdx * n + YIdx) * ML + pos];
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

                rho = H[0][1] / (sqrt(fabs(H[0][0])) * sqrt(fabs(H[1][1])));
                Z = fabs(0.5 * (log(fabs((1 + rho))) - log(fabs(1 - rho))));
                pcorrs[(XIdx * n + YIdx) * PCORR_MAX_DEGREE + NbrIdxPointer] = Z;
            }
        }
    }
}