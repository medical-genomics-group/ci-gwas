#include <mps/cuPC-S-sparse.h>
#include <mps/cu_test_helpers.h>
#include <mps/gpuerrors.h>

void test_cal_Indepl0(float *C, int *M, int *P, int *W, int *G, double *Th)
{
    double *C_cuda;
    int *G_cuda;
    double *pMax_cuda;

    // num phen
    int p = *P;
    // num markers
    int m = *M;
    // corr width
    int w = *W;
    // height (nrows) of the corr matrix: num_markers + num_phen
    int nr = p + m;
    // width (ncols) of the corr matrix: corr_width + num_phen
    int nc = w + p;

    int max_marker_degree = 2 * w + p;
    int max_phen_degree = m + p;
    int mixed_matrix_size = max_marker_degree * m + max_phen_degree * p;

    HANDLE_ERROR(cudaMalloc((void **)&C_cuda, nc * nr * sizeof(double)));
    HANDLE_ERROR(cudaMalloc((void **)&G_cuda, mixed_matrix_size * sizeof(int)));
    HANDLE_ERROR(cudaMalloc((void **)&pMax_cuda, mixed_matrix_size * sizeof(int)));
    // copy correlation matrix from CPU to GPU
    HANDLE_ERROR(cudaMemcpy(C_cuda, C, nc * nr * sizeof(float), cudaMemcpyHostToDevice));
    CudaCheckError();

    dim3 BLOCKS_PER_GRID;
    dim3 THREADS_PER_BLOCK;

    if ((nc * nr) < 1024)
    {
        BLOCKS_PER_GRID = dim3(1, 1, 1);
        THREADS_PER_BLOCK = dim3(nr, nc, 1);
        cal_Indepl0<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(C_cuda, G_cuda, Th[0],
                                                            pMax_cuda, m, p, w);
        CudaCheckError();
    }
    else
    {
        BLOCKS_PER_GRID = dim3(ceil(((double)nr) / 32.0), ceil(((double)nc) / 32.0), 1);
        THREADS_PER_BLOCK = dim3(32, 32, 1);
        cal_Indepl0<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(C_cuda, G_cuda, Th[0],
                                                            pMax_cuda, m, p, w);
        CudaCheckError();
    }

    // Copy Graph G from GPU to CPU
    HANDLE_ERROR(cudaMemcpy(G, G_cuda, mixed_matrix_size * sizeof(int), cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaFree(C_cuda));
    HANDLE_ERROR(cudaFree(G_cuda));
    HANDLE_ERROR(cudaFree(pMax_cuda));
}
