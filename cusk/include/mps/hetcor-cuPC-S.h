#ifndef CUPC_S_H_
#define CUPC_S_H_
//===============================> Definition <===============================
#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z
#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
// Note that EdgePerBlockLX * ParTestPerEdgeL1 < 1024

#define ParGivenL1 64
#define NumOfBlockForEachNodeL1 2
#define ParGivenL2 64
#define NumOfBlockForEachNodeL2 2
#define ParGivenL3 64
#define NumOfBlockForEachNodeL3 2
#define ParGivenL4 64
#define NumOfBlockForEachNodeL4 2
#define ParGivenL5 64
#define NumOfBlockForEachNodeL5 2
#define ParGivenL6 64
#define NumOfBlockForEachNodeL6 2
#define ParGivenL7 64
#define NumOfBlockForEachNodeL7 2
#define ParGivenL8 64
#define NumOfBlockForEachNodeL8 2
#define ParGivenL9 64
#define NumOfBlockForEachNodeL9 2
#define ParGivenL10 64
#define NumOfBlockForEachNodeL10 2
#define ParGivenL11 64
#define NumOfBlockForEachNodeL11 2
#define ParGivenL12 64
#define NumOfBlockForEachNodeL12 2
#define ParGivenL13 64
#define NumOfBlockForEachNodeL13 2
#define ParGivenL14 64
#define NumOfBlockForEachNodeL14 2
#define ML 14

#define PCORR_MAX_DEGREE 100

//==========================> Function Declaration <==========================
__global__ void cal_Indepl0(float *C, int *G, float *N, float th, int n);
__global__ void cal_Indepl1(
    float *C, int *G, float *N, int *GPrime, int *mutex, float th, int n
);
__global__ void cal_Indepl3(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl4(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl5(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl2(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl6(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl7(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl8(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl9(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl10(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl11(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl12(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl13(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__global__ void cal_Indepl14(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th
);
__device__ float mean_ess(float *N, int var_ixs[], int l, int n);

extern "C" void hetcor_skeleton(float *C, int *P, int *G, float *N, float *Th, int *l, const int *maxlevel);

#endif /* CUPC_S_H_ */
