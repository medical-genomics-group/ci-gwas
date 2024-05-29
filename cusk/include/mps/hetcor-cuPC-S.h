/*
 * cuPC_S.h
 *
 *  Created on: Apr 16, 2019
 *      Author: behrooz
 */

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
__device__ void pseudoinversel2(float M2[][2], float M2Inv[][2]);
__device__ void pseudoinversel3(float M2[][3], float M2Inv[][3]);
__device__ void pseudoinversel4(
    float M2[][4], float M2Inv[][4], float v[][4], float *rv1, float *w, float res1[][4]
);
__device__ void pseudoinversel5(
    float M2[][5], float M2Inv[][5], float v[][5], float *rv1, float *w, float res1[][5]
);
__device__ void pseudoinversel6(
    float M2[][6], float M2Inv[][6], float v[][6], float *rv1, float *w, float res1[][6]
);
__device__ void pseudoinversel7(
    float M2[][7], float M2Inv[][7], float v[][7], float *rv1, float *w, float res1[][7]
);
__device__ void pseudoinversel8(
    float M2[][8], float M2Inv[][8], float v[][8], float *rv1, float *w, float res1[][8]
);
__device__ void pseudoinversel9(
    float M2[][9], float M2Inv[][9], float v[][9], float *rv1, float *w, float res1[][9]
);
__device__ void pseudoinversel10(
    float M2[][10], float M2Inv[][10], float v[][10], float *rv1, float *w, float res1[][10]
);
__device__ void pseudoinversel11(
    float M2[][11], float M2Inv[][11], float v[][11], float *rv1, float *w, float res1[][11]
);
__device__ void pseudoinversel12(
    float M2[][12], float M2Inv[][12], float v[][12], float *rv1, float *w, float res1[][12]
);
__device__ void pseudoinversel13(
    float M2[][13], float M2Inv[][13], float v[][13], float *rv1, float *w, float res1[][13]
);
__device__ void pseudoinversel14(
    float M2[][14], float M2Inv[][14], float v[][14], float *rv1, float *w, float res1[][14]
);
__global__ void marginal_pMax(float *C, int *G, float th, float *pMax, int n);
__global__ void check_sepsets_l1(
    float *C, int *G, int *GPrime, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l2(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l3(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l4(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l5(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l6(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l7(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l8(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l9(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l10(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l11(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l12(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l13(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);
__global__ void check_sepsets_l14(
    float *C, int *G, int *GPrime, int *SepSet, float *pcorrs, int *unfinished_prime, int n
);

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
);

extern "C" void hetcor_skeleton(float *C, int *P, int *G, float *N, float *Th, int *l, const int *maxlevel);

__global__ void Initialize(int *Mat, int n);
__global__ void scan_compact(int *G_Compact, const int *G, const int n, int *nprime);
__global__ void SepSet_initialize(int *SepSet, int size);

__device__ void BINOM(int n, int k, int *out);
__device__ void IthCombination(int out[], int N, int P, int L);
__device__ float PYTHAG(float a, float b);
__device__ void inverse(float M2[][3], float M2Inv[][3]);

#endif /* CUPC_S_H_ */
