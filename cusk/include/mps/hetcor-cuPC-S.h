__global__ void cal_Indepl0_ess(float *C, int *G, float *N, int n, float th);
__global__ void cal_Indepl1_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl3_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl4_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl5_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl2_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl6_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl7_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl8_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl9_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl10_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl11_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl12_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl13_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__global__ void cal_Indepl14_ess(
    float *C, int *G, float *N, int *GPrime, int *mutex, int n, float th, int *time_index
);
__device__ float mean_ess(float *N, int var_ixs[], int l, int n);
__device__ bool valid_time_conditioning(int a, int b, int *S, int *time_index, int l);
extern "C" void hetcor_skeleton(float *C, int *P, int *G, float *N, float *Th, int *l, const int *maxlevel, const int *time_index);