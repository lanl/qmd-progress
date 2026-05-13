#include <cublas_v2.h>


__global__ void doubleToFloat(double*,
                     float*,
                     int);

__global__ void floatToDouble(float*,
                     double*,
                     int);

template<typename T>
__global__ void setToIdentityMatrix(T *X, int N)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N * N)
    {
        if (i % (N + 1) == 0)
        {
            X[i] = 1.0;
        }
        else
        {
            X[i] = 0.0;
        }
    };
};

template<typename T>
__global__ void shift_and_scale(int N,
                    T *X, 
                    const T shift,
                    const T scale)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N * N)
    {
        X[i] = scale * X[i] + ((i % (N + 1) == 0) ? shift : 0.0);
    };
}

__global__ void shift_v1(float*, int, float);

__global__ void shift_v2(float*, int, float);

__global__ void add_with_shift(int,
                                float*, 
                                float*, 
                                const float, 
                                const float, 
                                const float); 

void gershgorin(const unsigned,
                const double *,
                double *,
                double *);

//template <typename T> 
//extern void gershgorin_v2(const unsigned,
//                          const T *,
//                          T *,
//                          T *);
template<typename T> 
void gershgorin_v2(const unsigned N,
                   const T *X,
                   T *h1,
                   T *hN)
{
    T minest, maxest;

    minest = 0.0;
    maxest = 0.0;

    #pragma acc parallel loop reduction(min : minest) reduction(max : maxest) deviceptr(X)
    for (size_t i = 0; i < N; ++i)
    {
        float sum = 0;
        #pragma acc loop reduction(+ : sum)
        for (size_t j = 0; j < N; ++j)
        {
            sum += abs(X[i * N + j]); // assuming row major, running sum
        }

        minest = min(2 * X[i * N + i] - sum, minest); // sum always non-neg
        maxest = sum;                                 // sum always non-neg
    }
    h1[0] = minest;
    hN[0] = maxest;
};

template<typename T> 
void openacc_trace(T* TrX, T* d_X, const int N)
{
    // Compute initial traces
    T trace = 0.0;

    #pragma acc parallel loop deviceptr(d_X) reduction(+:trace)
    for(int i=0; i < N; i++){
        trace += d_X[i*N + i];
    }
    TrX[0] = trace;
}

template<typename T> 
void openacc_trace2(T* TrX2, T* d_X, const int N){
    T trace2 = 0.0;

    #pragma acc parallel loop deviceptr(d_X) reduction(+:trace2)
    for(int i=0; i<N*N; i++){
        T x = d_X[i];
        trace2 += x * x;
    }
    TrX2[0] = trace2;
}

cudaError_t GPUSTrace(const unsigned,
                      const float *,
                      float *); // Assumed to be on the device

cudaError_t GPUDTrace(const unsigned,
                      const double *,
                      double *); // Assumed to be on the device

cudaError_t GPUSTrace2(const unsigned,
                       const float *,
                       float *B); // Assumed to be on the device

cudaError_t computeSnp1(const unsigned,
                        const float *,
                        const float *,
                        const float *,
                        float *); // Assumed to be on the device

void computeSigma(unsigned Nocc, const float *TrXn, const float *TrX2n, float *Sig);

cudaError_t doRefinement(double *, double *, const int, const int, cublasHandle_t);
