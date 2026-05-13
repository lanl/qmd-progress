#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include <cuda.h>
#include <cuda_fp16.h>
#include <cublas_v2.h>



__global__ void doubleToFloat(double *X,
                     float *Y,
                     int N)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N * N)
    {
        Y[i] = float(X[i]);
    }
}

__global__ void floatToDouble(float *X,
                     double *Y,
                     int N)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N * N)
    {
        Y[i] = double(X[i]);
    }
}

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

__global__ void shift_v1(float *X, int N, float shift)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N * N)
    {
        if (i % (N + 1) == 0)
        {
            X[i] += shift;
        }
    };
};

__global__ void shift_v2(float *X, int N, float shift)
{
    int i = ( threadIdx.x + blockIdx.x * blockDim.x ) * (N + 1);

    if (i < N * N)
    {
            X[i] += shift;
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

__global__ void add_with_shift(int N,
                    float *X, 
                    float *Y, 
                    const float a, 
                    const float b, 
                    const float c) 
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i < N * N)
    {
        if (i % (N + 1) == 0)
        {
            X[i] = a * X[i] + b * Y[i] + c;
        }
        else
        {
            X[i] = a * X[i] + b * Y[i];
        }

        // X[i] = a * X[i] + b * Y[i] + ((i % (N + 1) == 0) ? c : 0.0f);
    };
}

void gershgorin(const unsigned N,
                const double *X,
                double *h1,
                double *hN)
{
    float sum, diag_elem, minest, maxest;

    for (size_t i = 0; i < N; ++i)
    {
        sum = 0.0;
        minest = 0.0;
        diag_elem = 0.0;
        maxest = 0.0;
        for (size_t j = 0; j < N; ++j)
        {
            if (i != j)
            {
                sum += abs(X[i * N + j]); // assuming row major, running sum
            }
            else
            {
                diag_elem = X[i * N + i];
            }
        }

        minest = diag_elem - sum; // sum always non-neg
        maxest = diag_elem + sum; // sum always non-neg

        if (minest < h1[0])
        {
            h1[0] = minest;
        }
        if (hN[0] < maxest)
        {
            hN[0] = maxest;
        }
    }
};

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

/*
    Compute trace of matrix X on the device using openacc pragma.

    Inputs:
    -------
    TrX:     host pointer to trace, T*
    d_X:     device pointer to square matrix, T*
    N:       size of matrix, int
*/
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

/**
    Kernal for computing the trace of a matrix A of size N.
    Code originally from LATTE.
**/

__global__ void MatrixFastTraceKernel(const unsigned N, const int size, const float *A, float *result)
{

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    const int grid_size = blockDim.x * gridDim.x;

    // Create intermediate sums
    extern __shared__ float sdata[];
    sdata[threadIdx.x] = (float)0.0;
    __syncthreads();

    int i = tid;
    while (i < size)
    {
        sdata[threadIdx.x] += A[i + i * N];
        i += grid_size;
    }
    __syncthreads();

    // Reduce the values in shared memory
    int blockSize = blockDim.x;

    switch (blockSize)
    {
    case 1024:
        if (threadIdx.x < 512)
            sdata[threadIdx.x] += sdata[threadIdx.x + 512];
        __syncthreads();
    case 512:
        if (threadIdx.x < 256)
            sdata[threadIdx.x] += sdata[threadIdx.x + 256];
        __syncthreads();
    case 256:
        if (threadIdx.x < 128)
            sdata[threadIdx.x] += sdata[threadIdx.x + 128];
        __syncthreads();
    case 128:
        if (threadIdx.x < 64)
            sdata[threadIdx.x] += sdata[threadIdx.x + 64];
        __syncthreads();
        break;
    }

    if (threadIdx.x < 32)
    {

        volatile float *s_ptr = sdata;

        s_ptr[threadIdx.x] += s_ptr[threadIdx.x + 32];
        s_ptr[threadIdx.x] += s_ptr[threadIdx.x + 16];
        s_ptr[threadIdx.x] += s_ptr[threadIdx.x + 8];
        s_ptr[threadIdx.x] += s_ptr[threadIdx.x + 4];
        s_ptr[threadIdx.x] += s_ptr[threadIdx.x + 2];
        s_ptr[threadIdx.x] += s_ptr[threadIdx.x + 1];
    }

    // write result for this block to global mem
    if (threadIdx.x == 0)
        result[blockIdx.x] = sdata[0];
}

float M_Trace(const unsigned N,
              const float *A)
{

    // Set number of threads
    unsigned NUM_THREADS = 1024;

    // Size is N/2
    int size = N >> 1;
    int blockCount = (int)ceil((float)size / (float)NUM_THREADS);
    int smemSize = NUM_THREADS * sizeof(float);
    float *device_trace;
    float *local_trace = (float *)malloc(blockCount * sizeof(float));

    float trace = (float)0.0;

    cudaMalloc(&device_trace, blockCount * sizeof(float));

    cudaSetDevice(0);

    MatrixFastTraceKernel<<<blockCount, NUM_THREADS, smemSize>>>(N, N, A, device_trace);
    // MatrixFastTraceKernel<<<blockCount,NUM_THREADS,smemSize>>>(N, size, A, device_trace);

    // Copy to local variable
    // cudaThreadSynchronize();
    cudaMemcpy(local_trace, device_trace, blockCount * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(device_trace);

    for (int i = 0; i < blockCount; i++)
    {
        trace += local_trace[i];
    }

    free(local_trace);

    return trace;
}

/*
    Kernel for computing the trace (Tr) of a given matrix A
    of SIZE n (Will only work for even N currently)

    !!WARNING, CURRENT ALGO REQUIRES N TO BE EVEN

    Block level reduction with grid level atomic add
    This may not be the best way to do this but it is a
    good quick solution. Keep in mind the read is non-
    contiguous which is why I read it into block level
    shared memory first. -JSS
*/
__global__ void GPUtracekernel(const unsigned N, const float *A, float *Tr)
{
    // Get global thread id
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Externally declared block shared memory
    extern __shared__ float s_X[];

    bool inbound(i < N);

    // Set Tr to 0.0f with first thread
    if (i == 0)
    {
        Tr[0] = 0.0f;
        // printf("GPU A[0]: %f %f \n",A[0],A[2500]);
    }

    // Wait for thread in block
    __syncthreads();

    // Load non-contiguous memory into block level shared memory with initial add
    if (inbound)
    {
        // printf("value %i %i: %f + %f\n",i,threadIdx.x,A[i+2*N*i],A[(i+N)+2*N*(i+N)]);
        s_X[threadIdx.x] = A[i + 2 * N * i] + A[(i + N) + 2 * N * (i + N)];
    }
    __syncthreads();

    // Determine reduction size
    uint16_t M = (N - blockIdx.x * blockDim.x < blockDim.x) ? N - blockIdx.x * blockDim.x : blockDim.x;

    uint16_t Nred(static_cast<uint16_t>(floorf(static_cast<float>(M) / 2.0f)));
    // printf("%i\n",Nred);

    // Sync all threads in the block - non-"center" case
    while (Nred != 1)
    {
        __syncthreads();

        if (threadIdx.x < Nred && inbound)
        {
            if (M % 2 == 1 && threadIdx.x == Nred - 1)
            {
                s_X[threadIdx.x] += s_X[threadIdx.x + Nred] + s_X[threadIdx.x + Nred + 1];
            }
            else
            {
                // printf("value %i - %i - %i : %f + %f\n",i,threadIdx.x,threadIdx.x + Nred, s_X[threadIdx.x], s_X[threadIdx.x + Nred]);
                s_X[threadIdx.x] += s_X[threadIdx.x + Nred];
            }
        }

        M = Nred;
        Nred = static_cast<uint16_t>(floorf(static_cast<float>(Nred) / 2.0f));
    }

    __syncthreads();
    if (threadIdx.x == 0)
    {
        if (M == 3)
            atomicAdd(&Tr[0], s_X[0] + s_X[1] + s_X[2]);
        else
        {
            // printf("value %i: %f\n",i,s_X[0]+s_X[1]);
            atomicAdd(&Tr[0], s_X[0] + s_X[1]);
        }
    }
    __syncthreads();

    // if (i == 0)
    //     printf("GPU Tr[0]: %f \n",Tr[0]);
    // Tr[0] = 1.0;
};

/*
    Launcher for computing the trace (Tr) of a given matrix A
    of size N
*/
cudaError_t
GPUSTrace(const unsigned N,
          const float *A,
          float *Tr) // Assumed to be on the device
{
    // Setup kernel launch
    unsigned MAX_THREADS = 1024;
    unsigned BLOCKS = ceil((N / 2) / float(MAX_THREADS));
    unsigned THREADS = MAX_THREADS;

    // Call trace kernel
    GPUtracekernel<<<BLOCKS, THREADS, THREADS * sizeof(float), 0>>>(N / 2, A, Tr);

    return cudaPeekAtLastError();
};

__global__ void GPUtracekernel2(const unsigned N, const float *A, float *Tr)
{
    float local_sum = 0.0;
    for (unsigned i = 0; i < N; i++)
        local_sum += A[i + i * N];
    Tr[0] = local_sum;
};

__global__ void GPUtracekernel3(const unsigned N, const double *A, double *Tr)
{
    double local_sum = 0.0;
    for (unsigned i = 0; i < N; i++)
        local_sum += A[i + i * N];
    Tr[0] = local_sum;
};

/*
    Launcher for computing the trace (Tr) of a given matrix A
    of size N
*/
cudaError_t GPUSTrace2(const unsigned N, const float *A, float *Tr // Assumed to be on the device
                       ,
                       cudaStream_t cuStrm)
{
    // Setup kernel launch
    // unsigned MAX_THREADS = 1024;
    unsigned BLOCKS = 1;
    unsigned THREADS = 1;

    // Split the floats into the high and low parts
    GPUtracekernel2<<<BLOCKS, THREADS, 0, 0>>>(N, A, Tr);
    // std::cout << "ERROR: " << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;

    return cudaPeekAtLastError();
};

/*
    Launcher for computing the trace (Tr) of a given matrix A
    of size N
*/
cudaError_t GPUDTrace(const unsigned N, const double *A, double *Tr // Assumed to be on the device
                      ,
                      cudaStream_t cuStrm)
{
    // Setup kernel launch
    // unsigned MAX_THREADS = 1024;
    unsigned BLOCKS = 1;
    unsigned THREADS = 1;

    // Split the floats into the high and low parts
    GPUtracekernel3<<<BLOCKS, THREADS, 0, 0>>>(N, A, Tr);
    // std::cout << "ERROR: " << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;

    return cudaPeekAtLastError();
};

/*
    Kernel for computing C = alpha * A + beta * B for array A and B
    of size N
*/
__global__ void computeSnp1kernel(const unsigned N, const float *Sig, const float *A, const float *__restrict__ B, float *C)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)
    {
        C[i] = Sig[0] * A[i] + (1 - Sig[0]) * B[i];
    }
};

/*
    Launcher for computing C = alpha * A + beta * B for array A and B
    of size N
*/
cudaError_t computeSnp1(const unsigned N, const float *Sig, const float *A, const float *B, float *C // Assumed to be on the device
                        ,
                        cudaStream_t cuStrm)
{
    // Setup kernel launch
    unsigned MAX_THREADS = 1024;
    unsigned BLOCKS = ceil(N / float(MAX_THREADS));
    unsigned THREADS = MAX_THREADS;

    // Split the floats into the high and low parts
    computeSnp1kernel<<<BLOCKS, THREADS>>>(N, Sig, A, B, C);

    // Copy trace back to device (blocking)

    return cudaPeekAtLastError();
};

/*
    Kernel for computing Sigma
*/
__global__ void computeSigmakernel(unsigned Nocc, const float *TrXn, const float *TrX2n, float *Sig)
{
    if (fabs(TrX2n[0] - static_cast<float>(Nocc)) < fabs(2.0f * TrXn[0] - TrX2n[0] - static_cast<float>(Nocc)))
    {
        Sig[0] = 1;
    }
    else
    {
        Sig[0] = -1;
    }
};

/*
    Kernel launcher for computing Sigma
*/
cudaError_t computeSigma(unsigned Nocc,
                         const float *TrXn,
                         const float *TrX2n,
                         float *Sig)
{
    unsigned BLOCKS = 1;
    unsigned THREADS = 1;

    // Split the floats into the high and low parts
    computeSigmakernel<<<BLOCKS, THREADS>>>(Nocc, TrXn, TrX2n, Sig);

    return cudaPeekAtLastError();
};

/*
    Kernel for computing Sigma in double
*/
__global__ void computeSigmakernel_double(unsigned Nocc, const double *TrXn, const double *TrX2n, double *Sig)
{
    if (fabs(TrX2n[0] - static_cast<double>(Nocc)) < fabs(2.0f * TrXn[0] - TrX2n[0] - static_cast<double>(Nocc)))
    {
        Sig[0] = 1;
    }
    else
    {
        Sig[0] = -1;
    }
};

/*
    Kernel launcher for computing Sigma in double
*/
cudaError_t computeSigma_double(unsigned Nocc, const double *TrXn, const double *TrX2n, double *Sig, cudaStream_t cuStrm)
{
    unsigned BLOCKS = 1;
    unsigned THREADS = 1;

    // Split the floats into the high and low parts
    computeSigmakernel_double<<<BLOCKS, THREADS>>>(Nocc, TrXn, TrX2n, Sig);

    return cudaPeekAtLastError();
};

void doRefinement(double *_dA,
                  double *D0_,
                  const int _N,
                  const int _Nocc,
                  cublasHandle_t handle)
{

    double *d_T02, *d_T04;
    int N = _N;

    cudaMalloc(&d_T02, N * N * sizeof(double));
    // cudaMalloc(&d_T04,N*N*sizeof(double));

    // T0^2 in double precision
    double alpha_dbl = 1.0, beta_dbl = 0.0;

    cublasDgemm(handle,
                CUBLAS_OP_N, CUBLAS_OP_N,
                N, N, N,
                &alpha_dbl,
                _dA, N,
                _dA, N,
                &beta_dbl,
                d_T02, N);

    cudaMemcpy(D0_, d_T02, N * N * sizeof(double), cudaMemcpyDeviceToDevice);
    // cudaMemcpy(d_T04, d_T02, N * N * sizeof(double), cudaMemcpyDeviceToDevice);

    // 2*T0^2 - T0^4 in double precision
    alpha_dbl = -1.0, beta_dbl = 2.0;

    cublasDgemm(handle,
                CUBLAS_OP_N, CUBLAS_OP_N,
                N, N, N,
                &alpha_dbl,
                d_T02, N,
                d_T02, N,
                &beta_dbl,
                // d_T04, N);  // this function computes D = 2.0*T2 - 1.0*T2*T2 in double precision
                D0_, N); // this function computes D = 2.0*T2 - 1.0*T2*T2 in double precision

    // Move D to device and host
    // cudaMemcpy(D0_, d_T04, N * N * sizeof(double), cudaMemcpyDeviceToDevice);

    cudaFree(d_T02);
    // cudaFree(d_T04);

    // Cuda Error Checking
    // return cudaPeekAtLastError();
};
