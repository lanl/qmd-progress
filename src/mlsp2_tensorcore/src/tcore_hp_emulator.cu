#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <cuda_fp16.h>
#include <cublas.h>

// device function for splitting a float into two halves
__device__ void split_single(const float x, half &hi, half &lo)
{
    hi = __float2half(x);
    float y = (x - __half2float(hi));
    lo = __float2half(y * 1024);
}

// global function for splitting a float matrix into two float halves
template <typename T>
__global__ void array_split_single(const float *AF, half *AH1, half *AH2, const unsigned N)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N)
    {
        half hi;
        half lo;

        split_single(AF[i], hi, lo);

        AH1[i] = hi;
        AH2[i] = lo;
    }
}

void tcoreSPGemmSymm(cublasHandle_t handle,
                     const unsigned N,
                     const float *A,
                     half *Ah,
                     half *Al,
                     float *B1,
                     float *B2,
                     float *B)
{
    // Setup kernel launch
    unsigned num_thds = 512;
    unsigned num_blks = int(ceil(float(N * N) / float(num_thds)));

    // Split the floats into the high and low parts
    array_split_single<half><<<num_blks, num_thds>>>(A, Ah, Al, N * N);

    float alpha = 1.0;
    float beta = 0.0;

    // Compute gemmEx for high
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_T, N, N, N,
                 &alpha,
                 Ah, CUDA_R_16F, N,
                 Ah, CUDA_R_16F, N,
                 &beta, B1, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // Compute gemmEx for low
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_T, N, N, N,
                 &alpha,
                 Ah, CUDA_R_16F, N,
                 Al, CUDA_R_16F, N,
                 &beta, B2, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    alpha = 1.0;
    beta = 1.0;
    cublasSgeam(handle,
                CUBLAS_OP_N, CUBLAS_OP_T,
                N, N,
                &alpha,
                B2, N,
                &beta,
                B2, N,
                B, N);

    // undo prior scaling of 2^10
    beta = powf(2, -10);
    cublasSgeam(handle,
                CUBLAS_OP_N, CUBLAS_OP_N,
                N, N,
                &alpha,
                B1, N,
                &beta,
                B, N,
                B, N);
};

void tcoreHPGemmSymm(cublasHandle_t handle,
                     const unsigned N,
                     const float *A,
                     half *Ah,
                     half *Al,
                     float *B,
                     const float a)
{
    // Setup kernel launch
    unsigned num_thds = 512;
    unsigned num_blks = int(ceil(float(N * N) / float(num_thds)));

    // Split the floats into the high and low parts
    array_split_single<half><<<num_blks, num_thds>>>(A, Ah, Al, N * N);

    float zero = 0.0;

    // Compute gemmEx for high
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &a,
                 Ah, CUDA_R_16F, N,
                 Ah, CUDA_R_16F, N,
                 &zero, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
};

void tcoreSPGemmSymmAlt(cublasHandle_t handle,
                     const unsigned N,
                     const float *A,
                     half *Ah,
                     half *Al,
                     float *Buff,
                     float *B,
                     const float a)
{
    // Setup kernel launch
    unsigned num_thds = 512;
    unsigned num_blks = int(ceil(float(N * N) / float(num_thds)));

    // Split the floats into the high and low parts
    array_split_single<half><<<num_blks, num_thds>>>(A, Ah, Al, N * N);

    float one  = 1.0;
    float zero = 0.0;
    float undo = powf(2, -10); // undo prior scaling of 2^10

    // Compute gemmEx for low
    // Buff <- a ( Ah * Al )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &a,
                 Ah, CUDA_R_16F, N,
                 Al, CUDA_R_16F, N,
                 &zero, Buff, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B <- Buff + Buff.T = a( Ah * Al + Al * Ah )
    cublasSgeam(handle,
                CUBLAS_OP_N, CUBLAS_OP_T,
                N, N,
                &one,
                Buff, N,
                &one,
                Buff, N,
                B, N);

    // Compute gemmEx for high
    // B <- a ( Ah * Ah ) + undo * B = a( Ah * Ah + 2^-10 (Ah * Al + Al * Ah) )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &a,
                 Ah, CUDA_R_16F, N,
                 Ah, CUDA_R_16F, N,
                 &undo, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
};

// OVERWRITES A
// Gives B <- a * A^2 + b * A
void tcoreSPGemmSymmAlt2(cublasHandle_t handle,
                     const unsigned N,
                     float *A,
                     half *Ah,
                     half *Al,
                     float *B,
                     const float a,
                     const float b 
                    )
{
    // Setup kernel launch
    unsigned num_thds = 512;
    unsigned num_blks = int(ceil(float(N * N) / float(num_thds)));

    // Split the floats into the high and low parts
    array_split_single<half><<<num_blks, num_thds>>>(A, Ah, Al, N * N);

    float one  = 1.0;
    float halv = 0.5;
    float undo = powf(2, -10) * a; // undo prior scaling of 2^10

    // Compute gemmEx for high
    // A <- a (Ah * Ah) + b ( Ai )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &a,
                 Ah, CUDA_R_16F, N,
                 Ah, CUDA_R_16F, N,
                 &b, A, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // Compute gemmEx for low
    // A <- 2^-10 a ( Ah * Al ) + a/2 ( Ah * Ah ) + b/2 ( Ai )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo,
                 Ah, CUDA_R_16F, N,
                 Al, CUDA_R_16F, N,
                 &halv, A, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B <- A + A.T = a ( Ah * Ah + 2^-10 Ah * Al + 2^-10 Al * Ah ) + b ( Ai )
    cublasSgeam(handle,
                CUBLAS_OP_N, CUBLAS_OP_T,
                N, N,
                &one,
                A, N,
                &one,
                A, N,
                B, N);
};

// NOTE: B is directly replaced, A -> A2, C -> C2
void tcoreBlockSquare(cublasHandle_t handle,
                     const unsigned N,
                     float *A,
                     float *B,
                     float *C,
                     float *A2,
                     float *C2,
                     half *Ah,
                     half *Al,
                     half *Bh,
                     half *Bl,
                     half *Ch,
                     half *Cl,
                     const float a,
                     const float b){

    // Setup kernel launch
    unsigned num_thds = 512;
    unsigned num_blks = int(ceil(float(N * N) / float(num_thds)));

    // Split the floats into the high and low parts
    array_split_single<half><<<num_blks, num_thds>>>(A, Ah, Al, N * N);
    array_split_single<half><<<num_blks, num_thds>>>(B, Bh, Bl, N * N);
    array_split_single<half><<<num_blks, num_thds>>>(C, Ch, Cl, N * N);

    float one    = 1.0;
    float half_a = 0.5 * a;
    float half_b = 0.5 * b;
    float undo_a = powf(2, -10) * a; // undo prior scaling of 2^10


/*CUBLASAPI cublasStatus_t CUBLASWINAPI cublasGemmEx(cublasHandle_t handle,
                                                   cublasOperation_t transa,
                                                   cublasOperation_t transb,
                                                   int m,
                                                   int n,
                                                   int k,
                                                   const void* alpha,
                                                   const void* A,
                                                   cudaDataType Atype,
                                                   int lda,
                                                   const void* B,
                                                   cudaDataType Btype,
                                                   int ldb,
                                                   const void* beta,
                                                   void* C,
                                                   cudaDataType Ctype,
                                                   int ldc,
                                                   cublasComputeType_t computeType,
                                                   cublasGemmAlgo_t algo);
*/
    // UPDATE TO A: A2 = a( AA + BB' ) + bA 

    // Compute gemmEx for high
    // A <- a/2 (Ah * Ah) + b/2 ( Ai )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &half_a,
                 Ah, CUDA_R_16F, N,
                 Ah, CUDA_R_16F, N,
                 &half_b, A, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // Compute gemmEx for low
    // A <- 2^-10 a ( Ah * Al ) + a/2 ( Ah * Ah ) + b/2 ( Ai )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Ah, CUDA_R_16F, N,
                 Al, CUDA_R_16F, N,
                 &one, A, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // A <- 2^-10 a ( Ah * Al ) + a/2 ( Ah * Ah + Bh * Bh' ) + b/2 ( Ai )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_T, N, N, N,
                 &half_a,
                 Bh, CUDA_R_16F, N,
                 Bh, CUDA_R_16F, N,
                 &one, A, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
    
    // A <- 2^-10 a ( Ah * Al + Bh * Bl' ) + a/2 ( Ah * Ah + Bh * Bh' ) + b/2 ( Ai )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_T, N, N, N,
                 &undo_a,
                 Bh, CUDA_R_16F, N,
                 Bl, CUDA_R_16F, N,
                 &one, A, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // A2 <- 2^-10 a ( Ah * Al + Al * Ah + Bh * Bl' + Bl * Bh' ) + a ( Ah * Ah + Bh * Bh' ) + b ( Ai )
    cublasSgeam(handle,
                CUBLAS_OP_N, CUBLAS_OP_T,
                N, N,
                &one,
                A, N,
                &one,
                A, N,
                A2, N);


    // UPDATE TO B: B = a( AB + BC ) + bB

    // B <- a Ah * Bh + bB 
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &a,
                 Ah, CUDA_R_16F, N,
                 Bh, CUDA_R_16F, N,
                 &b, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B += a * 2^-10 Ah * Bl
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Ah, CUDA_R_16F, N,
                 Bl, CUDA_R_16F, N,
                 &one, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B += a * 2^-10 Al * Bh
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Al, CUDA_R_16F, N,
                 Bh, CUDA_R_16F, N,
                 &one, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B += a * Bh * Ch
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &a,
                 Bh, CUDA_R_16F, N,
                 Ch, CUDA_R_16F, N,
                 &one, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B += a * 2^-10 Bh * Cl
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Bh, CUDA_R_16F, N,
                 Cl, CUDA_R_16F, N,
                 &one, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // B += a * 2^-10 Bl * Ch
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Bl, CUDA_R_16F, N,
                 Ch, CUDA_R_16F, N,
                 &one, B, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
    

    // UPDATE TO C: C2 = a( B'B + CC ) + bC

    // Compute gemmEx for high
    // C <- a/2 (Ch * Ch) + b/2 ( Ci )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &half_a,
                 Ch, CUDA_R_16F, N,
                 Ch, CUDA_R_16F, N,
                 &half_b, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // Compute gemmEx for low
    // C <- 2^-10 a ( Ch * Cl ) + a/2 ( Ch * Ch ) + b/2 ( Ci )
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Ch, CUDA_R_16F, N,
                 Cl, CUDA_R_16F, N,
                 &one, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // C <- 2^-10 a ( Ch * Cl ) + a/2 ( Ch * Ch + Bh' * Bh ) + b/2 ( Ci )
    cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, N, N, N,
                 &half_a,
                 Bh, CUDA_R_16F, N,
                 Bh, CUDA_R_16F, N,
                 &one, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
    
    // C <- 2^-10 a ( Ch * Cl + Bh' * Bl ) + a/2 ( Ch * Ch + Bh' * Bh ) + b/2 ( Ci )
    cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, N, N, N,
                 &undo_a,
                 Bh, CUDA_R_16F, N,
                 Bl, CUDA_R_16F, N,
                 &one, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // C2 <- 2^-10 a ( Ch * Cl + Cl * Ch + Bh' * Bl + Bl' * Bh ) + a ( Ch * Ch + Bh' * Bh ) + b ( Ci )
    cublasSgeam(handle,
                CUBLAS_OP_N, CUBLAS_OP_T,
                N, N,
                &one,
                C, N,
                &one,
                C, N,
                C2, N);
    
}


void tcoreSPGemmAsymm(cublasHandle_t handle,
                     const unsigned N,
                     const float *A,
                     half *Ah,
                     half *Al,
                     const float *B,
                     half *Bh,
                     half *Bl,
                     float *C){

    // Setup kernel launch
    unsigned num_thds = 512;
    unsigned num_blks = int(ceil(float(N * N) / float(num_thds)));

    // Split the floats into the high and low parts
    array_split_single<half><<<num_blks, num_thds>>>(A, Ah, Al, N * N);
    array_split_single<half><<<num_blks, num_thds>>>(B, Bh, Bl, N * N);

    float one = 1.0;
    float zero = 0.0;
    float undo = powf(2, -10); // undo prior scaling of 2^10

    // C <- Ah * Bl
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &one,
                 Ah, CUDA_R_16F, N,
                 Bl, CUDA_R_16F, N,
                 &zero, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);

    // C <- Ah * Bl + Al * Bh
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &one,
                 Al, CUDA_R_16F, N,
                 Bh, CUDA_R_16F, N,
                 &one, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
    
    // C <- Ah * Bh + 2^-10 (Ah * Bl + Al * Bh)
    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N,
                 &one,
                 Ah, CUDA_R_16F, N,
                 Bh, CUDA_R_16F, N,
                 &undo, C, CUDA_R_32F, N,
                 CUBLAS_COMPUTE_32F_FAST_16F, CUBLAS_GEMM_DEFAULT);
    
};

void tcoreSPGemmSymm1(cublasHandle_t handle
                     ,const unsigned N
                     ,const float* A
                     ,const float* B
                     ,half*  Ah
                     ,half*  Al
                     ,half*  Bh
                     ,half*  Bl
                     ,float* C1
                     ,float* C2
                     ,float* C)
{
    // Setup kernel launch
    unsigned MAX_THREADS = 1024;
    unsigned BLOCKS = ceil(N*N/float(MAX_THREADS));
    unsigned THREADS = MAX_THREADS;

    // Split the floats into the high and low parts
    array_split_single<half><<<BLOCKS, THREADS>>>(A, Ah, Al, N*N);

    // Split the floats into the high and low parts
    array_split_single<half><<<BLOCKS, THREADS>>>(B, Bh, Bl, N*N);

    // Set the math mode to allow cuBLAS to use Tensor Cores:
    cublasStatus_t cublasStat = cublasSetMathMode(handle, CUBLAS_TENSOR_OP_MATH);

    float alpha (1.0f);
    float beta  (0.0f);
    float gamma = powf(2,-10);

    // Compute gemm for high
    cublasStat = cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha,
                              Ah, CUDA_R_16F, N,
                              Bh, CUDA_R_16F, N,
                              &beta, C1, CUDA_R_32F, N, CUDA_R_32F, CUBLAS_GEMM_DEFAULT_TENSOR_OP);

    // Compute gemms for low
    cublasStat = cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha,
                              Ah, CUDA_R_16F, N,
                              Bl, CUDA_R_16F, N,
                              &beta, C2, CUDA_R_32F, N, CUDA_R_32F, CUBLAS_GEMM_DEFAULT_TENSOR_OP);

    cublasStat = cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha,
                              Al, CUDA_R_16F, N,
                              Bh, CUDA_R_16F, N,
                              &alpha, C2, CUDA_R_32F, N, CUDA_R_32F, CUBLAS_GEMM_DEFAULT_TENSOR_OP);

    // add the high gemm and low gemm together
    cublasStat = cublasSgeam(handle,
                             CUBLAS_OP_N, CUBLAS_OP_N,
                             N, N,
                             &alpha,
                             C1, N,
                             &gamma,
                             C2, N,
                             C2, N);

    // compute C + C^T 
    cublasStat = cublasSgeam(handle,
                             CUBLAS_OP_N, CUBLAS_OP_T,
                             N, N,
                             &alpha,
                             C2, N,
                             &alpha,
                             C2, N,
                             C, N);

};
