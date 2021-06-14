#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <regex>
#include <typeinfo>
#include <cuda.h>
#include <cublas_v2.h>
#include <cuda_fp16.h>
#include <random>
#include <cuda_runtime.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <chrono>
#include "tcore_hp_emulator.cuh"
#include "linalg_tools.cuh"

#ifdef  SP2TCFortran
extern "C"
{
    void prg_sp2_tensorcore(
    int,
    float *,
    double *,
    float,
    float,
    int,
    int,
    char,
    float,
    int);
}
#endif
void
CPU_float_to_double(
    float *S,
    double *T,
    int N)
{
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            T[i + j * N] = double (
    S[i + j * N]);
        }
    }
};


__global__ 
void 
FtoD(
    float *X, 
    double *Y, 
    int N) 
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  while (i < N * N) {
    Y[i] = double(X[i]);
    i += blockDim.x * gridDim.x; // add total number of threads to i
  }
}


__global__ 
void 
build_identity_gpu(
     float* X,
     int N)
{  
  int i = threadIdx.x + blockIdx.x * blockDim.x; 
  
  while (i < N * N) {
    if ( i % (N+1) == 0) {
      X[i] = 1.0f;
    } 
    else {
      X[i] = 0.0f;
    }
    i += blockDim.x * gridDim.x;  // add total number of threads to i
}
}

void
build_identity(
    const unsigned N,
    float *X)
{

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (i == j)
            {
                X[i * N + j] = 1;
            }
            else
            {
                X[i * N + j] = 0;
            }
        }
    }
};


/** SP2 Tensor core routine.
 * 
 *
 *  Returns the density matrix as computed using the Tensor 
 *  core version of SP2.
 * 
 *
 * \param N Number of orbitals (Size of the Hamiltonian)
 * \param H Pointer to the Hamiltonian array.
 * \param D Pointer to the Density matrix array. 
 * \param eps 
 */
void
prg_sp2_tensorcore(
    int N,
    float *H,
    double *D,
    float eps,
    float bndfil,
    int minsp2iter,
    int maxsp2iter,
    char sp2conv,
    float idemtol,
    int verbose)
{

    // Matrix size
    int Nocc = int(bndfil * N);

    int Stopp = 0;
    int iter = 0;

    // Prior estimate for spectral bounds
    float h1 =-27.547849409093747; //-27.04732512686311;//-27.229953288476242;
    float hN = 35.533175992743217; //52.378957263912767; //31.431533156948738;

    std::vector < float >Idemp_Error;

    // Set GPU
    int device = 0;
    cudaSetDevice(device);

    // Get device id
    // cudaGetDevice(&device);
    
    // Cublas Handle
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasStatus_t cublasStat =
        cublasSetMathMode(handle, CUBLAS_DEFAULT_MATH);

    // Define grid size for __global__ function calls
    int numThreads = 128;
    int numBlocks = N * N / 80 / 128 + 1; 
    
    // Declare Memory,
    double *d_T, *d_T2, *d_T4, *TrT, *TrT2, *d_Idd;
    float *S, *S2, *Id, *sbuf1, *sbuf2, *TrS, *TrSOld, 
        *TrS2, *Sig;
    half *hbuf1, *hbuf2;
    int *v_sgn;

    float b = hN / (hN - h1);
    float a = -1 / (hN - h1);

    // Allocate Memory
    v_sgn = (int *) malloc(N * sizeof(int));
    cudaMalloc(&d_T, N * N * sizeof(double));
    cudaMalloc(&d_T2, N * N * sizeof(double));
    cudaMalloc(&d_T4, N * N * sizeof(double));

    // Allocate cuda managed memory
    cudaMallocManaged(&S, N * N * sizeof(float));
    cudaMallocManaged(&S2, N * N * sizeof(float));
    cudaMallocManaged(&Id, N * N * sizeof(float));
    cudaMallocManaged(&TrS, sizeof(float));
    cudaMallocManaged(&TrS2, sizeof(float));
    cudaMallocManaged(&TrT, sizeof(double));
    cudaMallocManaged(&TrT2, sizeof(double));
    cudaMallocManaged(&TrSOld, sizeof(float));
    cudaMallocManaged(&Sig, sizeof(float));
    cudaMallocManaged(&d_Idd, N * N * sizeof(double));
    
    // Allocate Buffers
    cudaMallocManaged(&sbuf1, N * N * sizeof(float));
    cudaMallocManaged(&sbuf2, N * N * sizeof(float));
    cudaMallocManaged(&hbuf1, N * N * sizeof(half));
    cudaMallocManaged(&hbuf2, N * N * sizeof(half));

    // Copy Hamiltonian to device
    cudaMemcpy(S, H, N * N * sizeof(float), cudaMemcpyHostToDevice);  

    // Build idenity on GPU
    //build_identity(N, Id);
    build_identity_gpu<<< numBlocks, numThreads >>>(Id, N);
    

    //cudaDeviceSynchronize();
    
    // Rescale and fold eigenspectrum within unit interval (0,1)
    // compute initial layer of the DNN, S0=W*H+B
    cublasStat = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                            N, N, N, 
                            &b,
                            Id, N, 
                            Id, N, 
                            &a, 
                            S, N);  

    // Compute initial trace -- further optimization possible
    linalgtools::GPUSTrace2(N, S, TrS);

    float alphaS = 1.0, betaS = 0.0;

    // SP2 DNN Loop
    while (Stopp == 0)
    {
        cudaMemPrefetchAsync(TrS, sizeof(float), device, NULL);
        cudaDeviceSynchronize();

        //S^2=X
        /*tcoretools::tcoreSPGemmSymm(handle, N, 
                                   S, 
                                   hbuf1, 
                                   hbuf2, 
                                   sbuf1, 
                                   sbuf2,
                                   S2);
        */
        /* //S_high^2=X
        tcoretools::tcoreSPGemmSymm_ZEROX1(handle, N, 
                                   S, 
                                   hbuf1, 
                                   S2);*/

        // full single-precision S^2
        cublasStat = cublasSgemm(handle,
                             CUBLAS_OP_N, CUBLAS_OP_N,
                             N, N, N,
                             &alphaS,
                             S, N,
                             S, N,
                             &betaS,
                             S2, N);

 
        // Trace of S^2
        linalgtools::GPUSTrace2(N, S2, TrS2);
        cudaMemPrefetchAsync(TrS2, sizeof(float), device, NULL);
        cudaDeviceSynchronize();
        Idemp_Error.push_back((TrS[0] - TrS2[0]));
        //printf("%f\n", TrS[0] - TrS2[0]);
 
        // Convergence Control
        if (TrS[0] - TrS2[0] <= 0)
        {
            break;
        };
        if (iter > 2 && v_sgn[iter - 1] != v_sgn[iter - 2]
            && Idemp_Error[iter] >=
            4.5 * Idemp_Error[iter - 2] * Idemp_Error[iter - 2])
        {
            break;
        };

        // Compute Sigma
        linalgtools::computeSigma(Nocc, TrS, TrS2, Sig);

        // Compute S_{n+1}
        linalgtools::computeSnp1(N * N, Sig, S2, S, S);
        cudaDeviceSynchronize();

        // Compute TrS
        TrS[0] = Sig[0] * TrS2[0] + (1 - Sig[0]) * TrS[0];
        
        // Update sign vector
        v_sgn[iter] = int(Sig[0]);

        cudaMemPrefetchAsync(TrS, sizeof(float), device, NULL);
        iter += 1;

    }
    cudaDeviceSynchronize();


/*
    // Uncomment to turn offf refinement
    FtoD<<<numBlocks,numThreads>>>(S, d_T, N);
    cudaMemcpy(D, d_T, N * N * sizeof(double), cudaMemcpyDeviceToHost);
    std::cout << "NO REFINEMENT" << std::endl; 
  */

    ///////////////////////////////////////////////////////
    ///////// compute refinement step via GPU ////////////
    //////////////////////////////////////////////////////
   
    // Convert matrices to double prec on device
    FtoD<<<numBlocks,numThreads>>>(S, d_T, N);
    FtoD<<<numBlocks,numThreads>>>(Id, d_Idd, N);

    // Compute T^2 in double prec
    double alpha_dbl = 1.0, beta_dbl = 0.0;
    cublasStat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 
                            N, N, N, 
                            &alpha_dbl, 
                            d_T, N, 
                            d_T, N, 
                            &beta_dbl, 
                            d_T2, N);      // this function computes T^2
                                          
    cudaMemcpy(d_T4, d_T2, N * N * sizeof(double), cudaMemcpyDeviceToDevice);


    alpha_dbl = -1.0, beta_dbl = 2.0;
    cublasStat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 
                            N, N, N, 
                            &alpha_dbl, 
                            d_T2, N, 
                            d_T2, N, 
                            &beta_dbl, 
                            d_T4, N);      // this function computes D = 2*T^2 - T^4
    
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    
    
    // Send purified density matrix estimate to output D
     cudaMemcpy(D, d_T4, N * N * sizeof(double), cudaMemcpyDeviceToHost);
 
    //Deallocations
    free(v_sgn);
    
    cudaFree(S);
    cudaFree(S2);
    cudaFree(Id);
    cudaFree(TrS);
    cudaFree(TrS2);
    cudaFree(TrT);
    cudaFree(TrT2);
    cudaFree(TrSOld);
    cudaFree(Sig);
    cudaFree(d_T);
    cudaFree(d_T2);
    cudaFree(d_T4);
    cudaFree(d_Idd);

    // Buffers
    cudaFree(sbuf1);
    cudaFree(sbuf2);
    cudaFree(hbuf1);
    cudaFree(hbuf2);
    
    // Destroy handle
    cublasDestroy(handle);


}
