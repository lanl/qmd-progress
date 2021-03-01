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
    float h1 = -30.0;//-27.229953288476242;
    float hN = 55.0; //31.431533156948738;

    std::vector < float >Idemp_Error;

    // Set GPU
    int device = 0;
    cudaSetDevice(device);

    // Cublas Handle
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasStatus_t cublasStat =
        cublasSetMathMode(handle, CUBLAS_TENSOR_OP_MATH);

    // Declare Memory,
    double *T, *d_T, *d_T2, *d_T4, *TrT, *TrT2, *Idd, *d_Idd;
    float *S, *S2, *Id, *sbuf1, *sbuf2, *TrS, *TrSOld, 
        *TrS2, *Sig;
    half *hbuf1, *hbuf2;
    int *v_sgn;

    float b = hN / (hN - h1);
    float a = -1 / (hN - h1);

    // Allocate Memory
    cudaMallocManaged(&S, N * N * sizeof(float));
    cudaMallocManaged(&S2, N * N * sizeof(float));
    v_sgn = (int *) malloc(N * sizeof(int));
    T = (double *) malloc(N * N * sizeof(double));
    cudaMalloc(&d_T, N * N * sizeof(double));
    cudaMalloc(&d_T2, N * N * sizeof(double));
    cudaMalloc(&d_T4, N * N * sizeof(double));
    Idd = (double *) malloc(N * N * sizeof(double));
    cudaMalloc(&d_Idd, N * N * sizeof(double));

    // Allocate cuda managed memory
    cudaMallocManaged(&Id, N * N * sizeof(float));
    cudaMallocManaged(&TrS, sizeof(float));
    cudaMallocManaged(&TrS2, sizeof(float));
    cudaMallocManaged(&TrT, sizeof(double));
    cudaMallocManaged(&TrT2, sizeof(double));
    cudaMallocManaged(&TrSOld, sizeof(float));
    cudaMallocManaged(&Sig, sizeof(float));
    
    // Allocate Buffers
    cudaMallocManaged(&sbuf1, N * N * sizeof(float));
    cudaMallocManaged(&sbuf2, N * N * sizeof(float));
    cudaMallocManaged(&hbuf1, N * N * sizeof(half));
    cudaMallocManaged(&hbuf2, N * N * sizeof(half));

    // Produce Hamiltonian and Identity matrix
    cudaMemcpy(S, H, N * N * sizeof(float), cudaMemcpyHostToDevice);    // Send H to S

    build_identity(N, Id);
    CPU_float_to_double(Id, D, N);
    
    // Get device id
    cudaGetDevice(&device);

    //compute initial layer of the DNN, W*S+B

    cublasStat = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                            N, N, N, 
                            &b,
                            Id, N, 
                            Id, N, 
                            &a, 
                            S, N);    //this function computes S = b*Id*Id + a*S = W*S + B


    // Compute initial trace
    linalgtools::GPUSTrace2(N, S, TrS);

    // SP2 DNN Loop
    while (Stopp == 0)
    {
        cudaMemPrefetchAsync(TrS, sizeof(float), device, NULL);
        cudaDeviceSynchronize();

        //S^2
        tcoretools::tcoreSPGemmSymm(handle, N, 
                                   S, 
                                   hbuf1, 
                                   hbuf2, 
                                   sbuf1, 
                                   sbuf2,
                                   S2);


        // Trace of S^2
        linalgtools::GPUSTrace2(N, S2, TrS2);
        cudaMemPrefetchAsync(TrS2, sizeof(float), device, NULL);
        cudaDeviceSynchronize();
        Idemp_Error.push_back((TrS[0] - TrS2[0]));

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
    //exit(0);

    //////////////////////////////////////////////////////
    ////////////// Refinement starts here ////////////////
    //////////////////////////////////////////////////////
    CPU_float_to_double(S, D, N);  // NO REFINEMENT
    std::cout << "NO REFINEMENT" << std::endl; 
    //CPU_float_to_double(S, T, N);
    //CPU_float_to_double(Id, Idd, N);
    //////////////////////////////////////////////////////
/*
    //////////////////////////////////////////////////////
    //// Send double precision object back to the GPU ///
    //////////////////////////////////////////////////////
    cudaMemcpy(d_T, T, N * N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Idd, Idd, N * N * sizeof(double), cudaMemcpyHostToDevice);
    //////////////////////////////////////////////////////

    // Compute T^2 in double prec since last update was only to S, not S^2
    double alpha_dbl = 1.0, beta_dbl = 0.0;
    cublasStat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 
                            N, N, N, 
                            &alpha_dbl, 
                            d_T, N, 
                            d_T, N, 
                            &beta_dbl, 
                            d_T2, N);        // this function computes T2 = alpha_dbl*T*T + beta_dbl*T2 = T^2 
                                             // in double precision
    cudaDeviceSynchronize();
    cudaMemcpy(d_T4, d_T2, N * N * sizeof(double), cudaMemcpyDeviceToDevice);


    //////////////////////////////////////////////////////
    ////////////// compute matrix D via GPU //////////////
    //////////////////////////////////////////////////////
    alpha_dbl = -1.0, beta_dbl = 2.0;
    cublasStat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 
                            N, N, N, 
                            &alpha_dbl, 
                            d_T2, N, 
                            d_T2, N, 
                            &beta_dbl, 
                            d_T4, N);      // this function computes D = 2.0*T2 - 1.0*T2*T2 in double precision
    cudaMemcpy(D, d_T4, N * N * sizeof(double), cudaMemcpyDeviceToHost);
    //////////////////////////////////////////////////////
*/
    //Deallocations
    cudaFree(S);
    cudaFree(S2);
    free(v_sgn);
    free(T);
    cudaFree(d_T);
    cudaFree(d_T2);
    cudaFree(d_T4);
    free(Idd);
    cudaFree(d_Idd);
    
    // cuda managed memory
    cudaFree(Id);
    cudaFree(TrS);
    cudaFree(TrS2);
    cudaFree(TrT);
    cudaFree(TrT2);
    cudaFree(TrSOld);
    cudaFree(Sig);

    // Buffers
    cudaFree(sbuf1);
    cudaFree(sbuf2);
    cudaFree(hbuf1);
    cudaFree(hbuf2);
    
    // Destroy handle
    cublasDestroy(handle);


}
