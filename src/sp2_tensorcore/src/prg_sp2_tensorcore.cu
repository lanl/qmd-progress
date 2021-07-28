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
#include <cusolverDn.h>

#define  DIAG_OFF
#define  NO_REFINEMENT
#ifdef   SP2TCFortran

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
 * \param bndfil Percentage of occupied orbitals 
 * \param minsp2iter Minimium number of tensor core SP2 iterations 
 * \param maxsp2iter Maximum number of tensor core SP2 iterations 
 * \param sp2conv  
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
        cublasSetMathMode(handle, CUBLAS_TENSOR_OP_MATH);

    // Define grid size for __global__ function calls
    int numThreads = 512;
    int numBlocks = N * N / numThreads + 1; 
    
    // Declare Memory,
    double *d_T, *d_T2, *d_T4, *d_Idd;
    float  *sbuf1, *sbuf2, *TrS, *TrS2, *Sig, 
           *d_TrS2, *d_Sig, *d_TrS, *d_S, *d_S2, *d_Id;
    half   *hbuf1, *hbuf2;
    int    *v_sgn;

    cudaEvent_t start,stop,start_loop,stop_loop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventCreate(&start_loop);
    cudaEventCreate(&stop_loop);
    float elapsedTime, elapsedTime_loop, elapsedTime_rescale; 

    // Allocate Memory
    v_sgn = (int *) malloc(maxsp2iter * sizeof(int));
    cudaMalloc(&d_T, N * N * sizeof(double));
    cudaMalloc(&d_T2, N * N * sizeof(double));
    cudaMalloc(&d_T4, N * N * sizeof(double));
    cudaMalloc(&d_TrS, sizeof(float));
    cudaMalloc(&d_TrS2, sizeof(float));
    cudaMalloc(&d_Sig, sizeof(float));
    cudaMalloc(&d_S, N * N * sizeof(float));
    cudaMalloc(&d_S2, N * N * sizeof(float));
    cudaMalloc(&d_Id, N * N * sizeof(float));
    cudaMalloc(&d_Idd, N * N * sizeof(double));

    // Allocate cuda managed memory
    //cudaMallocManaged(&Id, N * N * sizeof(float));
    cudaMallocManaged(&TrS, sizeof(float));
    cudaMallocManaged(&TrS2, sizeof(float));
    cudaMallocManaged(&Sig, sizeof(float));
    
    // Allocate Buffers
    cudaMallocManaged(&sbuf1, N * N * sizeof(float));
    cudaMallocManaged(&sbuf2, N * N * sizeof(float));
    cudaMallocManaged(&hbuf1, N * N * sizeof(half));
    cudaMallocManaged(&hbuf2, N * N * sizeof(half));


    cudaEventRecord(start, 0);

    // Copy Hamiltonian to device
    cudaMemcpy(d_S, H, N * N * sizeof(float), cudaMemcpyHostToDevice);  
    cudaMemcpy(sbuf1, d_S, N * N * sizeof(float), cudaMemcpyDeviceToDevice);  
    
    
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    std::cout << "Time to transfer Hamiltonian = " << elapsedTime << " ms " << std::endl;
    
    


    
    //
    //
    // Estimate sprectral bounds
    //
    //
    cudaEventRecord(start, 0);
   
    float h1, hN;
    
    #ifdef DIAG_ON
    float *Eig;
    Eig = (float*) malloc(N * sizeof(float));
    linalgtools::getEigs(N, sbuf1, Eig);
    h1 = Eig[0]*1.01; 
    hN = Eig[N-1]*1.01;
    printf("h1 = %f \n", h1);
    printf("hN = %f \n", hN);
    #endif

    #ifdef DIAG_OFF  
    h1 =-27.547849409093747; //-27.04732512686311;//-27.229953288476242;
    hN = 25.0; //35.533175992743217; //52.378957263912767; //31.431533156948738;
    #endif


    float b = hN / (hN - h1);
    float a = -1 / (hN - h1);

    
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    std::cout << "Time for eigenvalues = " << elapsedTime << " ms " << std::endl;
    
    ////////////////////////////////////////
    
    
    //
    //
    // Begin DNN-SP2
    //
    //

    cudaEventRecord(start_loop, 0);    

    #ifdef SP2_SINGLE
    float alphaS = 1.0, betaS = 0.0;
    #endif
    
    // Build idenity on GPU
    build_identity_gpu<<< numBlocks, numThreads >>>(d_Id, N);
    
    // Rescale and fold eigenspectrum within unit interval (0,1)
    // compute initial layer of the DNN, S0=W*H+B
    cublasStat = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                            N, N, N, 
                            &b,
                            d_Id, N, 
                            d_Id, N, 
                            &a, 
                            d_S, N);  

    // Compute initial trace 
    linalgtools::GPUSTrace(N, d_S, d_TrS);
    cudaMemcpy(TrS, d_TrS, sizeof(float), cudaMemcpyDeviceToHost);

   
    // SP2 DNN Loop
    while (Stopp == 0)
    {
       
        //S^2=X
        tcoretools::tcoreSPGemmSymm(handle, N, 
                                   d_S, 
                                   hbuf1, 
                                   hbuf2, 
                                   sbuf1, 
                                   sbuf2,
                                   d_S2);
/*        //num+=1;
        //std::cout << num << std::endl;
   
        //S_high^2=X
        tcoretools::tcoreSPGemmSymm_ZEROX1(handle, N, 
                                   d_S, 
                                   hbuf1, 
                                   d_S2);
*/
        
        #ifdef SP2_SINGLE
        // full single-precision S^2
        cublasStat = cublasSgemm(handle,
                             CUBLAS_OP_N, CUBLAS_OP_N,
                             N, N, N,
                             &alphaS,
                             d_S, N,
                             d_S, N,
                             &betaS,
                             d_S2, N);
        #endif

        // Trace of S^2
        linalgtools::GPUSTrace(N, d_S2, d_TrS2);
        cudaMemcpy(TrS2, d_TrS2, sizeof(float), cudaMemcpyDeviceToHost);
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
        linalgtools::computeSigma(Nocc, d_TrS, d_TrS2, d_Sig);

        // Compute S_{n+1}
        linalgtools::computeSnp1(N * N, d_Sig, d_S2, d_S, d_S);
        
        // Copy trace and sigma data to CPU
        cudaMemcpy(TrS2, d_TrS2, sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(Sig, d_Sig, sizeof(float), cudaMemcpyDeviceToHost);
        
        // Compute TrS
        TrS[0] = Sig[0] * TrS2[0] + (1 - Sig[0]) * TrS[0];
        cudaMemcpy(d_TrS, TrS, sizeof(float), cudaMemcpyHostToDevice);
        
        // Update sign vector
        v_sgn[iter] = int(Sig[0]);

        iter += 1;
    }

    cudaEventRecord(stop_loop, 0);
    cudaEventSynchronize(stop_loop);
    cudaEventElapsedTime(&elapsedTime_loop, start_loop, stop_loop);
    

    //
    //
    // Print loop time and tflops
    //
    //

    double TFLOPS = 2*double(N)*double(N)*double(N)*(iter+1.5)/ \
                    ((elapsedTime_loop)/double(1e3))/double(1e12); 

    std::cout << "Time for SP2 loop = " << elapsedTime_loop << " ms " << std::endl;
    std::cout << "Loop FLOP rate = " << TFLOPS << " TFLOPS" <<std::endl;

    printf("Number of iters = %d\n", iter);

    ///////////////////////////////////////////////





    //
    //
    // End DNN-SP2
    //
    //

    cudaEventRecord(start, 0);
   
    #ifdef NO_REFINEMENT

    FtoD<<<numBlocks,numThreads>>>(d_S, d_T, N);
    cudaMemcpy(D, d_T, N * N * sizeof(double), cudaMemcpyDeviceToHost);
    std::cout << "NO REFINEMENT" << std::endl; 
    
    #endif  

    #ifdef REFINEMENT
     
    ///////////////////////////////////////////////////////
    ///////// compute refinement step via GPU ////////////
    //////////////////////////////////////////////////////
   
    // Convert matrices to double prec on device
    FtoD<<<numBlocks,numThreads>>>(d_S, d_T, N);
    FtoD<<<numBlocks,numThreads>>>(d_Id, d_Idd, N);

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
    
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    std::cout << "Finalize and transfer density matrix = " << elapsedTime << " ms " << std::endl;
  
    #endif

    //Deallocations
    free(v_sgn);
    
    cudaFree(TrS);
    cudaFree(TrS2);
    cudaFree(Sig);
    cudaFree(d_T);
    cudaFree(d_T2);
    cudaFree(d_T4);
    cudaFree(d_Id);
    cudaFree(d_Idd);
    cudaFree(d_S);
    cudaFree(d_S2);

    // Buffers
    cudaFree(sbuf1);
    cudaFree(sbuf2);
    cudaFree(hbuf1);
    cudaFree(hbuf2);
    
    // Destroy handle
    cublasDestroy(handle);


}
