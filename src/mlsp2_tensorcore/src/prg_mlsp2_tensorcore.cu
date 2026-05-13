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
#include <error_check.cuh>
#include "nvToolsExt.h"

#include <mlsp2.cuh>
#include <structs.h>

#define  DIAG_OFF
#define  NO_REFINEMENT

extern "C"
{
    void prg_mlsp2_tensorcore(
    int,
    double *,
    double *,
    double,
    double,
    int);
}

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
 * \param beta 
 * \param mu
 * \param verbose
 */
void
prg_mlsp2_tensorcore(
    int N,
    double *H,
    double *D,
    double beta,
    double mu,
    int verbose)
{

    std::cout <<  H[0] << std::endl;
    mlsp2_full(H, D, N, beta, mu, fp16_fp32, yes);

}


double model26_no_c[4 * 26] = {
-1.193067993736704, 1.2808163708951992, -2.0138064431153593, 1.737960453274357, -2.2432847067890003, -1.532607107366334, 1.3706555030686003, 2.0866943361854076, -4.226242242971189, -2.1423110498624207, 1.531862467247144, 1.9247548697829673, -3.9195673016046926, -2.0694801324608476, 1.4592758562389494, -1.801859259883661, 1.2248811540799969, 1.3964770115785683, -1.850619031018666, 1.2360916969650926, -1.3001631153274709, -1.608974630365624, 1.803207899542053, -1.200760280811683, 0.5468330861974116, 0.3488684055678595,
2.3804416902423853, -0.5900164628932307, 2.62992913023114, -0.25142304945791194, 2.9047520517320713, 2.2812744371019855, -0.644086018749064, 0.06301322491382681, 2.9830683285944253, 2.1962521358450577, -0.06493497205130346, -0.032520353452927475, 2.834189036901693, 2.1586513285414006, 0.03839243974588228, 1.9997711274733334, 0.21609368766628487, 0.21865194215881015, 1.9278961170374873, 0.23022056726418688, 1.420053950206337, 1.7964133930938535, 0.5526179671661151, 1.788010071849649, 0.4132593140354784, 0.43669028900306267,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
5.96068593395463e-12, -3.617270023664845e-11, 1.3051116284391382e-11, 2.5592521843964343e-11, -1.677066815136476e-11, -1.2620290229094606e-11, 2.9889491881707607e-12, 8.817716613316653e-11, -6.585311617630012e-11, 9.91465343479627e-12, 2.6084430682798928e-11, 1.0845293445046395e-06, 2.1348214277971534e-07, -0.00028676968523134295, -5.488307622373344e-05, -0.008480064600676285, -0.009231362135425173, -0.15112292817822895, -0.15026466549692039, -0.06877291057052441, -0.16447579024170086, -0.7699487281424027, 0.6949647639348486, 1.3934328465214243, -0.17174843386097668, 0.4366902890002094
};

double model26[4 * 26] = {
-1.353996666080607, 1.608677030263155, -2.3516367358921175, 1.810031371507046, -2.732204264661287, -1.8739874213782346, 1.4022610794248502, 1.7130894007631559, -3.252595736449542, -1.878046620274372, 1.3104640080496117, 1.5636089958977137, -3.0389837384371026, -1.6993773295730852, 1.2569406703914685, -1.3902819892914342, 1.1724154920520065, 1.3514480459756721, -1.212734900653146, 0.8610663421144649, -1.3691733390623337, -1.452079661340262, 1.3472629229112754, -1.1142383533588098, 0.9411010738180853, 0.8321693583250559,
2.437540893737664, -0.2987544634415441, 2.61381390888074, -0.037330360494104514, 2.8258863420021245, 2.1387983184086994, 0.04890736627204146, -0.1034037621506357, 2.6286496267174138, 2.0748803496087422, 0.05496499701340592, -0.008606990671650539, 2.4776868627993966, 2.0765652847572458, 0.05777573303514125, 1.9662437062782698, 0.02982302801604828, 0.057730643800221824, 1.861660915086903, 0.010487472035461482, 1.6990777043351213, 1.993171071701545, 0.02198074427424086, 1.8969376233131163, -0.1554719335561439, 0.5698259888759222,
-0.20384668041719653, -0.06746741346599809, -0.0020811153244366387, -0.05442179841399564, 0.008167975559885774, 0.03573114610928145, -0.04237220339751432, -0.024380949423571694, 0.03165347918592593, 0.03449163838903798, 0.003310087064751345, -0.013010654680499112, 0.026368712019503367, 0.0344246980219247, 0.049025969756312696, 0.019105034898269455, 0.04530045685625948, 0.03753306163777001, 0.01242573869144004, 0.04701135000272372, -0.007732969926664974, -0.016611158872197006, 0.19956152561687585, -0.07280733281872145, 0.0658425383539312, -0.42348060531994464,
1.2911561475361923e-11, -7.006084032292678e-12, -6.572840812953572e-13, 4.3185035417757944e-11, -1.6624560557980038e-11, 4.542623791264443e-12, 2.4246485246358385e-11, 7.873381643056254e-11, -6.312379678783267e-11, -6.649574392039178e-11, 2.825980776677432e-10, -1.424083924407828e-07, -1.3301005714900714e-07, 0.00489644656369143, -0.0023020549125814128, 0.016062593008344738, -0.015390213336897312, -0.12767959311255042, -0.1100276149057915, -0.0494684380432657, -0.17612285306119255, -0.5688902388963853, 0.596448730900158, 1.5309152952910559, -0.2059339861429303, 0.569825988876055
};

double beta0 = 1500.0;
double mu0 = 1.0/3.0;

void mlsp2_full(double *GPU_hamiltonian,
                double *GPU_densityMatrix,
                int N,
                double beta,
                double mu,
                precision_t precision,
                refine_t refinement)
{    
    double emin, emax;
    
    gershgorin_v2(N, GPU_hamiltonian, &emin, &emax);
    std::cout << "emin: " << emin << " emax: " << emax << std::endl;

    std::cout << "H(0) " << GPU_hamiltonian[0] << " N " << N << " beta " << beta << "mu "<< mu << std::endl;
    

    // Quick & Dirty rescaling from mlsp2.py.
    // Can be reduced to 1 shift_and_scale call by reduction:
    // s&s(m,b,x) = mx+b
    // s&s(m,b,s&s(M,B,x)) = m(Mx+B)+b = mMx + mB+b = s&s(mM,mB+b,x)

    // leaving as is now for debugging purposes.

    // primed variables given by Eqs. 27 - 29

    double shift = emax/(emax-emin);
    double scale = -1/(emax-emin);

    int numthds = 512;
    int numblks = int(ceil(double(N * N) / double(numthds)));

    // shift_and_scale<<<numblks, numthds>>>(N, GPU_hamiltonian, shift, scale);

    double mu_prime   = (emax - mu)/(emax - emin);
    double beta_prime = (emax - emin) * beta;
    
    // condition for validity by Eq. 43
    // assert beta_prime <= (2/3) * beta0
    // assert emin < mu < emax

    // flip given by Eq. 44 if mu' > 0.5
    bool mu_switch = (mu_prime > 0.5);
    if(mu_switch){
        // -x + 1
        shift = -shift + 1.0;
        scale = -scale;
        // shift_and_scale<<<numblks, numthds>>>(N, GPU_hamiltonian, shift, scale);
        mu_prime = 1.0 - mu_prime;
    }

    // H0 given by Eq. 30
    scale *= beta_prime/beta0;
    shift *= beta_prime/beta0;
    shift += -mu_prime*beta_prime/beta0 + mu0;

    shift_and_scale<<<numblks, numthds>>>(N, GPU_hamiltonian, shift, scale);
    
    mlsp2(model26, GPU_hamiltonian, GPU_densityMatrix, 26, N, precision, refinement);

    if (!mu_switch){
        shift =  1.0;
        scale = -1.0;
        shift_and_scale<<<numblks, numthds>>>(N, GPU_densityMatrix, shift, scale);
    }
}

void mlsp2_newton(double *GPU_hamiltonian,
                  double *GPU_densityMatrix,
                  int N,
                  double beta,
                  double mu_guess,
                  double f_occ,
                  double epsilon,
                  int max_iter,
                  precision_t precision,
                  refine_t refinement)
{
    double error = 10. * epsilon;
    int i = 0;

    double Tr, Tr2, d_mu;

    std::cout << 'BEGIN NEWTON LOOP' << std::endl;
    while (abs(error) > epsilon && i++ < max_iter){

        // generate initial guess at D
        mlsp2_full(GPU_hamiltonian,
                   GPU_densityMatrix,
                   N, beta, mu_guess,
                   precision, refinement);

        // Tr(D), Tr(D^2)
        openacc_trace (&Tr,  GPU_densityMatrix, N);
        openacc_trace2(&Tr2, GPU_densityMatrix, N);

        // occupation error
        error = double(N) * f_occ - Tr;

        // Guess given by Newton using f'(H) = -beta * D(I-D)
        d_mu = error/(Tr - Tr2)/beta;
        d_mu = min(d_mu,  0.5 * abs(mu_guess));
        d_mu = max(d_mu, -0.5 * abs(mu_guess));
        
        std::cout << mu_guess << ' ' << d_mu << " iter: " << i << std::endl;

        mu_guess += d_mu;

    }
}

void mlsp2(double *model,
           double *GPU_hamiltonian,
           double *GPU_densityMatrix,
           int numLayers,
           int N,
           precision_t precision,
           refine_t refinement)
{
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    float one_f = 1.0, zero_f = 0.0, neg_one_f = -1.0;

    // Cublas Handle
    nvtxRangePushA("build cublas handle");
    cublasHandle_t handle;
    CUBLAS_CHECK_ERR(cublasCreate(&handle));
    nvtxRangePop();

    nvtxRangePushA("declare memory");
    // Declare Memory
    float *GPU_Si, *GPU_Si_squared, *GPU_identityMatrix, *sbuf1, *sbuf2, *GPU_accumulationMatrix;
    half *hbuf1, *hbuf2;

    nvtxRangePop();

    // Allocate some host memory
    nvtxRangePushA("initialize memory");

    // Allocate device memory
    CUDA_CHECK_ERR(cudaMalloc(&GPU_Si, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_Si_squared, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_identityMatrix, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_accumulationMatrix, N * N * sizeof(float)));
    // Initialize the accumulation matrix to zero
    CUDA_CHECK_ERR(cudaMemset(GPU_accumulationMatrix, 0, N * N * sizeof(float)));

    // Allocate Buffers
    CUDA_CHECK_ERR(cudaMalloc(&sbuf1, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&sbuf2, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&hbuf1, N * N * sizeof(half)));
    CUDA_CHECK_ERR(cudaMalloc(&hbuf2, N * N * sizeof(half)));
    nvtxRangePop();

    // Define blk,thd grid size
    int numthds = 512;
    int numblks = int(ceil(double(N * N) / double(numthds)));

    // Initialize Hamiltonian and identity

    // build Identity on dev
    setToIdentityMatrix<<<numblks, numthds>>>(GPU_identityMatrix, N);

    // cast GPU_hamiltonian from double to float
    doubleToFloat<<<numblks, numthds>>>(GPU_hamiltonian, GPU_Si, N);
    CUDA_CHECK_ERR(cudaMemcpy(sbuf1, GPU_Si, N * N * sizeof(float), cudaMemcpyDeviceToDevice));

    printf("%d\n", N );

    float a, b, c, d;

    nvtxRangePop();

    nvtxRangePushA("Main loop");
    for (int iter = 0; iter < numLayers; ++iter)
    {

        a = model[iter];
        b = model[numLayers + iter];
        c = model[2 * numLayers + iter];
        d = model[3 * numLayers + iter];

        nvtxRangePushA("TC matmul");
        if (precision == fp32)
        {
            CUBLAS_CHECK_ERR(cublasSgemm(handle,
                                         CUBLAS_OP_N, CUBLAS_OP_N,
                                         N, N, N,
                                         &one_f,
                                         GPU_Si, N,
                                         GPU_Si, N,
                                         &zero_f,
                                         GPU_Si_squared, N));


            // cublasSideMode_t side = CUBLAS_SIDE_LEFT;
            // cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

            // CUBLAS_CHECK_ERR(cublasSsymm(
            //     handle, 
            //     side, uplo,
            //     N, N,
            //     &one_f,
            //     GPU_Si, N,
            //     GPU_Si, N,
            //     &zero_f,
            //     GPU_Si_squared, N));
                
        }
        else if (precision == fp16_fp32)
        {
            // tcoreSPGemmSymm(handle,
            //                 N,
            //                 GPU_Si,
            //                 hbuf1, hbuf2,
            //                 sbuf1, sbuf2,
            //                 GPU_Si_squared);
            

            // if ( abs(d) > 1e-10) 
            // {
                tcoreSPGemmSymmAlt(handle,
                                   N,
                                   GPU_Si,
                                   hbuf1, hbuf2,
                                   sbuf1,
                                   GPU_Si_squared, one_f);
            // }
            // else 
            // {
            //     tcoreHPGemmSymm(handle,
            //                        N,
            //                        GPU_Si,
            //                        hbuf1, hbuf2,
            //                        GPU_Si_squared, one_f);
            // }

        };
        nvtxRangePop();

        // Accumulate d * GPU_S₀ into the accumulation matrix
        CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                     CUBLAS_OP_N, CUBLAS_OP_N,
                                     N, N,
                                     &d, // Scale factor for GPU_S₀
                                     GPU_Si, N,
                                     &one_f, // Accumulate into GPU_accumulationMatrix
                                     GPU_accumulationMatrix, N,
                                     GPU_accumulationMatrix, N));

        //std::cout << a << std::endl;

        CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                     CUBLAS_OP_N, CUBLAS_OP_N,
                                     N, N,
                                     &a,
                                     GPU_Si_squared, N,
                                     &b,
                                     GPU_Si, N,
                                     GPU_Si, N));
           
        // shift_v1<<<numblks, numthds>>>(GPU_Si, N, c);

        int num_blks = int(ceil(double(N) / double(numthds)));
        shift_v2<<<num_blks, numthds>>>(GPU_Si, N, c);

        // CUBLAS_CHECK_ERR(cublasSgeam(handle,
        //                              CUBLAS_OP_N, CUBLAS_OP_N,
        //                              N, N,
        //                              &one_f, // Keep GPU_S₀ as-is (scaled by 1)
        //                              GPU_Si, N,
        //                              &c,
        //                              GPU_identityMatrix, N,
        //                              GPU_Si, N));

    }
    nvtxRangePop();

    // Accumulate d * GPU_S₀ into the accumulation matrix
    CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                CUBLAS_OP_N, CUBLAS_OP_N,
                                N, N,
                                &one_f, // Scale factor for GPU_S₀
                                GPU_Si, N,
                                &one_f, // Accumulate into GPU_accumulationMatrix
                                GPU_accumulationMatrix, N,
                                GPU_accumulationMatrix, N));

    // Subtract GPU_accumulationMatrix from GPU_identityMatrix and store in GPU_densityMatrix
    CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                 CUBLAS_OP_N, CUBLAS_OP_N,
                                 N, N,
                                 &neg_one_f, // Scale factor for GPU_accumulationMatrix (-1)
                                 GPU_accumulationMatrix, N,
                                 &one_f, // Scale factor for GPU_identityMatrix (1)
                                 GPU_identityMatrix, N,
                                 GPU_accumulationMatrix, N)); // Store result in GPU_densityMatrix

    floatToDouble<<<numblks, numthds>>>(GPU_accumulationMatrix, GPU_densityMatrix, N);

    // Save GPU_densityMatrix to disk as plain text before cleaning up resources
    // double *host_densityMatrix = (double *)malloc(N * N * sizeof(double));
    // CUDA_CHECK_ERR(cudaMemcpy(host_densityMatrix, GPU_densityMatrix, N * N * sizeof(double), cudaMemcpyDeviceToHost));

    // FILE *file = fopen("density_matrix.txt", "w");
    // if (file != NULL)
    // {
    //     for (int i = 0; i < N; i++)
    //     {
    //         for (int j = 0; j < N; j++)
    //         {
    //             fprintf(file, "%lf ", host_densityMatrix[i * N + j]);
    //         }
    //         fprintf(file, "\n"); // New line at the end of each row
    //     }
    //     fclose(file);
    // }
    // else
    // {
    //     std::cerr << "Error: Could not open file for writing" << std::endl;
    // }
    // free(host_densityMatrix);

    //nvtxRangePushA("Handle destroy");
    // Destroy CUBLAS handle
    //CUBLAS_CHECK_ERR(cublasDestroy(handle));
    //nvtxRangePop();

    //nvtxRangePushA("cudaFree");
    // Free device memory
/*    CUDA_CHECK_ERR(cudaFree(GPU_Si));
    CUDA_CHECK_ERR(cudaFree(GPU_Si_squared));
    CUDA_CHECK_ERR(cudaFree(GPU_identityMatrix));
    CUDA_CHECK_ERR(cudaFree(GPU_accumulationMatrix));
    CUDA_CHECK_ERR(cudaFree(sbuf1));
    CUDA_CHECK_ERR(cudaFree(sbuf2));
    CUDA_CHECK_ERR(cudaFree(hbuf1));
    CUDA_CHECK_ERR(cudaFree(hbuf2));*/
    //nvtxRangePop();

    // Record the stop event
    cudaEventRecord(stop);

    // Synchronize and measure elapsed time
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "MLSP2 time from CUDA = " << milliseconds / 1000.0 << " seconds" << std::endl;

    // Destroy CUDA events
    CUDA_CHECK_ERR(cudaEventDestroy(start));
    CUDA_CHECK_ERR(cudaEventDestroy(stop));
}

// essentially DNN-SP2 form
void mlsp2_alt(double *model,
           double *GPU_hamiltonian,
           double *GPU_densityMatrix,
           int numLayers,
           int N,
           precision_t precision,
           refine_t refinement)
{

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    float one_f = 1.0, zero_f = 0.0, neg_one_f = -1.0;

    // Cublas Handle
    nvtxRangePushA("build cublas handle");
    cublasHandle_t handle;
    CUBLAS_CHECK_ERR(cublasCreate(&handle));
    nvtxRangePop();

    nvtxRangePushA("declare memory");
    // Declare Memory
    float *GPU_Si, *GPU_Si_squared, *GPU_identityMatrix, *sbuf1, *sbuf2, *GPU_accumulationMatrix;
    half *hbuf1, *hbuf2;

    nvtxRangePop();

    // Allocate some host memory
    nvtxRangePushA("initialize memory");

    // Allocate device memory
    CUDA_CHECK_ERR(cudaMalloc(&GPU_Si, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_identityMatrix, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_accumulationMatrix, N * N * sizeof(float)));
    // Initialize the accumulation matrix to zero
    CUDA_CHECK_ERR(cudaMemset(GPU_accumulationMatrix, 0, N * N * sizeof(float)));

    // Allocate Buffers
    CUDA_CHECK_ERR(cudaMalloc(&hbuf1, N * N * sizeof(half)));
    CUDA_CHECK_ERR(cudaMalloc(&hbuf2, N * N * sizeof(half)));
    CUDA_CHECK_ERR(cudaMalloc(&sbuf1, N * N * sizeof(float)));
    nvtxRangePop();

    // Define blk,thd grid size
    int numthds = 512;
    int numblks = int(ceil(double(N * N) / double(numthds)));

    // Initialize Hamiltonian and identity

    // build Identity on dev
    setToIdentityMatrix<<<numblks, numthds>>>(GPU_identityMatrix, N);

    // cast GPU_hamiltonian from double to float
    doubleToFloat<<<numblks, numthds>>>(GPU_hamiltonian, GPU_Si, N);

    printf("%d\n", N );

    float a, b, c, d, acc_err;

    a = 0.0; b = 0.0, c=0.0;
    acc_err = 1.0;

    nvtxRangePop();

    nvtxRangePushA("Main loop");
    for (int iter = 0; iter < numLayers; ++iter)
    {

        // carry over from previous layer
        b = c - a*b*b;

        a  = model[                iter];
        b += model[    numLayers + iter] / (2.*a);
        c  = model[2 * numLayers + iter];
        d  = model[3 * numLayers + iter];

        // so turns out this is just completely unstable. Probably going to have to train the DNN-SP2 form independently.
        std::cout << a << b << c << d << acc_err << std::endl; 

        // Accumulate d * GPU_S₀ into the accumulation matrix
        CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                     CUBLAS_OP_N, CUBLAS_OP_N,
                                     N, N,
                                     &d, // Scale factor for GPU_S₀
                                     GPU_Si, N,
                                     &one_f, // Accumulate into GPU_accumulationMatrix
                                     GPU_accumulationMatrix, N,
                                     GPU_accumulationMatrix, N));

        // Si <- Si + b/2a I
        CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                     CUBLAS_OP_N, CUBLAS_OP_N,
                                     N, N,
                                     &one_f,
                                     GPU_Si, N,
                                     &b,
                                     GPU_identityMatrix, N,
                                     GPU_Si, N));


        // Si <- a ( Si^2 + 2 (b/2a) Si + (b/2a)^2 I )
        // i.e. a Si^2 + b Si + b^2/4a I  

        // Si is missing a factor of (c - b^2/4a)
        // thus, for the next layer, b_new <- b_new + c - b^2/4a 

        // accumulator is missing d * (what Si is missing)
        acc_err += d * ( c - a * b * b);

        nvtxRangePushA("TC matmul");
        if (precision == fp32)
        {
            CUBLAS_CHECK_ERR(cublasSgemm(handle,
                                         CUBLAS_OP_N, CUBLAS_OP_N,
                                         N, N, N,
                                         &a,
                                         GPU_Si, N,
                                         GPU_Si, N,
                                         &zero_f,
                                         GPU_Si, N));
        }
        else if (precision == fp16_fp32)
        {
            tcoreSPGemmSymmAlt(handle,
                            N,
                            GPU_Si,
                            hbuf1, hbuf2,
                            sbuf1,
                            GPU_Si, a);
        };
        nvtxRangePop();

    }
    nvtxRangePop();

    // Accumulate d * GPU_S₀ into the accumulation matrix
    CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                CUBLAS_OP_N, CUBLAS_OP_N,
                                N, N,
                                &one_f, // Scale factor for GPU_S₀
                                GPU_Si, N,
                                &one_f, // Accumulate into GPU_accumulationMatrix
                                GPU_accumulationMatrix, N,
                                GPU_accumulationMatrix, N));

    // Subtract GPU_accumulationMatrix from GPU_identityMatrix and store in GPU_densityMatrix
    CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                 CUBLAS_OP_N, CUBLAS_OP_N,
                                 N, N,
                                 &neg_one_f, // Scale factor for GPU_accumulationMatrix (-1)
                                 GPU_accumulationMatrix, N,
                                 &acc_err, // Scale factor for GPU_identityMatrix
                                 GPU_identityMatrix, N,
                                 GPU_accumulationMatrix, N)); // Store result in GPU_densityMatrix

    floatToDouble<<<numblks, numthds>>>(GPU_accumulationMatrix, GPU_densityMatrix, N);

    // Save GPU_densityMatrix to disk as plain text before cleaning up resources
    // double *host_densityMatrix = (double *)malloc(N * N * sizeof(double));
    // CUDA_CHECK_ERR(cudaMemcpy(host_densityMatrix, GPU_densityMatrix, N * N * sizeof(double), cudaMemcpyDeviceToHost));

    // FILE *file = fopen("density_matrix.txt", "w");
    // if (file != NULL)
    // {
    //     for (int i = 0; i < N; i++)
    //     {
    //         for (int j = 0; j < N; j++)
    //         {
    //             fprintf(file, "%lf ", host_densityMatrix[i * N + j]);
    //         }
    //         fprintf(file, "\n"); // New line at the end of each row
    //     }
    //     fclose(file);
    // }
    // else
    // {
    //     std::cerr << "Error: Could not open file for writing" << std::endl;
    // }
    // free(host_densityMatrix);

    //nvtxRangePushA("Handle destroy");
    // Destroy CUBLAS handle
    //CUBLAS_CHECK_ERR(cublasDestroy(handle));
    //nvtxRangePop();

    //nvtxRangePushA("cudaFree");
    // Free device memory
/*    CUDA_CHECK_ERR(cudaFree(GPU_Si));
    CUDA_CHECK_ERR(cudaFree(GPU_Si_squared));
    CUDA_CHECK_ERR(cudaFree(GPU_identityMatrix));
    CUDA_CHECK_ERR(cudaFree(GPU_accumulationMatrix));
    CUDA_CHECK_ERR(cudaFree(sbuf1));
    CUDA_CHECK_ERR(cudaFree(sbuf2));
    CUDA_CHECK_ERR(cudaFree(hbuf1));
    CUDA_CHECK_ERR(cudaFree(hbuf2));*/
    //nvtxRangePop();

    // Record the stop event
    cudaEventRecord(stop);

    // Synchronize and measure elapsed time
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "MLSP2 time from CUDA = " << milliseconds / 1000.0 << " seconds" << std::endl;

    // Destroy CUDA events
    CUDA_CHECK_ERR(cudaEventDestroy(start));
    CUDA_CHECK_ERR(cudaEventDestroy(stop));
}

void mlsp2_alt2(double *model,
           double *GPU_hamiltonian,
           double *GPU_densityMatrix,
           int numLayers,
           int N,
           precision_t precision,
           refine_t refinement,
           bool flip = true)
{

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    float one_f = 1.0, zero_f = 0.0, neg_one_f = -1.0;

    // Cublas Handle
    nvtxRangePushA("build cublas handle");
    cublasHandle_t handle;
    CUBLAS_CHECK_ERR(cublasCreate(&handle));
    nvtxRangePop();

    nvtxRangePushA("declare memory");
    // Declare Memory
    float *GPU_Si, *GPU_Si_squared, *GPU_identityMatrix, *sbuf1, *sbuf2, *GPU_accumulationMatrix;
    half *hbuf1, *hbuf2;

    nvtxRangePop();

    // Allocate some host memory
    nvtxRangePushA("initialize memory");

    // Allocate device memory
    CUDA_CHECK_ERR(cudaMalloc(&GPU_Si, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_Si_squared, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_identityMatrix, N * N * sizeof(float)));
    CUDA_CHECK_ERR(cudaMalloc(&GPU_accumulationMatrix, N * N * sizeof(float)));
    // Initialize the accumulation matrix to zero
    CUDA_CHECK_ERR(cudaMemset(GPU_accumulationMatrix, 0, N * N * sizeof(float)));

    // Allocate Buffers
    CUDA_CHECK_ERR(cudaMalloc(&hbuf1, N * N * sizeof(half)));
    CUDA_CHECK_ERR(cudaMalloc(&hbuf2, N * N * sizeof(half)));
    nvtxRangePop();

    // Define blk,thd grid size
    int numthds = 512;
    int numblks = int(ceil(double(N * N) / double(numthds)));

    // Initialize Hamiltonian and identity

    // build Identity on dev
    setToIdentityMatrix<<<numblks, numthds>>>(GPU_identityMatrix, N);

    // cast GPU_hamiltonian from double to float
    doubleToFloat<<<numblks, numthds>>>(GPU_hamiltonian, GPU_Si, N);
    
    printf("%d\n", N );

    float a, b, c, d;

    nvtxRangePop();

    nvtxRangePushA("Main loop");
    for (int iter = 0; iter < numLayers; ++iter)
    {

        a = model[                iter];
        b = model[    numLayers + iter];
        c = model[2 * numLayers + iter];
        d = model[3 * numLayers + iter];

        // Accumulate d * GPU_S₀ into the accumulation matrix
        if (abs(d) > 1e-8) {
        CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                     CUBLAS_OP_N, CUBLAS_OP_N,
                                     N, N,
                                     &d, // Scale factor for GPU_S₀
                                     GPU_Si, N,
                                     &one_f, // Accumulate into GPU_accumulationMatrix
                                     GPU_accumulationMatrix, N,
                                     GPU_accumulationMatrix, N));
        }

        nvtxRangePushA("TC matmul");
        if (precision == fp32)
        {
            // Si^2 <- Si * Si
            CUBLAS_CHECK_ERR(cublasSgemm(handle,
                                         CUBLAS_OP_N, CUBLAS_OP_N,
                                         N, N, N,
                                         &one_f,
                                         GPU_Si, N,
                                         GPU_Si, N,
                                         &zero_f,
                                         GPU_Si_squared, N));

            // Si <- a * Si^2 + b * Si
            CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                        CUBLAS_OP_N, CUBLAS_OP_N,
                                        N, N,
                                        &a, 
                                        GPU_Si_squared, N,
                                        &b,
                                        GPU_Si, N,
                                        GPU_Si, N));

        }
        else if (precision == fp16_fp32)
        {
            // Si^2 <- a * Si * Si + b * Si
            tcoreSPGemmSymmAlt2(handle,
                            N,
                            GPU_Si,
                            hbuf1, hbuf2,
                            GPU_Si_squared,
                            a, b);
            
            // swap overwritten Si with Si^2
            float* temp = GPU_Si_squared;
            GPU_Si_squared = GPU_Si;
            GPU_Si = temp;
        };
        nvtxRangePop();
        
        // Si <- Si + c * I
        if (abs(c) > 1e-10) {
            CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                        CUBLAS_OP_N, CUBLAS_OP_N,
                                        N, N,
                                        &one_f,
                                        GPU_Si, N,
                                        &c,
                                        GPU_identityMatrix, N,
                                        GPU_Si, N));
            
            std::cout << "nonzero c" << std::endl;
        }

    }
    nvtxRangePop();

    // Accumulate d * GPU_S₀ into the accumulation matrix
    CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                CUBLAS_OP_N, CUBLAS_OP_N,
                                N, N,
                                &one_f, // Scale factor for GPU_S₀
                                GPU_Si, N,
                                &one_f, // Accumulate into GPU_accumulationMatrix
                                GPU_accumulationMatrix, N,
                                GPU_accumulationMatrix, N));

    // Subtract GPU_accumulationMatrix from GPU_identityMatrix and store in GPU_densityMatrix
    if (flip) {
    CUBLAS_CHECK_ERR(cublasSgeam(handle,
                                 CUBLAS_OP_N, CUBLAS_OP_N,
                                 N, N,
                                 &neg_one_f, // Scale factor for GPU_accumulationMatrix (-1)
                                 GPU_accumulationMatrix, N,
                                 &one_f, // Scale factor for GPU_identityMatrix (1)
                                 GPU_identityMatrix, N,
                                 GPU_accumulationMatrix, N)); // Store result in GPU_densityMatrix
    }

    floatToDouble<<<numblks, numthds>>>(GPU_accumulationMatrix, GPU_densityMatrix, N);

    // Save GPU_densityMatrix to disk as plain text before cleaning up resources
    // double *host_densityMatrix = (double *)malloc(N * N * sizeof(double));
    // CUDA_CHECK_ERR(cudaMemcpy(host_densityMatrix, GPU_densityMatrix, N * N * sizeof(double), cudaMemcpyDeviceToHost));

    // FILE *file = fopen("density_matrix.txt", "w");
    // if (file != NULL)
    // {
    //     for (int i = 0; i < N; i++)
    //     {
    //         for (int j = 0; j < N; j++)
    //         {
    //             fprintf(file, "%lf ", host_densityMatrix[i * N + j]);
    //         }
    //         fprintf(file, "\n"); // New line at the end of each row
    //     }
    //     fclose(file);
    // }
    // else
    // {
    //     std::cerr << "Error: Could not open file for writing" << std::endl;
    // }
    // free(host_densityMatrix);

    nvtxRangePushA("Handle destroy");
    // Destroy CUBLAS handle
    CUBLAS_CHECK_ERR(cublasDestroy(handle));
    nvtxRangePop();

    nvtxRangePushA("cudaFree");
    // Free device memory
    CUDA_CHECK_ERR(cudaFree(GPU_Si));
    CUDA_CHECK_ERR(cudaFree(GPU_Si_squared));
    CUDA_CHECK_ERR(cudaFree(GPU_identityMatrix));
    CUDA_CHECK_ERR(cudaFree(GPU_accumulationMatrix));
    CUDA_CHECK_ERR(cudaFree(sbuf1));
    CUDA_CHECK_ERR(cudaFree(sbuf2));
    CUDA_CHECK_ERR(cudaFree(hbuf1));
    CUDA_CHECK_ERR(cudaFree(hbuf2));
    nvtxRangePop();

    // Record the stop event
    cudaEventRecord(stop);

    // Synchronize and measure elapsed time
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "MLSP2 time from CUDA = " << milliseconds / 1000.0 << " seconds" << std::endl;

    // Destroy CUDA events
    CUDA_CHECK_ERR(cudaEventDestroy(start));
    CUDA_CHECK_ERR(cudaEventDestroy(stop));
}

