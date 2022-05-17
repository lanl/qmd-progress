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
#include <ctime>
#include <cmath>
#include <vector>
#include <chrono>
#include "tcore_hp_emulator.cuh"
#include "linalg_tools.cuh"
#include "prg_sp2_tensorcore.cuh"


// Test Hamiltonian
void produce_hamiltonian (const unsigned N, float *X) {
    for(int i=0; i<N; ++i) {
        for(int j=i; j<N; ++j) {
            X[i+j*N] = exp(-0.5f*abs((float)(i-j)))*sin((float)(i+1));
            X[j+i*N] = X[i+j*N];
        }
    }
};


int main(int argc, char *argv[])
{

    // Matrix size
    size_t N = atoi(argv[1]);
    size_t Nocc = atoi(argv[2]);

    int Stopp = 0;
    int Kvot = 0;
    int iter = 0;
    int Pur_Start = 0;
    float eps = 1e-16;


    std::vector<float> Idemp_Error;
    
    std::cout << "Mat Size: " << N << std::endl;
    std::cout << "Occupied orbitals: " << Nocc << std::endl;

    // Set GPU
    int device = 0;
    cudaSetDevice(device);

    // Cublas Handle
    cublasHandle_t handle;
    cublasCreate(&handle);
    

    cublasStatus_t cublasStat = cublasSetMathMode(handle, CUBLAS_TENSOR_OP_MATH);

    float *H;
    H = (float*) malloc(N*N*sizeof(float));
    double *D; 
    D = (double*) malloc(N*N*sizeof(double));
    

    // Produce Hamiltonian and Identity matrix 
    std::cout << "Loading Hamiltonian..." << std::endl;
    produce_hamiltonian(N,H);
    
    // Get device id
    cudaGetDevice(&device); 

    float idemtol=1e-16;
    char sp2conv;
    int verbose=0;
    float bndfil=float(Nocc)/float(N);
    std::cout << bndfil << std::endl; 
    prg_sp2_tensorcore(N,H,D,eps,bndfil,1,1000,sp2conv,idemtol,verbose);
    
}



