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
#include "../include/tcore_hp_emulator.cuh"
#include "../include/linalg_tools.cuh"
#include "../include/prg_sp2_tensorcore.cuh"



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

    float eps = 1e-16;
 
    std::cout << "Mat Size: " << N << std::endl;
    std::cout << "Occupied orbitals: " << Nocc << std::endl;

    
    // set cublas math mode (i.e. turn on/off tensorcores)
    //cublasStatus_t cublasStat = cublasSetMathMode(handle, CUBLAS_TENSOR_OP_MATH);

    // Declre and allocate memory for the Hamiltonian and Density matrix
    float *H; double *D;
    H = (float*) malloc(N*N*sizeof(float));
    D = (double*) malloc(N*N*sizeof(double));
    
    // Produce Hamiltonian
    std::cout << "Loading Hamiltonian..." << std::endl;
    produce_hamiltonian(N,H);
    
    // Get device id
    //cudaGetDevice(&device); 

    float idemtol=1e-16;
    char sp2conv;
    int verbose=0;
    float bndfil=float(Nocc);//float(Nocc)/float(N);
    std::cout << N << Nocc << bndfil << std::endl;
    prg_sp2_tensorcore(N,H,D,eps, bndfil,1,1000,sp2conv,idemtol,verbose);
}



