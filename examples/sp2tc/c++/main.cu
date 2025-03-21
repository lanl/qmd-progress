#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <typeinfo>
#include <cuda.h>
#include <cublas_v2.h>
#include <cuda_fp16.h>
#include <cmath>
#include <vector>
#include "tcore_hp_emulator.cuh"
#include "linalg_tools.cuh"
#include "prg_sp2_tensorcore.cuh"


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

    std::vector<float> Idemp_Error;
    
    std::cout << "Matrix Size: " << N << std::endl;
    std::cout << "Occupied orbitals: " << Nocc << std::endl;

    // Set GPU
    int device = 0;
    cudaSetDevice(device);

    // Allocate memory
    float *H;
    H = (float*) malloc(N*N*sizeof(float));
    double *D; 
    D = (double*) malloc(N*N*sizeof(double));
    
    // Produce Hamiltonian 
    produce_hamiltonian(N,H);

    float idemtol=1e-6;
    char sp2conv,prec;
    int verbose=0;
    float bndfil=float(Nocc)/float(N);
    
    prg_sp2_tensorcore(N,H,D,eps,bndfil,1,1000,sp2conv,idemtol,verbose);
 
    std::cout << D[0] << std::endl; 
}



