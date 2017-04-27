#if ! _OPENMP
#error Missing OpenMP
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int size;
    int rank;
    int counter;
    int provided;

    //MPI_Init(&argc, &argv);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("MPI size %d\n", size);
        printf("MPI_THREAD_MULTIPLE = %d\n", MPI_THREAD_MULTIPLE);
        printf("provided = %d\n", provided);
    }

#pragma omp parallel
    {
#pragma omp master
        printf("Running on %d threads\n", omp_get_num_threads());
    }

    printf("MPI rank %d\n", rank);
#pragma omp parallel private(counter)
    {
        printf("Thread %d starting up\n", omp_get_thread_num());
        while (1) {
            counter++;
            counter = counter%100;
        }
    }

    MPI_Finalize();
}
