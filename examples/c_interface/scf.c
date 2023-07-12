#include <stdlib.h>
#include "bml.h"

void
main(
    )
{
    // SCF example using c interface of progress lib
    // TBA
    int dp = sizeof(double);    // equivalent of kind(1.0d0)
    int i, nel, norb, j;
    int **hindex = NULL;        // allocatable in Fortran, use malloc or calloc in C when actual size is known
    char (
    *intKind)[4] = NULL;        // 3 char + null termination '\0'

    double bndfil, scferror;
    double *charges_old = NULL;
    double **coul_forces_k = NULL;
    double **coul_forces_r = NULL;
    double *coul_pot_k = NULL;
    double *coul_pot_r = NULL;
    double **dqin = NULL;
    double **dqout = NULL;

    bml_matrix_t ham0_bml, ham_bml, orthoh_bml, orthop_bml;
    bml_matrix_t over_bml, rho_bml, zmat_bml;

    latte_type lt;
    sp2data_type sp2;
    system_type sy;
    tbparams_type tb;

    intpairs_type **intPairsH = NULL;
    intpairs_type **intPairsS = NULL;

    double **onsitesH = NULL;
    double **onsitesS = NULL;

    char (
    *TypeA)[3] = NULL;          // 2 char + null termination '\0'
    char (
    *TypeB)[3] = NULL;          // 2 char + null termination '\0'


}
