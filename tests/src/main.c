
#include "bml.h"
#include "prg_progress_mod.h"
#include <stdio.h>

void main(){
    int norb, mdim, verbose, i, j;
    bml_matrix_t *ham;
    bml_matrix_t *rho, *rho1_bml;
    bml_matrix_t *rho_ortho_bml;
    bml_matrix_t *zmat_bml, *zk1_bml, *zk2_bml;
    bml_matrix_t *zk3_bml, *zk4_bml, *zk5_bml, *zk6_bml;
    bml_matrix_t *nonortho_ham_bml;
    bml_matrix_t *over_bml, *orthox_bml, *g_bml;
    bml_matrix_t *aux_bml, *pcm_bml, *pcm_ref_bml, *inv_bml[10];
    
    // graph_partitioning_t gp;
    // system_type mol;

    bml_matrix_type_t matrix_type;
    bml_matrix_precision_t precision;
    bml_distribution_mode_t distrib_mode = sequential;

    char bml_type[21]; // 20 char + null termination '\0'
    char sp2conv_sue[21]; // 20 char + null termination '\0'
    char test[51]; // 50 char + null termination '\0'
    char dummy[10][21]; // 20 char + null termination '\0' each

    double threshold, gthreshold, idempotency;
    double sp2tol, idempotency_tol;
    double bndfil, mu, tscale, tracelimit, beta;
    /*
    double* ham = NULL;    // ham(:,:)
    double* zmat = NULL;   // zmat(:,:)
    double* nonortho_ham = NULL;   // nonortho_ham(:,:)
    double* over = NULL;   // over(:,:)
    double* rho_ortho = NULL;   // rho_ortho(:,:)
    double* trace = NULL;   // trace(:)
    double* rho = NULL;   // rho(:,:)
    double* eigenvalues = NULL;   // eigenvalues(:)
    double* row = NULL;   // row(:)
    double error_calc, error_tol, errlimit;
    double eps, beta0, nocc, kbt;
    double mineval, maxeval, occerrlimit;
    double drho, drho_ref;
    double* gbnd = NULL;   // gbnd(:)

    int minsp2iter, icount, nodesPerPart, occsteps;
    int norecs, nsiter, occiter, numparts;
    int maxsp2iter, npts, sp2all_timer, sp2all_timer_init;
    int* pp = NULL;   // pp(:)
    int* signlist = NULL;   // signlist(:)

    double* vv = NULL;   // vv(:)

    char sp2conv[11]; // 10 char + null termination '\0'

    //tbparams_type tbparams;

    char (*intKind)[4] = NULL; // 3 char + null termination '\0'
    char (*TypeA)[3] = NULL; // 2 char + null termination '\0'
    char (*TypeB)[3] = NULL; // 2 char + null termination '\0'

    double** onsitesH = NULL;   // onsitesH(:,:)
				// 
    */
    printf("Hello World\n");
    matrix_type = dense;
    precision = double_real;
    norb = 600;
    threshold = 1.0e-9;
    bndfil = 0.666666666666666666;
    idempotency_tol = 1.0e-8;

    rho = bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
    ham = bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
    bml_read_bml_matrix(ham, "hamiltonian.mtx");

    bml_matrix_t *ham_t = NULL;
    ham_t = bml_transpose_new(ham);
    bml_add(ham, ham_t, 0.5, 0.5, 0.0);

    LOG_INFO("(A + A_t)/2 = \n");
    int max_col=10;
    int max_row = 10;
    //bml_print_bml_matrix(ham, 0, max_row, 0, max_col);

    prg_build_density_T0(ham, rho, threshold, bndfil);
    LOG_INFO("rho = \n");
    bml_print_bml_matrix(rho, 0, max_row, 0, max_col);
    
    //bml_scale(0.5, rho);
    //prg_check_idempotency(rho, threshold,idempotency);

}
