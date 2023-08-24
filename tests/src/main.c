
#include "string.h"
#include "bml.h"
#include "prg_progress_mod.h"
#include <stdio.h>

int
main(
    int argc,
    char *argv[])
{
    if (argc < 2)
    {
        printf("Please provide an argument.\n");
        return 1;
    }
    char *test = argv[1];

    int norb, mdim, verbose, i, j;
    bml_matrix_t *ham;
    bml_matrix_t *rho, *rho1;
    bml_matrix_t *rho_ortho_bml;
    bml_matrix_t *zmat, *zk1_bml, *zk2_bml;
    bml_matrix_t *zk3_bml, *zk4_bml, *zk5_bml, *zk6_bml;
    bml_matrix_t *nonortho_ham_bml;
    bml_matrix_t *over_bml, *g_bml;
    bml_matrix_t *aux_bml, *pcm, *pcm_ref, *inv_bml[10];
    bml_matrix_t *orthox_bml;

    bml_matrix_type_t matrix_type;
    bml_matrix_precision_t precision;
    bml_distribution_mode_t distrib_mode = sequential;

    char bml_type[21];          // 20 char + null termination '\0'
    char sp2conv_sue[21];       // 20 char + null termination '\0'
    char dummy[10][21];         // 20 char + null termination '\0' each
    char sp2conv[10] = "Rel";

    double threshold, gthreshold;
    double idempotency;
    double sp2tol, idempotency_tol;
    double bndfil, tscale, tracelimit, beta;
    double error_calc, error_tol, errlimit;
    double maxeval, mineval, occerrlimit;
    double beta0;

    matrix_type = dense;
    precision = double_real;
    norb = 600;
    threshold = 1.0e-9;
    bndfil = 0.666666666666666666;
    idempotency_tol = 1.0e-8;

    double kbt = 0.01;
    double mu = 0.0;
    double scale_factor = 0.5;
    double *eigenvalues = malloc(norb * sizeof(double));
    int max_col = 10;
    int max_row = 10;
    int minsp2iter = 25;
    int maxsp2iter = 100;

    verbose = 1;
    sp2tol = 1.e-10;

    prg_progress_init();

    if (strcmp(test, "prg_density_c") == 0)
    {
        LOG_INFO
            ("Testing the construction of the density matrix from density_mod \n");
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);

        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        bml_matrix_t *ham_t = NULL;
        ham_t = bml_transpose_new(ham);
        bml_add(ham, ham_t, 0.5, 0.5, 0.0);

        prg_build_density_T0(norb, ham, rho, threshold, bndfil, eigenvalues);
        bml_scale(&scale_factor, rho, rho);

        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_build_density_T0: %.15e\n",
                 idempotency);
        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_density_T_c") == 0)
    {
        LOG_INFO
            ("Testing the construction of the density matrix from density_mod \n");

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");
        prg_build_density_T(norb, ham, rho, threshold, bndfil, kbt, mu,
                            eigenvalues);
        bml_scale(&scale_factor, rho, rho);

        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_build_density_T: %.15e\n", idempotency);

        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_density_T_fermi_c") == 0)
    {
        LOG_INFO
            ("Testing the construction of the density matrix at KbT > 0 and at mu = Ef from density_mod \n");

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        double drho;
        double ef = -0.10682896819759;

        prg_build_density_T_fermi(ham, rho, threshold, kbt, ef, 1, drho);
        bml_scale(&scale_factor, rho, rho);

        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_build_density_T_fermi: %.15e\n", idempotency);
        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }

    }
    else if (strcmp(test, "prg_density_T_fermi_grad_c") == 0)
    {
        LOG_INFO
            ("Testing the occupation gradient w.r.t. mu at KbT > 0 and at mu = Ef from density_mod \n");

        double drho_ref = 2.39454635999e-4;
        double ef = -0.10682896819759;
        double drho;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");
        prg_build_density_T_fermi(ham, rho, threshold, kbt, ef, 1, drho);

        if (abs(drho - drho_ref) > 1.0e-5)
        {
            printf("Difference is too high %f  %f\n", drho, drho_ref);
            exit(EXIT_FAILURE);
        }
    }

    //Diagonalize H and build \rho gradient w.r.t chemical potential mu
    else if (strcmp(test, "prg_density_cheb_fermi_c") == 0)
    {
        LOG_INFO
            ("Testing the construction of the density matrix at KbT > 0 and at mu = Ef from chebyshev_mod \n");

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        double ef = -0.10682896819759;
        double drho;
        double athr = 1.0;
        double fermitol = 0.001;
        int ncoeffs = 200;
        int npts = 2000;

        rho1 =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        prg_build_density_T_fermi(ham, rho, threshold, kbt, ef, 1, drho);
        prg_build_density_cheb_fermi(ham, rho1, athr, threshold, ncoeffs, kbt,
                                     mu, bndfil, 1, fermitol, 1, npts, 0,
                                     verbose);

        bml_add(rho1, rho, 1.0, -1.0, 0.0);
        error_calc = bml_fnorm(rho1);
        if (error_calc > 0.1)
        {
            printf("Error in Chebyshev expansion = %f\n", error_calc);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_implicit_fermi_c") == 0)
    {
        LOG_INFO
            ("Testing the construction of the density matrix at KbT > 0 and at mu = Ef from implicit_fermi_mod \n");
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        mu = 0.2;
        beta = 4.0;

        rho1 =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);

        prg_implicit_fermi(ham, rho1, 10, 2, 10.0, mu, beta, 0, 1, 1.0, threshold, 10e-8);
        prg_test_density_matrix(ham, rho, beta, mu, 10.0, 1, 1.0, threshold);

        bml_add(rho1, rho, 1.0, -1.0, threshold);

        error_calc = bml_fnorm(rho1);
        if (error_calc > 0.1)
        {
            printf("Error in Implicit Fermi expansion = %f", error_calc);
            exit(EXIT_FAILURE);
        }

    }

    else if (strcmp(test, "prg_implicit_fermi_save_inverse_c") == 0)
    {
        int norecs = 10;
        int occiter;
        int nsiter;
        double nocc = 10.0;

        mu = 0.2;
        beta = 4.0;

        for (i = 0; i < norecs; i++){
            inv_bml[i] =
                bml_identity_matrix(matrix_type, precision, norb, norb, distrib_mode);
        }
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        rho1 =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);

        prg_implicit_fermi_save_inverse(inv_bml, ham, rho, norecs, nocc,
                                        &mu, beta, 1e-4, threshold, 1e-5, 1, occiter, nsiter);

        prg_test_density_matrix(ham, rho1, beta, mu, nocc, 1, 1.0e-4, threshold);
        printf("mu= %f \n", mu);

        bml_scale(&scale_factor, rho, rho);
        bml_add(rho1, rho, 1.0, -1.0, threshold);

        error_calc = bml_fnorm(rho1);
        if (error_calc > 0.1)
        {
            printf("Error in Implicit Fermi expansion save inverse= %f \n", error_calc);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_basic_c") == 0)
    {
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");
        prg_sp2_basic(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                      sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_basic: %.15e\n", idempotency);

        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", &idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg1_dense_c") == 0)
    {
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");
        prg_sp2_alg1(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                     sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_dense: %.15e\n", idempotency);
        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_dense_c") == 0)
    {
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");
        prg_sp2_alg2(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                     sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_dense: %.15e\n", idempotency);
        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg1_ellpack_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-6;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        bml_read_bml_matrix(ham, "poly.512.mtx");
        prg_sp2_alg1(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                     sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_ellpack: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_ellpack_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-6;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        bml_read_bml_matrix(ham, "poly.512.mtx");
        prg_sp2_alg2(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                     sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_ellpack: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_ellpack_poly_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-2;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-5;
        sp2tol = 1.0e-7;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        bml_read_bml_matrix(ham, "poly.512.mtx");

        prg_sp2_alg2(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_ellpack_poly: %.15e\n",
                 idempotency);

        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg1_seq_dense_c") == 0)
    {
        matrix_type = dense;
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        int icount = 0;

        prg_sp2_alg1_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_seq_dense: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_seq_dense_c") == 0)
    {
        matrix_type = dense;
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        int icount = 0;

        prg_sp2_alg2_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv, verbose);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_seq_dense: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg1_seq_ellpack_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-6;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        int icount = 0;

        prg_sp2_alg1_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_seq_ellpack: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_seq_ellpack_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-6;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        int icount = 0;

        prg_sp2_alg2_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv, verbose);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_seq_ellpack: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg1_seq_inplace_dense_c") == 0)
    {
        matrix_type = dense;
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        double* gbnd = NULL; //malloc(2 * sizeof(double));
        int icount = 0;

        prg_sp2_alg1_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv);

        bml_copy(ham, rho);
        gbnd = bml_gershgorin(rho);
        prg_prg_sp2_alg1_seq_inplace(rho, threshold, maxsp2iter, pp, &icount,
                                     vv, gbnd[0], gbnd[1]);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_seq_inplace_dense: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_seq_inplace_dense_c") == 0)
    {
        matrix_type = dense;
        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        double *gbnd = malloc(2 * sizeof(double));
        //double* gbnd = NULL; //malloc(2 * sizeof(double));
        int icount = 0;

        prg_sp2_alg2_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv, verbose);

        bml_copy(ham, rho);
        gbnd = bml_gershgorin(rho);
        prg_prg_sp2_alg2_seq_inplace(rho, threshold, maxsp2iter, pp, &icount,
                                     vv, gbnd[0], gbnd[1]);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_seq_inplace_dense: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg1_seq_inplace_ellpack_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-6;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        double *gbnd = malloc(2 * sizeof(double));
        int icount = 0;

        prg_sp2_alg1_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv);

        bml_copy(ham, rho);
        gbnd = bml_gershgorin(rho);
        prg_prg_sp2_alg1_seq_inplace(rho, threshold, maxsp2iter, pp, &icount,
                                     vv, gbnd[0], gbnd[1]);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_seq_inplace_ellpack: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_alg2_seq_inplace_ellpack_c") == 0)
    {
        matrix_type = ellpack;
        idempotency_tol = 1.0e-6;
        bndfil = 0.5;
        norb = 6144;
        mdim = 600;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        int *pp = malloc(maxsp2iter * sizeof(int));
        double *vv = malloc(maxsp2iter * sizeof(double));
        double *gbnd = malloc(2 * sizeof(double));
        int icount = 0;

        prg_sp2_alg2_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, &icount, vv, verbose);

        bml_copy(ham, rho);
        gbnd = bml_gershgorin(rho);
        prg_prg_sp2_alg2_seq_inplace(rho, threshold, maxsp2iter, pp, &icount,
                                     vv, gbnd[0], gbnd[1]);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, &idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg2_seq_inplace_ellpack: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }

    else if (strcmp(test, "prg_sp2_fermi_dense_c") == 0)
    {
        matrix_type = dense;
        mdim = -1;
        threshold = 1.0e-9;
        sp2tol = 1.0e-10;
        occerrlimit = 1.0e-9;
        tracelimit = 1.0e-12;
        mineval = 0.0;
        maxeval = 0.0;
        beta0 = 20.1;
        tscale = 1.0;

        int norecs = 30;
        int occsteps = 0;
        double drho;
        double eps = 1.0e-4;
        double nocc = bndfil * norb;
        int *signlist = malloc(norecs * sizeof(int));
        memset(signlist, 0, norecs * sizeof(int));

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        bml_read_bml_matrix(ham, "hamiltonian_ortho.mtx");

        orthox_bml =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);

        prg_sp2_fermi_init_norecs(ham, &norecs, nocc, tscale, threshold,
                                  occerrlimit, tracelimit, orthox_bml, &mu,
                                  &beta0, &mineval, &maxeval, &signlist[0], 10);

        prg_sp2_fermi_init(ham, norecs, nocc, tscale, threshold, occerrlimit,
                           tracelimit, orthox_bml, &mu, &beta0, &mineval, &maxeval,
                           &signlist[0]);

        prg_sp2_fermi(ham, occsteps, norecs, nocc, &mu, &beta0, &mineval, &maxeval,
                      &signlist[0], threshold, eps, tracelimit, orthox_bml);

        LOG_INFO("BETA %.15e\n", beta0);
        LOG_INFO("ETEMP (eV) %.15e\n", 1.0/beta0);
        LOG_INFO("NORECS USED %10d\n", norecs);
        LOG_INFO("CHEMPOT %.15e\n", mu);

        kbt = 1.0 / beta0;

        prg_build_density_T_fermi(ham, rho, threshold, kbt, mu, 0, drho);
        bml_add(orthox_bml, rho, 2.0, -1.0, threshold);
        error_calc = bml_fnorm(orthox_bml);

        if (error_calc > 0.01)
        {
            printf("Error in sp2 Fermi %f is too high\n", error_calc);
            exit(EXIT_FAILURE);
        }

    }

    // more tests TBA
    else if(strcmp(test, "prg_equal_partition_c") == 0){
        //Create equal partitions
        //prg_equalPartition(gp, 6, 72)
    }

    else if(strcmp(test, "prg_file_partition_c") == 0){
        // Create partition from a file
    }

    else if(strcmp(test, "prg_subgraphsp2_equal_c") == 0){
        // Subgraph SP2 using equal size parts
    }

    else if(strcmp(test, "prg_deorthogonalize_dense_c") == 0){
        //Deorthogonalization of the density matrix
    }

    else if(strcmp(test, "prg_orthogonalize_dense_c") == 0){
        // Orthogonalization of the Hamiltonian
    }

    else if(strcmp(test, "prg_pulaycomponent0_c") == 0){
        error_tol = 1.e-9;
        matrix_type = dense;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        zmat =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        pcm =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);
        pcm_ref =
            bml_zero_matrix(matrix_type, precision, norb, norb, distrib_mode);

        bml_read_bml_matrix(zmat,"zmatrix.mtx");
        bml_read_bml_matrix(ham,"hamiltonian.mtx");
        bml_read_bml_matrix(rho,"density.mtx");
        bml_read_bml_matrix(pcm_ref,"pcm.mtx");

        prg_PulayComponent0(rho, ham, pcm, threshold, mdim, verbose);
        bml_add(pcm, pcm_ref, 1.0, -1.0, threshold);
        error_calc = bml_fnorm(pcm);
        LOG_INFO("Pulay Component error %.12e\n", error_calc);

        if (error_calc > error_tol)
        {
            printf("Error in PulayComponent0 %f is too high\n", error_calc);
            exit(EXIT_FAILURE);
        }

    }

    else if(strcmp(test, "prg_buildzdiag_c") == 0){
        // Building inverse overlap factor matrix (Lowdin method)
    }

    else if(strcmp(test, "prg_buildzsparse_c") == 0){
        // Building inverse overlap factor matrix (Lowdin method)
    }

    else if(strcmp(test, "prg_system_parse_write_xyz_c") == 0){
    }

    else if(strcmp(test, "prg_system_parse_write_pdb_c") == 0){
    }

    else if(strcmp(test, "prg_system_parse_write_dat_c") == 0){
    }

    else if(strcmp(test, "prg_twolevel_model_c") == 0){
    }

    else if(strcmp(test, "canon_response_c") == 0){
    }

    else if(strcmp(test, "prg_build_zmatGP_c") == 0){
    }

    else if(strcmp(test, "load_tbparms_latte_c") == 0){
    }

    else if(strcmp(test, "load_bintTBparamsH_c") == 0){
    }

    else
    {
        LOG_INFO("ERROR: unknown test \n");
        exit(EXIT_FAILURE);
    }

    prg_progress_shutdown();
    free(eigenvalues);

    return 0;
}
