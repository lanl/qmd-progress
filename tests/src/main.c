
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
    bml_matrix_t *zmat_bml, *zk1_bml, *zk2_bml;
    bml_matrix_t *zk3_bml, *zk4_bml, *zk5_bml, *zk6_bml;
    bml_matrix_t *nonortho_ham_bml;
    bml_matrix_t *over_bml, *orthox_bml, *g_bml;
    bml_matrix_t *aux_bml, *pcm_bml, *pcm_ref_bml, *inv_bml[10];

    bml_matrix_type_t matrix_type;
    bml_matrix_precision_t precision;
    bml_distribution_mode_t distrib_mode = sequential;

    char bml_type[21];          // 20 char + null termination '\0'
    char sp2conv_sue[21];       // 20 char + null termination '\0'
    char dummy[10][21];         // 20 char + null termination '\0' each
    char sp2conv[10] = "Rel";

    double threshold, gthreshold, idempotency;
    double sp2tol, idempotency_tol;
    double bndfil, tscale, tracelimit, beta;
    double error_calc, error_tol, errlimit;

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

        prg_check_idempotency(rho, threshold, idempotency);
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

        prg_check_idempotency(rho, threshold, idempotency);
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

    // more tests TBA
    else if (strcmp(test, "prg_implicit_fermi_c") == 0)
    {
        //TBA
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
        //call prg_implicit_fermi(ham, rho1, 10, 2, 10.0_dp, mu, beta, 0, 1, 1.0_dp, threshold, 10e-8_dp)
        //call prg_test_density_matrix(ham, rho, beta, mu, 10.0_dp, 1, 1.0_dp, threshold)
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
        prg_check_idempotency(rho, threshold, idempotency);
        LOG_INFO("Idempotency for prg_sp2_basic: %.15e\n", idempotency);

        if (idempotency > 1.0e-5)
        {
            printf("Idempotency is too high %f\n", idempotency);
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
        prg_check_idempotency(rho, threshold, idempotency);
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
        prg_check_idempotency(rho, threshold, idempotency);
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
        prg_check_idempotency(rho, threshold, idempotency);
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
        prg_check_idempotency(rho, threshold, idempotency);
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
        mdim = 288;
        threshold = 1.0e-5;
        sp2tol = 1.0e-10;

        rho =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        ham =
            bml_zero_matrix(matrix_type, precision, norb, mdim, distrib_mode);
        bml_read_bml_matrix(ham, "poly.512.mtx");

        prg_sp2_alg2(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, sp2tol, verbose);
        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, idempotency);
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

        double *pp = malloc(100 * sizeof(double));
        double *vv = malloc(100 * sizeof(double));
        int icount = 0;

        prg_sp2_alg1_genseq(ham, rho, threshold, bndfil, minsp2iter, maxsp2iter,
                            sp2conv, sp2tol, pp, icount, vv);

        bml_scale(&scale_factor, rho, rho);
        prg_check_idempotency(rho, threshold, idempotency);
        LOG_INFO("Idempotency for prg_sp2_alg1_seq_dense: %.15e\n",
                 idempotency);
        if (idempotency > idempotency_tol)
        {
            printf("Idempotency is too high %f\n", idempotency);
            exit(EXIT_FAILURE);
        }
    }


    //else if(strcmp(test, "prg_xxxx_c") == 0){
    //}
    else
    {
        LOG_INFO("ERROR: unknown test \n");
        exit(EXIT_FAILURE);
    }

    prg_progress_shutdown();
    free(eigenvalues);

    return 0;
}
