
#ifndef PROGRESS_INTERFACE_H
#define PROGRESS_INTERFACE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>
#include <bml.h>

    extern void prg_version(
        );
    extern void prg_progress_init(
        );
    extern void prg_progress_shutdown(
        );

//-------- prg_densitymatrix_mod headers ------------------------------------
    void prg_build_density_T0(
    int norbs,
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    double *eigenvalues_out);

    void prg_build_density_T(
    int norbs,
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    double kbt,
    double ef,
    double *eigenvalues_out);

    void prg_build_density_T_fulldata(
    int norbs,
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    double kbt,
    double ef,
    double *eigenvalues_out,
    bml_matrix_t * evects_bml,
    double *fvals_out);

    void prg_build_density_T_ed(
    int norbs,
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    bml_matrix_t * evects_bml,
    double threshold,
    double bndfil,
    double kbt,
    double ef,
    double *evals,
    double *dvals,
    int **hindex,
    int llsize,
    int verbose);

    void prg_get_evalsDvalsEvects(
    int norbs,
    bml_matrix_t * ham_bml,
    double threshold,
    int **hindex,
    int llsize,
    double *evals,
    double *dvals,
    bml_matrix_t * evects_bml,
    int verbose);

    void prg_build_density_fromEvalsAndEvects(
    int norbs,
    bml_matrix_t * evects_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    double kbt,
    double ef,
    int verbose);

    void prg_build_density_T_fermi(
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double kbt,
    double ef,
    int verbose,
    double *drho);

    void prg_build_atomic_density(
    bml_matrix_t * rhoat_bml,
    int norb,
    char *bml_type);

    void prg_get_flevel(
    double kbt,
    double bndfil,
    double tol,
    double Ef,
    int *err);

    void prg_get_flevel_nt(
    double kbt,
    double bndfil,
    double tol,
    double ef,
    int *err,
    int verbose);

    void prg_get_eigenvalues(
    int norbs,
    bml_matrix_t * ham_bml,
    double *eigenvalues,
    int verbose);

    void prg_check_idempotency(
    bml_matrix_t * mat_bml,
    double threshold,
    double idempotency);

    void prg_toEigenspace(
    bml_matrix_t * mat_bml,
    bml_matrix_t * matEig_bml,
    bml_matrix_t * evects_bml,
    double threshold,
    int verbose);

    void prg_toCanonicalspace(
    bml_matrix_t * mat_bml,
    bml_matrix_t * matCan_bml,
    bml_matrix_t * evects_bml,
    double threshold,
    int verbose);

    void Canon_DM_PRT(
    double T,
    double mu0,
    int m,
    int HDIM);

//-----prg_charges_mod headers -----------------------------------------
    void prg_get_charges(
    int nats,
    int norbs,
    bml_matrix_t * rho_bml,
    bml_matrix_t * over_bml,
    int **hindex,
    double *charges,
    double *numel,
    int *spindex,
    int mdimin,
    double threshold);

    void prg_get_hscf(
    int nats,
    bml_matrix_t * ham0_bml,
    bml_matrix_t * over_bml,
    bml_matrix_t * ham_bml,
    int *spindex,
    int **hindex,
    double *hubbardu,
    double *charges,
    double *coulomb_pot_r,
    double *coulomb_pot_k,
    int mdimin,
    double threshold);

    void prg_get_hscf_v2(
    int nats,
    bml_matrix_t * ham0_bml,
    bml_matrix_t * over_bml,
    bml_matrix_t * ham_bml,
    int *spindex,
    int **hindex,
    double *hubbardu,
    double *charges,
    double *coulomb_pot_r,
    double *coulomb_pot_k,
    int mdimin,
    double threshold);

//------prg_chebyshev_mod headers --------------------------------------
//void prg_parse_cheb(char* filename);

    void prg_build_density_cheb(
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    double athr,
    double threshold,
    int ncoeffs,
    double kbt,
    double ef,
    double bndfil,
    int jon,
    int verbose);

    void prg_build_density_cheb_fermi(
    bml_matrix_t * ham_bml,
    bml_matrix_t * rho_bml,
    double athr,
    double threshold,
    int ncoeffs,
    double kbt,
    double ef,
    double bndfil,
    int getef,
    double fermitol,
    int jon,
    int npts,
    int trkfunc,
    int verbose);

//-----prg_dos_mod -----------------------------------------------------

    void prg_write_tdos(
    int nstates,
    double *eigenvals,
    double gamma,
    int npts,
    double emin,
    double emax,
    char *filename);

// prg_ewald_mod
    void Ewald_Real_Space_Single_latte(
    double COULOMBV,
    int I,
    double RXYZ,
    double Box,
    int Nr_elem,
    double DELTAQ,
    int J,
    double U,
    int Element_Pointer,
    int Nr_atoms,
    double COULACC,
    int HDIM,
    int Max_Nr_Neigh);

    void Ewald_Real_Space_Single(
    double COULOMBV,
    double FCOUL,
    int I,
    double RX,
    double RY,
    double RZ,
    double LBox,
    double DELTAQ,
    int J,
    double U,
    char *Element_Type,
    int Nr_atoms,
    double COULACC,
    double TIMERATIO,
    int HDIM,
    int Max_Nr_Neigh);

    void Ewald_Real_Space_Matrix_latte(
    double E,
    double RXYZ,
    double Box,
    double U,
    int Element_Pointer,
    int Nr_atoms,
    double COULACC,
    int nebcoul,
    int totnebcoul,
    int HDIM,
    int Max_Nr_Neigh,
    int Nr_Elem);

    void Ewald_Real_Space_latte(
    double COULOMBV,
    int I,
    double RXYZ,
    double Box,
    double DELTAQ,
    double U,
    int Element_Pointer,
    int Nr_atoms,
    double COULACC,
    int nebcoul,
    int totnebcoul,
    int HDIM,
    int Max_Nr_Neigh,
    int Nr_Elem);

    void Ewald_Real_Space_Test(
    double COULOMBV,
    int I,
    double RX,
    double RY,
    double RZ,
    double LBox,
    double DELTAQ,
    double U,
    char *Element_Type,
    int Nr_atoms,
    double COULACC,
    double nnRx,
    double nnRy,
    double nnRz,
    int nrnnlist,
    int nnType,
    int Max_Nr_Neigh);

    void Ewald_Real_Space(
    double COULOMBV,
    double FCOUL,
    int I,
    double RX,
    double RY,
    double RZ,
    double LBox,
    double DELTAQ,
    double U,
    char *Element_Type,
    int Nr_atoms,
    double COULACC,
    double TIMERATIO,
    double nnRx,
    double nnRy,
    double nnRz,
    int nrnnlist,
    int nnType,
    int HDIM,
    int Max_Nr_Neigh);
    void prg_print_matrix(
    char *matname,
    double amat,
    int i1,
    int i2,
    int j1,
    int j2);

// prg_sp2_fermi headers
    void prg_sp2_fermi_init_norecs(
    bml_matrix_t * h_bml,
    int nsteps,
    double nocc,
    double tscale,
    double threshold,
    double occErrLimit,
    double traceLimit,
    bml_matrix_t * x_bml,
    double mu,
    double beta,
    double h1,
    double hN,
    int sgnlist,
    int verbose);

    void prg_sp2_fermi(
    bml_matrix_t * h_bml,
    int osteps,
    int nsteps,
    double nocc,
    double mu,
    double beta,
    double h1,
    double hN,
    int sgnlist,
    double threshold,
    double eps,
    double traceLimit,
    bml_matrix_t * x_bml);

    void prg_sp2_entropy_function(
    double mu,
    double h1,
    double hN,
    int nsteps,
    int sgnlist,
    double GG,
    double ee);


// prg_sp2_mode headers
    void prg_sp2_basic(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    int minsp2iter,
    int maxsp2iter,
    char *sp2conv,
    double idemtol,
    int verbose);

    void prg_sp2_basic_tcore(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    bml_matrix_t * rhofull_bml,
    double threshold,
    double bndfil,
    int minsp2iter,
    int maxsp2iter,
    char *sp2conv,
    double idemtol,
    int verbose);

    void prg_sp2_alg2(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    int minsp2iter,
    int maxsp2iter,
    char *sp2conv,
    double idemtol,
    int verbose);

    void prg_sp2_alg2_genseq(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    int minsp2iter,
    int maxsp2iter,
    char *sp2conv,
    double idemtol,
    int pp,
    int icount,
    double vv,
    int verbose);

    void prg_sp2_alg2_seq(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    int pp,
    int icount,
    double vv,
    int verbose);

    void prg_prg_sp2_alg2_seq_inplace(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    int pp,
    int icount,
    double vv,
    double mineval,
    double maxeval,
    int verbose);

    void prg_sp2_alg1(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    double bndfil,
    int minsp2iter,
    int maxsp2iter,
    char *sp2conv,
    double idemtol,
    int verbose);

    void prg_sp2_alg1_seq(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    int pp,
    int icount,
    double vv);

    void prg_prg_sp2_alg1_seq_inplace(
    double threshold,
    int pp,
    int icount,
    double vv,
    double mineval,
    double maxeval);

    void prg_sp2_submatrix(
    double threshold,
    int pp,
    int icount,
    double vv,
    double mineval,
    double maxeval,
    int core_size);

    void prg_sp2_submatrix_inplace(
    double threshold,
    int pp,
    int icount,
    double vv,
    double mineval,
    double maxeval,
    int core_size);



// prg_timer_mod headers
    void timer_prg_init(
        );

    void prg_timer_shutdown(
        );

    void prg_timer_start(
    int itimer,
    char *tag);

    void prg_timer_stop(
    int itimer,
    int verbose);

    void prg_timer_collect(
        );

    void prg_timer_results(
        );

    void prg_print_date_and_time(
    char *tag);


#ifdef __cplusplus
}
#endif

#endif                          /* PROGRESS_INTERFACE_H */
