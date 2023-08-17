
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
    double drho);

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
    double* idempotency);

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

//------prg_implicit_fermi_mod headers----------------------------------

    void prg_implicit_fermi_save_inverse(
    bml_matrix_t** Inv_bml,
    bml_matrix_t* h_bml,
    bml_matrix_t* p_bml,
    int nsteps,
    double nocc,
    double* mu,
    double beta,
    double occErrLimit,
    double threshold,
    double tol,
    int SCF_IT,
    int occiter,
    int totns);

    void prg_implicit_fermi(
    bml_matrix_t* h_bml,
    bml_matrix_t* p_bml,
    int nsteps,
    int k,
    double nocc,
    double mu,
    double beta,
    int method,
    int osteps,
    double occErrLimit,
    double threshold,
    double tol);

    void prg_implicit_fermi_zero(
    bml_matrix_t* h_bml,
    bml_matrix_t* p_bml,
    int nsteps,
    double mu,
    int method,
    double threshold,
    double tol);

    void prg_implicit_fermi_first_order_response(
    bml_matrix_t* H0_bml,
    bml_matrix_t* H1_bml,
    bml_matrix_t* P0_bml,
    bml_matrix_t* P1_bml,
    bml_matrix_t** Inv_bml,
    int nsteps,
    double mu0,
    double beta,
    double nocc,
    double threshold);

    void prg_implicit_fermi_response(
    bml_matrix_t* H0_bml,
    bml_matrix_t* H1_bml,
    bml_matrix_t* H2_bml,
    bml_matrix_t* H3_bml,
    bml_matrix_t* P0_bml,
    bml_matrix_t* P1_bml,
    bml_matrix_t* P2_bml,
    bml_matrix_t* P3_bml,
    int nsteps,
    double mu0,
    double mu,
    double beta,
    double nocc,
    double occ_tol,
    double lin_tol,
    int order,
    double threshold);

    void prg_finite_diff(
    bml_matrix_t* H0_bml,
    bml_matrix_t* H_list,
    double mu0,
    double mu_list,
    double beta,
    int order,
    double lambda,
    double h,
    double threshold);

    void prg_test_density_matrix(
    bml_matrix_t* ham_bml,
    bml_matrix_t* p_bml,
    double beta,
    double mu,
    double nocc,
    int osteps,
    double occErrLimit,
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
    int* pp,
    int* icount,
    double* vv,
    int verbose);

    void prg_sp2_alg2_seq(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    int* pp,
    int* icount,
    double* vv,
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

    void prg_sp2_alg1_genseq(
    bml_matrix_t* h_bml,
    bml_matrix_t* rho_bml,
    double threshold,
    double bndfil,
    int minsp2iter,
    int maxsp2iter,
    char* sp2conv,
    double idemtol,
    int* pp,
    int* icount,
    double* vv);

    void prg_sp2_alg1_seq(
    bml_matrix_t * h_bml,
    bml_matrix_t * rho_bml,
    double threshold,
    int* pp,
    int* icount,
    double* vv);

    void prg_prg_sp2_alg2_seq_inplace(
    bml_matrix_t * rho_bml,
    double threshold,
    int ppsize,
    int* pp,
    int* icount,
    double* vv,
    double mineval,
    double maxeval);
    //int verbose);

    void prg_prg_sp2_alg1_seq_inplace(
    bml_matrix_t * rho_bml,
    double threshold,
    int ppsize,
    int* pp,
    int* icount,
    double* vv,
    double mineval,
    double maxeval);

    void prg_sp2_submatrix(
    double threshold,
    int* pp,
    int* icount,
    double* vv,
    double mineval,
    double maxeval,
    int core_size);

    void prg_sp2_submatrix_inplace(
    bml_matrix_t * rho_bml,
    double threshold,
    int* pp,
    int* icount,
    double* vv,
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

     void prg_get_pulayforce(
     int nats,
     bml_matrix_t* zmat_bml,
     bml_matrix_t* ham_bml,
     bml_matrix_t* rho_bml,
     bml_matrix_t* dSx_bml,
     bml_matrix_t* dSy_bml,
     bml_matrix_t* dSz_bml,
     int **hindex,
     double **FPUL,
     double threshold);


//-----prg_genz_mod headers -------------------------------------------
    //void prg_parse_ZSP(char* filename);

    void prg_init_ZSPmat(
    int igenz, bml_matrix_t* zk1_bml, bml_matrix_t* zk2_bml, bml_matrix_t* zk3_bml, bml_matrix_t* zk4_bml, bml_matrix_t* zk5_bml, bml_matrix_t* zk6_bml, int norb, char* bml_type, char* bml_element_type);

    void prg_buildZdiag(bml_matrix_t* smat_bml, bml_matrix_t* zmat_bml, double threshold, int mdimin, char* bml_type, int verbose);

    void prg_genz_sp_initialz0(bml_matrix_t* smat_bml, bml_matrix_t* zmat_bml, int norb, int mdim, char* bml_type_f, double threshold);

    void prg_genz_sp_initial_zmat(bml_matrix_t* smat_bml, bml_matrix_t* zmat_bml, int norb, int mdim, char* bml_type_f, double threshold);

    void prg_genz_sp_ref(bml_matrix_t* smat_bml, bml_matrix_t* zmat_bml, int nref, int norb, char* bml_type, double threshold);

// prg_parallel_mod

    void prg_initParallel();

    void prg_shutdownParallel();

    void prg_barrierParallel();

    void sendReceiveParallel(double sendBuf, int sendLen, int dest, double recvBuf, int recvLen, int source);

    void isendParallel(double* sendBuf, int sendLen, int dest);

    void sendParallel(double* sendBuf, int sendLen, int dest);

    void prg_iprg_recvParallel(double* recvBuf, int recvLen, int rind);

    void prg_recvParallel(double* recvBuf, int recvLen);

    void sumIntParallel(int* sendBuf, int* recvBuf, int icount);

    void sumRealParallel(double* sendBuf, double* recvBuf, int icount);

    void maxIntParallel(int* sendBuf, int* recvBuf, int icount);

    void maxRealParallel(double* sendBuf, double* recvBuf, int icount);

    void minIntParallel(double* sendBuf, double* recvBuf, int icount);

    void minRealParallel(double* sendBuf, double* recvBuf, int icount);

    void prg_minRealReduce(double rvalue);

    void prg_maxRealReduce(double rvalue);

    void prg_maxIntReduce2(int value1, int value2);

    void prg_sumIntReduce2(int value1, int value2);

    //void prg_sumRealReduce();

//void prg_sumRealReduce2();
//
//void prg_sumRealReduce3();
//
//void prg_sumRealReduceN(int N);
//
//void prg_sumIntReduceN(int valueVec, int N);
//
//void minRankRealParallel(int icount);
//
//void maxRankRealParallel(int icount);
//
//void prg_bcastParallel(int blen, int root);
//
//void allGatherRealParallel(double sendBuf, int sendLen, double recvBuf, int recvLen);
//
//void allGatherIntParallel(int sendBuf, int sendLen, int recvBuf, int recvLen);
//
//void allGatherVRealParallel(double sendBuf, int sendLen, double recvBuf, int recvLen, int recvDispl);
//
//void allGatherVIntParallel(int sendBuf, int sendLen, int recvBuf, int recvLen, int recvDispl);
//
//void prg_allSumRealReduceParallel(double buf, int buflen);
//
//void prg_allSumIntReduceParallel(int buf, int buflen);
//
//void prg_allGatherParallel();

    void prg_wait();

//-----prg_graph_mod headers (TBA)-------------------------------------------

//-----prg_graphsolver_mod headers (TBA)-------------------------------------
    void prg_build_densityGP_T0(
    bml_matrix_t* ham_bml,
    bml_matrix_t* g_bml,
    bml_matrix_t* rho_bml,
    double threshold,
    double bndfil,
    double Ef,
    int nparts,
    int verbose);

    void prg_build_zmatGP(
    bml_matrix_t* over_bml,
    bml_matrix_t* g_bml,
    bml_matrix_t* zmat_bml,
    double threshold,
    int nparts,
    int verbose);



#ifdef __cplusplus
}
#endif

#endif                          /* PROGRESS_INTERFACE_H */
