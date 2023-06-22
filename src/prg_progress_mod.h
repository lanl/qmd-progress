
#ifndef PROGRESS_INTERFACE_H
#define PROGRESS_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <bml.h>

extern void prg_version();
extern void prg_progress_init();
extern void prg_progress_shutdown();

//-------- prg_densitymatrix_mod headers ------------------------------------
void prg_build_density_T0(int norbs, bml_matrix_t* ham_bml, bml_matrix_t* rho_bml, double threshold, double bndfil, double* eigenvalues_out);

void prg_build_density_T(int norbs, bml_matrix_t* ham_bml, bml_matrix_t* rho_bml, double threshold, double bndfil, double kbt, double ef, double* eigenvalues_out);

void prg_build_density_T_fulldata(int norbs, bml_matrix_t* ham_bml, bml_matrix_t* rho_bml, double threshold, double bndfil, double kbt, double ef, double* eigenvalues_out, bml_matrix_t* evects_bml, double* fvals_out);

void prg_build_density_T_ed(int norbs, bml_matrix_t* ham_bml, bml_matrix_t* rho_bml, bml_matrix_t* evects_bml, double threshold, double bndfil, double kbt, double ef, double* evals, double* dvals, int** hindex,int llsize, int verbose);

void prg_get_evalsDvalsEvects(int norbs, bml_matrix_t* ham_bml, double threshold, int** hindex, int llsize, double* evals, double* dvals, bml_matrix_t* evects_bml, int verbose);

void prg_build_density_fromEvalsAndEvects(int norbs, bml_matrix_t* evects_bml, bml_matrix_t* rho_bml, double threshold, double bndfil, double kbt, double ef, int verbose);

void prg_build_density_T_Fermi(bml_matrix_t* ham_bml, bml_matrix_t* rho_bml, double threshold, double kbt, double ef, int verbose, double drho);

void prg_build_atomic_density(bml_matrix_t* rhoat_bml, int norb, char* bml_type);

void prg_get_flevel(double kbt, double bndfil, double tol, double Ef, int* err);

void prg_get_flevel_nt(double kbt, double bndfil, double tol, double ef, int* err, int verbose);

void prg_get_eigenvalues(int norbs, bml_matrix_t* ham_bml, double* eigenvalues, int verbose);

void prg_check_idempotency(bml_matrix_t* mat_bml, double threshold, double idempotency);

void prg_toEigenspace(bml_matrix_t* mat_bml, bml_matrix_t* matEig_bml, bml_matrix_t* evects_bml, double threshold, int verbose);

void prg_toCanonicalspace(bml_matrix_t* mat_bml, bml_matrix_t* matCan_bml, bml_matrix_t* evects_bml, double threshold, int verbose);

void Canon_DM_PRT(double T, double mu0, int m, int HDIM);void prg_write_tdos(double gamma, int npts, double emin, double emax);void Ewald_Real_Space_Single_latte(double COULOMBV, int I, int Nr_elem, int J, int Nr_atoms, double COULACC, int HDIM, int Max_Nr_Neigh);

//-----prg_charges_mod headers (TBA)-----------------------------------------

//-----prg_graph_mod headers (TBA)-------------------------------------------

//-----prg_graphsolver_mod headers (TBA)-------------------------------------



#ifdef __cplusplus
}
#endif

#endif /* PROGRESS_INTERFACE_H */

