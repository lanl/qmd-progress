
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

void prg_build_density_T0(bml_matrix_t* ham_bml, bml_matrix_t* rho_bml, double threshold, double bndfil);

// old interface

//void prg_timer_start_c(int* itimer, char* tag);
//void prg_timer_stop_c(int* itimer, char* verbose);

#ifdef __cplusplus
}
#endif

#endif /* PROGRESS_INTERFACE_H */

