
#include "prg_progress_mod.h"

int
main(
    )
{
    //prg_progress_init_c();   // initialize progress
    //prg_version_c();   // print version
    //prg_progress_shutdown_c();   // shutdown progress

    // new interface
    prg_progress_init();        // initialize progress
    prg_version();              // print version
    prg_progress_shutdown();    // shutdown progress

}
