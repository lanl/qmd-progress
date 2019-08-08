!> Standalone graph partitioning test program.
!!
program gptest

  use bml 
  use prg_progress_mod
  use prg_parallel_mod
  use prg_timer_mod
  use prg_graphsp2parser_mod
  use prg_sp2_mod
  use prg_graph_mod
  use prg_subgraphLoop_mod
  use prg_homolumo_mod

  implicit none     

  integer, parameter :: dp = kind(1.0d0)

  integer ::  norb
  type(gsp2data_type) :: gsp2
  type(bml_matrix_t) :: h_bml, rho_bml, g_bml, fullg_bml
  type(graph_partitioning_t) :: gp

  integer :: N, M, i, icount
  integer :: pp(100)
  real(dp) :: vv(100)
  real(dp) :: mineval, maxeval, trace_mult, trace_multg
  real(dp) :: ehomo, elumo, egap, frobNorm
  real(dp), allocatable :: gbnd(:)

  !< Start progress
  call prg_progress_init()

  if (printRank() .eq. 1) then
      write(*,*) "gptest start ..."
  endif

  !< parsing input file.
  call prg_parse_gsp2(gsp2,"input.in")

  !< Allocate hamiltonian and density matrices
  call bml_zero_matrix(gsp2%bml_type, BML_ELEMENT_REAL, dp, gsp2%ndim, gsp2%mdim, h_bml)
  call bml_zero_matrix(gsp2%bml_type, BML_ELEMENT_REAL, dp, gsp2%ndim, gsp2%mdim, rho_bml)
  call bml_read_matrix(h_bml, gsp2%hamfile)

  N = bml_get_N(h_bml)
  M = bml_get_M(h_bml)

  !< Create the same size matrix for the connectivity graph
  call bml_zero_matrix(BML_MATRIX_ELLPACK, BML_ELEMENT_REAL, dp, N, M, g_bml);

  !< Run regular SP2 
  call prg_timer_start(sp2_timer)
  call prg_sp2_alg2_genseq(h_bml, g_bml, gsp2%threshold, gsp2%bndfil, gsp2%minsp2iter, &
         gsp2%maxsp2iter, gsp2%sp2conv, gsp2%sp2tol, pp, icount, vv)
  call prg_timer_stop(sp2_timer)

  !< List fnorm per iteration
  if (printRank() .eq. 1) then
    write(*,*)
    do i = 1, icount
      write(*,*) "iter = ", i, " vv = ", vv(i)
    enddo
    write(*,*)
  endif

  !< Calculate Trace[HP]
  trace_mult = bml_trace_mult(h_bml, g_bml)
#ifdef DO_MPI_BLOCK
  if (getNRanks() > 1) then
    call prg_sumRealReduce(trace_mult)
  endif
#endif
  if (printRank() .eq. 1) then
    write(*,*) "Trace[HP] for SP2 = ", trace_mult
    write(*,*) "Band energy per atom = ", trace_mult/gsp2%natoms
    write(*,*)
  endif

  !< Calculate Homo-Lumo gap
  allocate(gbnd(2))
  call bml_gershgorin(h_bml, gbnd)
  mineval = gbnd(1)
  maxeval = gbnd(2)
  deallocate(gbnd)
  call prg_homolumogap(vv, icount, pp, mineval, maxeval, ehomo, elumo, egap)
  if (printRank() .eq. 1) then
    write(*,*) "Gershgorin: mineval = ", mineval, " maxeval = ", maxeval
    write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
         " egap = ", egap
    write(*,*)
  endif

  !< Save copy of SP2 density matrix
  call bml_copy_new(g_bml, fullg_bml)

  !< Start graph partitioning SP2 

  !< Threshold the graph
  call bml_threshold(g_bml, gsp2%gthreshold)

  !< Create graph partitioning of equal parts
  call prg_equalPartition(gp, gsp2%nodesPerPart, N)
  gp%mineval = mineval
  gp%maxeval = maxeval

  !< Calculate SP2 sequence
  call prg_sp2sequence(gp%pp, gp%maxIter, gp%mineval, gp%maxeval, ehomo, elumo, &
        gsp2%errlimit)
  if (printRank() .eq. 1) then
    write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter
    write(*,*)
  endif

  !< Run subgraph SP2 loop
  call prg_timer_start(graphsp2_timer)
  call prg_subgraphSP2Loop(h_bml, g_bml, rho_bml, gp, gsp2%threshold)
  call prg_timer_stop(graphsp2_timer)

  !< Calculate homo-lumo gap
  call prg_homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, &
         egap)
  if (printRank() .eq. 1) then
    write(*,*)
    write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
         " egap = ", egap
    write(*,*)
  endif

  !< Calculate trace[HP]
  trace_multg = bml_trace_mult(h_bml, rho_bml)
#ifdef DO_MPI_BLOCK
  if (getNRanks() > 1) then
    call prg_sumRealReduce(trace_multg)
  endif
#endif
  if (printRank() .eq. 1) then
    write(*,*) "Trace[HP] for subgraph SP2 = ", trace_multg
    write(*,*) "Band energy per atom = ", trace_multg/gsp2%natoms
    write(*,*)
  endif

  !< Calculate frobenius norm
  frobNorm = bml_fnorm2(fullg_bml, rho_bml)
  if (printRank() .eq. 1) then
    write(*,*) "Frobenius norm = ", frobNorm
    write(*,*) "Frobenius norm/atom = ", frobNorm/gsp2%natoms
    write(*,*) "Error in band energy = ", &
        trace_mult/gsp2%natoms - trace_multg/gsp2%natoms
    write(*,*)
  endif

  call bml_deallocate(fullg_bml)
  call bml_deallocate(g_bml)
  call bml_deallocate(h_bml)
  call bml_deallocate(rho_bml)

  call prg_destroyGraphPartitioning(gp)

  call prg_progress_shutdown()

  call exit(0)

end program gptest
