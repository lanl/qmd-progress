!> Standalone graph partitioning test program.
!!
program gptest

  use bml 
  use progress_mod
  use parallel_mod
  use timer_mod
  use graph_sp2parser_mod
  use sp2_mod
  use graph_mod
  use subgraphLoop_mod
  use homolumo_mod

  implicit none     

  integer, parameter :: dp = kind(1.0d0)

  integer ::  norb
  type(gsp2data_type) :: gsp2
  type(bml_matrix_t) :: h_bml, rho_bml, g_bml, fullg_bml
  type(graph_partitioning_t) :: gp

  integer :: N, M, i, icount
  integer :: pp(100)
  real(dp) :: vv(100)
  real(dp) :: mineval, maxeval, traceMult, traceMultg
  real(dp) :: ehomo, elumo, egap, frobNorm
  real(dp), allocatable :: gbnd(:)

  !< Start progress
  call progress_init()

  if (printRank() .eq. 1) then
      write(*,*) "gptest start ..."
  endif

  !< parsing input file.
  call parse_gsp2(gsp2,"input.in")

  !< Allocate hamiltonian and density matrices
  call bml_zero_matrix(gsp2%bml_type, BML_ELEMENT_REAL, dp, gsp2%ndim, gsp2%mdim, h_bml)
  call bml_zero_matrix(gsp2%bml_type, BML_ELEMENT_REAL, dp, gsp2%ndim, gsp2%mdim, rho_bml)
  call bml_read_matrix(h_bml, gsp2%hamfile)

  N = bml_get_N(h_bml)
  M = bml_get_M(h_bml)

  !< Create the same size matrix for the connectivity graph
  call bml_zero_matrix(BML_MATRIX_ELLPACK, BML_ELEMENT_REAL, dp, N, M, g_bml);

  !< Run regular SP2 
  call timer_start(sp2_timer)
  call sp2_alg2_genseq(h_bml, g_bml, gsp2%threshold, gsp2%bndfil, gsp2%minsp2iter, &
         gsp2%maxsp2iter, gsp2%sp2conv, gsp2%sp2tol, pp, icount, vv)
  call timer_stop(sp2_timer)

  !< List fnorm per iteration
  if (printRank() .eq. 1) then
    write(*,*)
    do i = 1, icount
      write(*,*) "iter = ", i, " vv = ", vv(i)
    enddo
    write(*,*)
  endif

  !< Calculate Trace[HP]
  traceMult = bml_traceMult(h_bml, g_bml)
#ifdef DO_MPI_BLOCK
  if (getNRanks() > 1) then
    call sumRealReduce(traceMult)
  endif
#endif
  if (printRank() .eq. 1) then
    write(*,*) "Trace[HP] for SP2 = ", traceMult
    write(*,*) "Band energy per atom = ", traceMult/gsp2%natoms
    write(*,*)
  endif

  !< Calculate Homo-Lumo gap
  allocate(gbnd(2))
  call bml_gershgorin(h_bml, gbnd)
  mineval = gbnd(1)
  maxeval = gbnd(2)
  deallocate(gbnd)
  call homolumogap(vv, icount, pp, mineval, maxeval, ehomo, elumo, egap)
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
  call equalPartition(gp, gsp2%nodesPerPart, N)
  gp%mineval = mineval
  gp%maxeval = maxeval

  !< Calculate SP2 sequence
  call sp2sequence(gp%pp, gp%maxIter, gp%mineval, gp%maxeval, ehomo, elumo, &
        gsp2%errlimit)
  if (printRank() .eq. 1) then
    write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter
    write(*,*)
  endif

  !< Run subgraph SP2 loop
  call timer_start(graphsp2_timer)
  call subgraphSP2Loop(h_bml, g_bml, rho_bml, gp, gsp2%threshold)
  call timer_stop(graphsp2_timer)

  !< Calculate homo-lumo gap
  call homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, &
         egap)
  if (printRank() .eq. 1) then
    write(*,*)
    write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
         " egap = ", egap
    write(*,*)
  endif

  !< Calculate trace[HP]
  traceMultg = bml_traceMult(h_bml, rho_bml)
#ifdef DO_MPI_BLOCK
  if (getNRanks() > 1) then
    call sumRealReduce(traceMultg)
  endif
#endif
  if (printRank() .eq. 1) then
    write(*,*) "Trace[HP] for subgraph SP2 = ", traceMultg
    write(*,*) "Band energy per atom = ", traceMultg/gsp2%natoms
    write(*,*)
  endif

  !< Calculate frobenius norm
  frobNorm = bml_fnorm2(fullg_bml, rho_bml)
  if (printRank() .eq. 1) then
    write(*,*) "Frobenius norm = ", frobNorm
    write(*,*) "Frobenius norm/atom = ", frobNorm/gsp2%natoms
    write(*,*) "Error in band energy = ", &
        traceMult/gsp2%natoms - traceMultg/gsp2%natoms
    write(*,*)
  endif

  call bml_deallocate(fullg_bml)
  call bml_deallocate(g_bml)
  call bml_deallocate(h_bml)
  call bml_deallocate(rho_bml)

  call destroyGraphPartitioning(gp)

  call progress_shutdown()

  call exit(0)

end program gptest
