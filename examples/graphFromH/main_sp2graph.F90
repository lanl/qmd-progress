program main_sp2graph

  use prg_graph_mod
  use prg_sp2_mod
  use prg_subgraphloop_mod
  use prg_homolumo_mod
  use prg_timer_mod
  use bml
  use omp_lib
  use prg_modelham_mod
  use prg_progress_mod

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: nodesPerPart
  integer :: minsp2iter, maxsp2iter
  real(dp) :: sthreshold, gthreshold, bndfil, idemtol
  real(dp) :: errlimit, threshold
  character(len=3) :: sp2conv
  type(bml_matrix_t) :: h_bml
  type(bml_matrix_t) :: rho_bml, aux_bml
  type(mham_type) ::  mham
  type(bml_matrix_t) :: g_bml
  type(graph_partitioning_t) :: gp
  integer :: N, M, i, icount, verbose
  integer :: pp(100)
  real(dp) :: vv(100)
  real(dp) :: mineval, maxeval, traceMult
  real(dp) :: ehomo, elumo, egap
  real(dp), allocatable :: gbnd(:)


  call prg_progress_init()

  !> General parameters
  threshold = 1.0d-5
  bndfil = 0.5_dp
  gthreshold = 1.0d-3
  errlimit = 1.0d-12
  nodesPerPart = 100
  minsp2iter = 10
  maxsp2iter = 50
  sp2conv="REL"
  idemtol=1.0D-5
  bndfil = 0.5_dp
  sthreshold=1.0D-5


  !Parsing input file.
  call prg_parse_mham(mham,"input.in") !Reads the input for modelham

  !Constructng the Hamiltonian               
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,mham%norbs,mham%norbs,h_bml)
  call prg_twolevel_model(mham%ea, mham%eb, mham%dab, mham%daiaj, mham%dbibj, &
       &mham%dec, mham%rcoeff, mham%reshuffle, mham%seed, h_bml, verbose)

  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,mham%norbs,mham%norbs,rho_bml)

  !! Create the same size matrix for the connectivity graph
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,mham%norbs,mham%norbs,g_bml);
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,mham%norbs,mham%norbs,aux_bml);

  !! Run regular SP2
  call prg_timer_start(sp2_timer)
  call prg_sp2_alg2_genseq(h_bml, g_bml, sthreshold, bndfil, minsp2iter, &
       maxsp2iter, sp2conv, idemtol, pp, icount, vv)

  call bml_print_matrix("g_bml",g_bml,0,10,0,10)
  call prg_timer_stop(sp2_timer)

  !! Calculate Homo-Lumo gap
  allocate(gbnd(2))
  call bml_gershgorin(h_bml, gbnd)
  mineval = gbnd(1)
  maxeval = gbnd(2)
  deallocate(gbnd)
  write(*,*) "Gershgorin: mineval = ", mineval, " maxeval = ", maxeval
  call prg_homolumogap(vv, icount, pp, mineval, maxeval, ehomo, elumo, egap)
  write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
       " egap = ", egap

  !! Start graph partitioning SP2

  !! Threshold the graph
  call bml_threshold(g_bml, gthreshold)

  !! Create graph partitioning of equal parts
  call prg_equalPartition(gp, nodesPerPart, mham%norbs)
  write(*,*)gp%reorder

  gp%mineval = mineval
  gp%maxeval = maxeval

  !! Calculate SP2 sequence
  call prg_sp2sequence(gp%pp, gp%maxIter, mineval, maxeval, ehomo, elumo, &
       errlimit)
  write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter


  !! Run subgraph SP2 loop
  call prg_timer_start(graphsp2_timer)
  call prg_subgraphSP2Loop(h_bml, g_bml, rho_bml, gp, sthreshold)
  call prg_timer_stop(graphsp2_timer)

  !! Calculate fnorm across subgraphs per iteration
  call prg_fnormGraph(gp)

  !! Calculate homo-lumo gap
  call prg_homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, &
       egap)
  write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
       " egap = ", egap

  call bml_add(g_bml,rho_bml,1.0d0,-1.0d0,threshold)
  write(*,*)"|rhoGP-rho|",bml_fnorm(g_bml)

  call bml_multiply(rho_bml, rho_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
  call bml_print_matrix("rhos_bml^2",aux_bml,0,10,0,10)
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)
  write(*,*)"|aux|",bml_fnorm(aux_bml)

  call bml_deallocate(g_bml)

  call prg_destroyGraphPartitioning(gp)
  
  call prg_progress_shutdown

end program main_sp2graph
