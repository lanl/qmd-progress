!> Variables used in the code
!!
module gpmdcov_vars

  use bml
 
 !Progress and LATTE lib modules.
  use prg_progress_mod
  use prg_system_mod
  use prg_ptable_mod
  use latteparser_latte_mod
  use huckel_latte_mod
  use tbparams_latte_mod
  use ham_latte_mod
  use coulomb_latte_mod
  use prg_charges_mod
  use prg_initmatrices_mod
  use prg_genz_mod
  use prg_nonortho_mod
  use prg_pulaymixer_mod
  use prg_dos_mod
  use prg_densitymatrix_mod
  use gpmdcov_neighbor_mod
  use gpmdcov_dispersion_mod
  !use neighborlist_latte_mod
  use ppot_latte_mod
  use hsderivative_latte_mod
  use slaterkosterforce_latte_mod
  use prg_pulaycomponent_mod
  use nonorthocoulombforces_latte_mod

 ! Graph partitioning modules
  use prg_parallel_mod
  use prg_timer_mod
  use prg_graphsp2parser_mod
  use prg_sp2parser_mod
  use prg_sp2_mod
  use prg_graph_mod
  use prg_subgraphLoop_mod
  use prg_homolumo_mod
  use prg_xlbo_mod
  use md_latte_mod
  use prg_partition_mod
  use prg_extras_mod

 ! gpmdcov modules
  use gpmdcov_parser_mod
 
  implicit none 
 
  integer, parameter                ::  dp = kind(1.0d0)
  integer                           ::  seed = 1
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable         ::  intKind(:)
  character(50)                     ::  inputfile, filename
  character(2)                      ::  auxchar
  integer                           ::  mdstep, Nr_SCF_It, i, icount, ierr, num
  integer                           ::  j, nel, norb, pp(100), nnodes, iii, kk
  integer                           ::  nparts, niter=500, npat, ipt, iptt
  integer                           ::  ii, jj, iscf, norb_core, totalNorbs
  integer                           ::  mdim, shift, nranks, norbsInRank
  integer, allocatable              ::  hindex(:,:), hnode(:), vectorint(:), norbsInEachCHAtRank(:), norbsInEachRank(:)
  integer, allocatable              ::  xadj(:), adjncy(:), CH_count(:), norbsInEachCH(:)
  integer, allocatable              ::  part(:), core_count(:), Halo_count(:,:)
  integer, allocatable              ::  partsInEachRank(:), reshuffle(:,:), npartsVect(:), displ(:),  PartsInRankI(:)
  integer, allocatable              ::  whichParts_guess_saved(:)
  real(dp)                          ::  C0, C1, C2, C3
  real(dp)                          ::  C4, C5, ECoul, ECoulU, ECoulK, ECoulR, EKIN, beta, kbt
  real(dp)                          ::  EPOT, ERep, Energy, Etot, nocc
  real(dp)                          ::  F2V, KE2T, MVV2KE, EVOVERV2P, M_prg_init
  real(dp)                          ::  TRRHOH, Temp, Time, alpha
  real(dp)                          ::  bndfil, bndfilTotal, cc, coulcut, dt
  real(dp)                          ::  dx, egap, egap_glob, ehomo, elumo
  real(dp)                          ::  kappa, scferror, traceMult, vv(100)
  real(dp)                          ::  sumCubes, maxCH, Ef, smooth_maxCH, pnorm=6
  real(dp)                          ::  dvdw, d, mls_i, Efstep, costperrank, costperrankmax, costperrankmin
  real(dp)                          ::  sparsity, entropy
  real(dp)                          ::  HOMO, LUMO
  real(dp)                          ::  savets
  real(dp), allocatable             ::  eigenvalues(:), evals(:), fvals(:), dvals(:)
  real(dp), allocatable             ::  evalsAll(:), fvalsAll(:), dvalsAll(:)
  real(dp), allocatable             ::  evalsInRank(:), fvalsInRank(:), dvalsInRank(:)
  real(dp), allocatable             ::  GFPUL(:,:), GFSCOUL(:,:), PairForces(:,:), DispForces(:,:)
  real(dp), allocatable             ::  SKForce(:,:), VX(:), VY(:), VZ(:), collectedforce(:,:), smdForce(:,:)
  real(dp), allocatable             ::  charges_old(:), coul_forces(:,:), coul_forces_k(:,:), coul_forces_r(:,:)
  real(dp), allocatable             ::  coul_pot_k(:), coul_pot_r(:), dqin(:,:), dqout(:,:)
  real(dp), allocatable             ::  eigenvals(:), gbnd(:), n(:), n_0(:)
  real(dp), allocatable             ::  n_1(:), n_2(:), n_3(:), n_4(:), acceprat(:)
  real(dp), allocatable             ::  n_5(:), onsitesH(:,:), onsitesS(:,:), rhoat(:)
  real(dp), allocatable             ::  origin(:), row(:), row1(:), auxcharge(:), auxcharge1(:)
  real(dp), allocatable             ::  g_dense(:,:),tch, Ker(:,:)
  real(dp), allocatable             ::  voltagev(:)
  integer, allocatable             ::  freeze_list(:)
  type(bml_matrix_t)                ::  aux_bml, dH0x_bml, dH0y_bml, dH0z_bml
  type(bml_matrix_t)                ::  dSx_bml, dSy_bml, dSz_bml, eigenvects
  type(bml_matrix_t)                ::  g_bml, ham0_bml, ham_bml
  type(bml_matrix_t)                ::  over_bml, rho_bml, rhoat_bml
  type(bml_matrix_t)                ::  rhoh_bml, zmat_bml, gch_bml
  type(bml_matrix_t)                ::  copy_g_bml, gcov_bml, aux1_bml
  type(graph_partitioning_t)        ::  gp
  type(graph_partitioning_t)        ::  gpat
  type(gsp2data_type)               ::  gsp2
  type(intpairs_type), allocatable  ::  intPairsH(:,:), intPairsS(:,:)
  type(latte_type)                  ::  lt
  type(neighlist_type)              ::  nl
  type(ppot_type), allocatable      ::  ppot(:,:)
  type(disppot_type), allocatable      ::  disppot(:,:)
  type(sp2data_type)                ::  sp2
  type(system_type)                 ::  sy
  type(system_type), allocatable    ::  syprt(:), oldsyprt(:)
  type(system_type), allocatable    ::  syprtk(:)
  type(tbparams_type)               ::  tb
  type(xlbo_type)                   ::  xl
  type(gpmd_type)                   ::  gpmdt
  type(estructout_type)             ::  estrout
  logical                           ::  first_part = .true.
  logical                           ::  converged = .false.
  logical                           ::  firstKernel = .true.
  logical                           ::  newPart = .false.
  logical                           ::  eig = .true.  
  logical                           ::  nlistSparse = .false.
  logical                           ::  reactionDone = .false.
  logical, save                     ::  lib_init = .false.
  logical, save                     ::  lib_mode = .false.
  logical, save                     ::  lib_init2 = .false.
  logical, save                     ::  err_status = .false.
  logical, save                     ::  vinit = .false.

  type(bml_matrix_t)                :: ZK1_bml, ZK2_bml, ZK3_bml
  type(bml_matrix_t)                :: ZK4_bml, ZK5_bml, ZK6_bml
  integer                           :: igenz !Couter to keep track of the times zmat is computed.
  logical                           :: tZSP, restart
  type(genZSPinp)                   :: zsp

  integer, allocatable              :: xadj_cov(:), adjncy_cov(:), CH_count_cov(:), partsSizes(:)
  integer, allocatable              :: part_cov(:), core_count_cov(:), Halo_count_cov(:,:)
  integer, allocatable              :: thisPartSize(:)
  integer                           :: vsize(2)
  integer                           :: nparts_cov, myRank, numRanks
  integer                           :: getKernel_cont=0,getKernel_byBlocks_cont=0
  integer                           :: applyKernel_cont=0,getKernel_byParts_cont=0 
  integer                           :: rankN_update_byParts_cont=0

end module gpmdcov_vars
