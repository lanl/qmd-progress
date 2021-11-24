!> High-level program to construct a model Hamiltonian
!!
program hmodel

  !BML lib.
  use bml

  !PROGRESS lib modes
  use prg_modelham_mod
  use prg_system_mod
  use prg_densitymatrix_mod
  use prg_dos_mod
  use prg_sp2_mod
  use prg_timer_mod
  use prg_extras_mod
  use prg_parallel_mod

  implicit none
  integer, parameter ::  dp = kind(1.0d0)
  integer ::  norbs,seed
  integer ::  verbose
  character(20) :: filename
  type(bml_matrix_t) ::  ham_bml,rho_bml,rhos_bml,evects_bml,aux_bml
  type(mham_type) ::  mham
  type(system_type) ::  sys
  real(dp) ::  threshold, bndfil
  real(dp), allocatable :: trace(:)
  real(dp), allocatable :: eigenvalues(:)
  real(dp) :: ef,sparsity,dec,mlsi,mlsf,bnorm
  character(20) :: bml_dmode

  call prg_initParallel()

  if (getNRanks().gt.1)then
    if (printRank() .eq. 1)print*,'BML_DMODE_DISTRIBUTED'
    bml_dmode = BML_DMODE_DISTRIBUTED
  else
    print*,'BML_DMODE_SEQUENTIAL'
    bml_dmode = BML_DMODE_SEQUENTIAL
  endif

  !Parsing input file.
  call getarg(1,filename)
  call prg_parse_mham(mham,trim(adjustl(filename))) !Reads the input for modelham

  !Number of orbitals/matrix size
  norbs=mham%norbs

  !Allocating bml matrices
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,ham_bml, &
       &               bml_dmode)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,rho_bml, &
       &               bml_dmode)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs, &
       &               evects_bml, bml_dmode)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,rhos_bml, &
       &               bml_dmode)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,aux_bml, &
       &               bml_dmode)

  seed = 1000 !Seed to reproduce the Hamiltonian build
  verbose = 1 !Verbosity level
  threshold = 1.0d-5 !Threshold value for the matrices through the whole code
  bndfil = 0.5d0 !Fraction of orbitals that will be filled

  allocate(eigenvalues(norbs))

  !Constructing the Hamiltonian
  call prg_twolevel_model(mham%ea, mham%eb, mham%dab, mham%daiaj, mham%dbibj, &
       &mham%dec, mham%rcoeff, mham%reshuffle, mham%seed, ham_bml, verbose)
  call bml_threshold(ham_bml, threshold)
  call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)

  sparsity = bml_get_sparsity(ham_bml, 1.0D-5)
  if (printRank() .eq. 1)write(*,*)"Sparsity Ham=",sparsity

  !Computing the density matrix with diagonalization
  if (printRank() .eq. 1)print*,'prg_build_density_T0'
  mlsi = mls()
  call prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil, eigenvalues)
  mlsf = mls()
  if (printRank() .eq. 1)write(*,*)"Time_for_prg_build_density_T0",mlsf-mlsi

  sparsity = bml_get_sparsity(rho_bml, 1.0D-5)
  if (printRank() .eq. 1)write(*,*)"Sparsity Rho=",sparsity

  !Getting the fermi level
  ef = (eigenvalues(int(norbs/2)+1) + eigenvalues(int(norbs/2)))/2
  eigenvalues = eigenvalues - ef

  !Writting the total DOS
  call prg_write_tdos(eigenvalues, 0.05d0, 10000, -20.0d0, 20.0d0, "tdos.dat")

  !Solving for Rho using SP2
  mlsi = mls()
  call prg_sp2_alg1(ham_bml,rhos_bml,threshold,bndfil,15,100 &
       ,"Rel",1.0D-10,20)
  mlsf = mls()
  if (printRank() .eq. 1)write(*,*)"Time_for_prg_sp2_alg1",mlsf-mlsi
  call bml_print_matrix("rho_bml",rho_bml,0,10,0,10)
  call bml_print_matrix("rhos_bml",rhos_bml,0,10,0,10)

  call bml_copy(rhos_bml,aux_bml)
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)
  bnorm=bml_fnorm(aux_bml)
  if (printRank() .eq. 1)write(*,*)"|DM_sp2-DM_diag|",bnorm

  call bml_multiply(rhos_bml, rhos_bml, aux_bml, 0.5_dp, 0.0_dp, threshold)
  call bml_print_matrix("rhos_bml^2",aux_bml,0,10,0,10)
  call bml_add(aux_bml,rhos_bml,1.0d0,-1.0d0,threshold)
  bnorm=bml_fnorm(aux_bml)
  if (printRank() .eq. 1)write(*,*)"|DM_sp2-DM_sp2^2|",bnorm

  call bml_multiply(ham_bml,rhos_bml,aux_bml,1.0_dp,0.0_dp,threshold)
  call bml_multiply(rhos_bml,ham_bml,aux_bml,1.0_dp,-1.0_dp,threshold)
  bnorm=bml_fnorm(aux_bml)
  if (printRank() .eq. 1)write(*,*)"|DM_sp2*H-H*DM_sp2|",bnorm

  !call bml_write_matrix(ham_bml, "hamiltonian.mtx")

  call prg_shutdownParallel()

end program hmodel
