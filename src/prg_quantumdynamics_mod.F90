!> A module to add in common quantum dynamical operations.
!! \brief This module contains routines that perform the following tasks: apply
!! a apply an excitation or  perturbation to the initial density matrix, compute
!! the comutator of two two matricies, calculate the sparsity of a real or
!! complex matrix, and time evolve a density matrix using Liouville-von Neumann
!! equation with the leap-frog method of integration.
!! @ingroup PROGRESS
!!
module prg_quantumdynamics_mod
  use bml
  implicit none

  private
  integer, parameter :: dp = kind(1.0d0)


  public::prg_kick_density, prg_kick_density_bml,&
  prg_commutate_bml, prg_get_sparsity_cplxmat, prg_getcharge,&
  prg_getdipole, prg_get_sparsity_realmat, prg_lvni_bml,&
  prg_excitation 

contains

  !> Provides perturbation to initial density matrix in the form of an electric
  !! field kick.
  !! This routine does:
  !! \f$\hat{\rho_{kick}} =
  !! \exp{\frac{-i}{\hbar}\hat{V}}\hat{\rho}\hat{S}\exp{\frac{i}{\hbar}\hat{V}}\hat{S^{-1}}\f$
  !! where \f$\hat{V}\f$ is the field disturbance.
  !! \param kick_direc the direction of the kick in the electric field
  !! \param kick_mag the magnitude of the kick in the electric field
  !! \param dens the initial density matrix to be kicked.
  !! \param norbs the number of orbitals in the density matrix
  !! \param S the overlap matrix
  !! \param SINV the inverse of the overlap matrix
  !! \param which_atom vector containing atom identification
  !! \param r direction vector for kick based on atom and kick_direc
  !! \param bmltype type of BML matrix desired for faster computation
  !! \param thresh threshold for BML matrix conversion
  !!
  subroutine prg_kick_density(kick_direc,kick_mag,dens,norbs,M,S,SINV,&
             which_atom,r,bmltype,thresh)
    implicit none
    integer                                     :: i
    integer, intent(in)                         :: kick_direc, norbs, M
    integer, allocatable, intent(in)            :: which_atom(:)
    real(dp), allocatable                       :: r(:,:)
    real(dp)                                    :: kick_mag, thresh
    complex(dp), allocatable, intent(inout)     :: dens(:,:)
    complex(dp), allocatable                    :: tmat1(:), tmat2(:), S(:,:), SINV(:,:)
    complex(dp)                                 :: telem
    type(bml_matrix_t)                          :: T1, T2, rho_bml, s_bml, sinv_bml
    character(len=*), intent(in)                :: bmltype

    allocate (tmat1(norbs))
    allocate (tmat2(norbs))

    do i=1,norbs
       telem = cmplx(0.0_dp,-kick_mag*r(which_atom(i),kick_direc))
       tmat1(i) = exp(telem)
       tmat2(i) = exp(-telem)
    enddo

    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,norbs,M, T1)
    call bml_set_diagonal(T1,tmat1,thresh)
    deallocate (tmat1)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,norbs,M, rho_bml)
    call bml_convert_from_dense(bmltype,dens,rho_bml,thresh,M)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,norbs,M, T2)
    call bml_multiply(T1,rho_bml,T2)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,norbs,M, s_bml)
    call bml_convert_from_dense(bmltype,S,s_bml,thresh,M)
    call bml_multiply(T2,s_bml,rho_bml)
    call bml_deallocate(s_bml)

    call bml_set_diagonal(T1,tmat2,thresh)
    deallocate (tmat2)
    call bml_multiply(rho_bml,T1,T2)
    call bml_deallocate(T1)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,norbs,M, sinv_bml)
    call bml_convert_from_dense(bmltype,SINV,sinv_bml,thresh,M)
    call bml_multiply(T2,sinv_bml,rho_bml)
    call bml_deallocate(sinv_bml)
    call bml_deallocate(T2)

    write(*,*)"#######Density matrix Sparsity######"
    write(*,*)"thr,sparsity ",1.0d-3,bml_get_sparsity(rho_bml,1.0d-3)
    write(*,*)"thr,sparsity ",1.0d-4,bml_get_sparsity(rho_bml,1.0d-4)
    write(*,*)"thr,sparsity ",1.0d-5,bml_get_sparsity(rho_bml,1.0d-5)
    write(*,*)"thr,sparsity ",1.0d-6,bml_get_sparsity(rho_bml,1.0d-6)

   call bml_convert_to_dense(rho_bml,dens)

  end subroutine prg_kick_density


  !> This computes the sparsity of a complex matrix given a threshold value
  !! This routine does:
  !! \f$ f = \frac{N_0}{N_{tot}}\f$ where \f$f\f$ is the sparsity, \f$N_0\f$ is
  !! the number of values less than the threshold, and \f$N_{tot}\f$ is the total
  !! number of values. The sparsity and threshold are printed to the screen.
  !! \param matrix_type the BML matrix type
  !! \param element_type the BML element type
  !! \param thresh the threshold for sparsity evaluation
  !! \param a_dense the dense complex matrix to be evaluated for sparsity
  !!
  subroutine prg_get_sparsity_cplxmat(matrix_type,element_type,thresh,a_dense)
    implicit none
    character(len=*), intent(in)         :: matrix_type ,element_type
    complex(dp),intent(in)               :: a_dense(:,:)
    type(bml_matrix_t)                   :: a
    real(dp),intent(in)                  :: thresh
    integer                              :: asize
    character(len=20)                    :: test_type

    asize=SIZE(a_dense,1)
    call bml_zero_matrix(matrix_type, element_type,dp,asize,asize, a)
    call bml_convert_from_dense(matrix_type, a_dense, a, thresh,asize)
    !write(*,*)"thr,sparsity ",thresh, bml_get_sparsity(a, thresh)

    call bml_deallocate(a)

  end subroutine prg_get_sparsity_cplxmat

  !> This computes the sparsity of a real matrix given a threshold value
  !! This routine does:
  !! \f$ f = \frac{N_0}{N_{tot}}\f$ where \f$f\f$ is the sparsity, \f$N_0\f$ is
  !! the number of values less than the threshold, and \f$N_{tot}\f$ is the total
  !! number of values. The sparsity and threshold are printed to the screen.
  !! \param matrix_type the BML matrix type
  !! \param element_type the BML element type
  !! \param thresh the threshold for sparsity evaluation
  !! \param a_dense the dense real matrix to be evaluated for sparsity
  !!
  subroutine prg_get_sparsity_realmat(matrix_type,element_type,thresh,a_dense)
    implicit none
    character(len=*), intent(in)         :: matrix_type ,element_type
    real(dp),intent(in)                  :: a_dense(:,:)
    type(bml_matrix_t)                   :: a
    real(dp),intent(in)                  :: thresh
    integer                              :: asize

    asize=SIZE(a_dense,1)
    call bml_zero_matrix(matrix_type, element_type,dp,asize,asize, a)
    call bml_convert_from_dense(matrix_type,a_dense,a,thresh,asize)
    !write(*,*)"thr,sparsity ", thresh,bml_get_sparsity(a,thresh)

    call bml_deallocate(a)

  end subroutine prg_get_sparsity_realmat

  
  !> Provides perturbation to initial density matrix in the form of an electric
  !! field kick given input matricies in BML format.
  !! This routine does:
  !! \f$\hat{\rho_{kick}} =
  !! \exp{\frac{-i}{\hbar}\hat{V}}\hat{\rho}\hat{S}\exp{\frac{i}{\hbar}\hat{V}}\hat{S^{-1}}\f$
  !! where \f$\hat{V}\f$ is the field disturbance.
  !! \param kick_direc the direction of the kick in the electric field
  !! \param kick_mag the magnitude of the kick in the electric field
  !! \param rho_bml the initial density matrix to be kicked in BML format.
  !! \param s_bml the overlap matrix
  !! \param sinv_bml the inverse of the overlap matrix
  !! \param M maximum number of nonzero values per row in BML matrix
  !! \param which_atom vector containing atom identification
  !! \param r position vector for kicked atom 
  !! \param matrix_type the type of BML format
  !! \param thresh the threshold for the BML matrix
  !!
  subroutine prg_kick_density_bml(kick_direc,kick_mag,rho_bml,s_bml,sinv_bml,M,which_atom,&
             r,matrix_type,thresh)
    implicit none
    integer                             :: i, norbs, M
    integer, intent(in)                 :: kick_direc
    integer, allocatable, intent(in)    :: which_atom(:)
    real(dp), allocatable               :: r(:,:)
    real(dp)                            :: kick_mag, thresh
    complex(dp), allocatable            :: tmat1(:), tmat2(:)
    complex(dp)                         :: telem
    type(bml_matrix_t)                  :: T1, T2, rho_bml,s_bml,sinv_bml
    character(len=*), intent(in)        :: matrix_type

    norbs=bml_get_N(rho_bml)
    allocate (tmat1(norbs))
    allocate (tmat2(norbs))

    do i=1,norbs
       telem = cmplx(0.0_dp,-kick_mag*r(which_atom(i),kick_direc))
       tmat1(i) = exp(telem)
       tmat2(i) = exp(-telem)
    enddo

    call bml_zero_matrix(matrix_type, BML_ELEMENT_COMPLEX,dp,norbs,M, T1)
    call bml_zero_matrix(matrix_type, BML_ELEMENT_COMPLEX,dp,norbs,M, T2)
    call bml_set_diagonal(T1,tmat1,thresh)
    deallocate (tmat1)

    call bml_multiply(T1,rho_bml,T2)
    call bml_multiply(T2,s_bml,rho_bml)
    call bml_set_diagonal(T1,tmat2,thresh)
    deallocate (tmat2)
    call bml_multiply(rho_bml,T1,T2)
    call bml_deallocate(T1)

    call bml_multiply(T2,sinv_bml,rho_bml)
    call bml_deallocate(T2)

  end subroutine prg_kick_density_bml


  !> This computes the commutator of two matricies.
  !! This routine does:
  !! \f$C = [A,B] = AB - BA\f$
  !! \param A matrix in first index of the commutator
  !! \param B matrix in second index of the commutator
  !! \param C output matrix that is the commutator of A and B
  !! \param thresh the BML threshold
  !!
  subroutine prg_commutate_bml(A,B,C,thresh)
    implicit none
    type(bml_matrix_t)  :: A, B, C
    real(dp), intent(in) :: thresh

    call bml_multiply(A,B,C,1.0_dp,1.0_dp,thresh)
    call bml_multiply(B,A,C,-1.0_dp,1.0_dp,thresh)

  end subroutine prg_commutate_bml


  !> Performs Liouville-von Neumann integration using leap-frog method.
  !! This routine does:
  !! \f$\hat{\rho}(t+\Delta t)=\hat{\rho}(t-\Delta t) +2\Delta t\frac{\partial
  !! \hat{\rho}(t)}{\partial t}\f$ where the time derivative of the density matrix
  !! is defined as follows:
  !! \f$\frac{\partial\hat{\rho}(t)}{\partial
  !! t}=\frac{-i}{\hbar}\left(S^{-1}\hat{H}(t)\hat{\rho}(t)-\hat{\rho}(t)\hat{H}(t)S^{-1}\right)\f$
  !! \param H the Hamiltonian matrix at time t
  !! \param sinv_bml the inverse overlap matrix
  !! \param dt the timestep for integration
  !! \param hbar the Dirac constant (generally taken to be 1 in simulation units)
  !! \param rho_old the density matrix at previous time-step
  !! \param rho_bml the density matrix at current time-step
  !! \param aux_bml the temp matrix used for value storage during computation
  !! \param matrix_type the type of BML matrix
  !! \param thresh the threshold for the BML matrix
  !!
  subroutine prg_lvni_bml(H, sinv_bml, dt, hbar, rho_old, rho_bml,aux_bml,matrix_type,M,thresh)
    implicit none
    integer                             :: norbs, M
    real(dp)                            :: hbar, dt
    type(bml_matrix_t)                  :: H, sinv_bml
    type(bml_matrix_t)                  :: rho_bml, rho_old, aux_bml,aux_1
    complex(dp)                         :: dR
    real(dp), intent(in)                :: thresh
    character(len=*), intent(in)        :: matrix_type
    
    norbs=bml_get_N(rho_bml)
    dR = 2.0_dp*dt*cmplx(0.0_dp,-1.0_dp/hbar)
    call bml_multiply(sinv_bml,H,aux_bml,1.0_dp,0.0_dp,thresh)
    call bml_zero_matrix(matrix_type,BML_ELEMENT_COMPLEX,dp,norbs,M,aux_1)
    call bml_multiply(aux_bml,rho_bml,aux_1,1.0_dp,0.0_dp,thresh)
    call bml_multiply(rho_bml,H,aux_bml,1.0_dp,0.0_dp,thresh)
    call bml_multiply(aux_bml,sinv_bml,aux_1,-1.0_dp,1.0d0)
    call bml_scale(dR,aux_1,aux_bml)
    call bml_deallocate(aux_1)
    call bml_add(aux_bml,rho_old,1.0_dp,1.0_dp,thresh)
    call bml_copy(rho_bml,rho_old)
    call bml_copy(aux_bml,rho_bml)
    
  end subroutine prg_lvni_bml


  !> Constructs the charges from the density matrix.
  !! \param rho_bml Density matrix in BML format.
  !! \param over_bml Overlap matrix in BML format.
  !! \param charges the array of charges.
  !! \param aux_bml the auxiliary matrix in BML format.
  !! \param spindex Start and end index for every atom in the system.
  !! \param z
  !! \param nats the number of atoms
  !! \param N
  !! \param thresh threshold for the BML matrix
  !!
  subroutine prg_getcharge(rho_bml,over_bml,charges,aux_bml,z,spindex,N,nats,thresh)
    implicit none
    integer                               ::  i, j, k, nats,norbs
    integer, allocatable,intent(in)       ::  N(:), spindex(:)
    complex(dp),allocatable               ::  auxd(:,:)
    real(dp), allocatable                 ::  charges(:)
    real(dp), intent(in)                  ::  thresh,z(:)
    type(bml_matrix_t), intent(in)        ::  over_bml, rho_bml
    type(bml_matrix_t)                    ::  aux_bml

    norbs=bml_get_N(rho_bml)
    allocate(auxd(norbs,norbs))
    auxd=cmplx(0.0_dp,0.0_dp)
    if(.not.allocated(charges)) allocate(charges(nats))
    call bml_multiply(rho_bml,over_bml,aux_bml,1.0_dp,0.0_dp,thresh)
    call bml_convert_to_dense(aux_bml,auxd)
    k=0
    do i = 1,nats
       charges(i)=0.0_dp
       do j = 1, N(spindex(i))
          k = k+1
          charges(i) = charges(i) + 2.0_dp*auxd(k,k)
       enddo
       charges(i) = -charges(i) + z(spindex(i)) 
    enddo
    deallocate(auxd)

  end subroutine prg_getcharge


  !> This routine computes the dipole moment of the system with units determined
  !! by the units of the coordinate matrix and charges given.
  !! \param charges Charge on each atom.
  !! \param r Coordinate matrix of the atoms.
  !! \param p Dipole moment vector.
  !!
  subroutine prg_getdipole(charges,r,p)
    implicit none
    integer                  ::  i, norbs
    real(dp), intent(in)     ::  charges(:), r(:,:)
    real(dp), intent(inout)  ::  p(3)

    norbs = SIZE(charges,1)
    p = 0.0_dp

    do i=1,norbs
       p=p+r(:,i)*charges(i)
    enddo

  end subroutine prg_getdipole


  !> Produce an excitation in the initially calculated density matrix to
  !simulate
  !! photo-excitation.
  !! This routine does:
  !! Manually induces electronic excitation based on the initial filling matrix
  !! by promoting an electron.
  !! \param fill_mat the initial filling matrix
  !! \param orbit_orig the origin orbital with an electron to be promoted
  !! \param orbit_exci the destination orbital of the promoted electron
  !!
  subroutine prg_excitation(fill_mat, orbit_orig, orbit_exci)
    implicit none
    integer, intent(inout)      :: fill_mat(:)
    integer, intent(in)         :: orbit_orig, orbit_exci

    fill_mat(orbit_orig) = fill_mat(orbit_orig) - 1.0d0
    fill_mat(orbit_exci) = fill_mat(orbit_exci) + 1.0d0

  end subroutine prg_excitation

end module prg_quantumdynamics_mod
