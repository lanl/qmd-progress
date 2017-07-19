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


  public::prg_kick_density, prg_kick_density_bml, prg_commutate,&
  prg_commutate_bml, prg_get_sparsity_cplxmat,prg_getcharge,&
  prg_get_sparsity_realmat, prg_lvni, prg_lvni_bml, prg_excitation,&
  prg_allocate_bml, prg_deallocate_bml

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
  !! \param N the number of orbitals in the density matrix
  !! \param S the overlap matrix
  !! \param SINV the inverse of the overlap matrix
  !! \param which_atom vector containing atom identification
  !! \param r direction vector for kick based on atom and kick_direc
  !! \param bmltype type of bml matrix desired for faster computation
  !! \param thresh threshold for bml matrix conversion
  !!
  subroutine prg_kick_density(kick_direc,kick_mag,dens,N,S,SINV,&
             which_atom,r,bmltype,thresh)
    implicit none
    integer                                     :: i
    integer, intent(in)                         :: kick_direc, N
    integer, allocatable, intent(in)            :: which_atom(:)
    real(dp), allocatable                       :: r(:,:)
    real(dp)                                    :: kick_mag, thresh
    complex(dp), allocatable, intent(inout)     :: dens(:,:)
    complex(dp), allocatable                    :: tmat1(:), tmat2(:), S(:,:), SINV(:,:)
    complex(dp)                                 :: telem
    type(bml_matrix_t)                          :: T1, T2, rho_bml, s_bml, sinv_bml
    character(len=*), intent(in)                :: bmltype

    allocate (tmat1(N))
    allocate (tmat2(N))

    do i=1,N
       telem = cmplx(0.0_dp,-kick_mag*r(which_atom(i),kick_direc))
       tmat1(i) = exp(telem)
       tmat2(i) = exp(-telem)
    enddo

    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,N,N, T1)
    call bml_add_diagonal(tmat1,T1)
    !call bml_convert_from_dense(bmltype,tmat1,T1,thresh,N)
    deallocate (tmat1)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,N,N, rho_bml)
    call bml_convert_from_dense(bmltype,dens,rho_bml,thresh,N)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,N,N, T2)
    call bml_multiply(T1,rho_bml,T2)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,N,N, s_bml)
    call bml_convert_from_dense(bmltype,S,s_bml,thresh,N)
    call bml_multiply(T2,s_bml,rho_bml)
    call bml_deallocate(s_bml)
    
    call bml_add_diagonal(tmat2,T1)
    !call bml_convert_from_dense(bmltype,tmat2,T1,thresh,N)
    deallocate (tmat2)
    call bml_multiply(rho_bml,T1,T2)
    call bml_deallocate(T1)
    call bml_zero_matrix(bmltype, BML_ELEMENT_COMPLEX,dp,N,N, sinv_bml)
    call bml_convert_from_dense(bmltype,SINV,sinv_bml,thresh,N)
    call bml_multiply(T2,sinv_bml,rho_bml)
    call bml_deallocate(sinv_bml)
    call bml_deallocate(T2)

    call bml_convert_to_dense(rho_bml,dens)

  end subroutine prg_kick_density


  !> This computes the commutator of two matricies.
  !! This routine does:
  !! \f$C = [A,B] = AB - BA\f$
  !! \param A matrix in first index of the commutator
  !! \param B matrix in second index of the commutator
  !! \param C output matrix that is the commutator of A and B
  !!
  subroutine prg_commutate(A,B,C)
    implicit none
    integer                               :: N
    real(dp)                              :: thresh
    complex(dp), allocatable, intent(in)  :: A(:,:), B(:,:)
    complex(dp), allocatable, intent(out) :: C(:,:)
    type(bml_matrix_t)                    :: A1, B1, C1
    thresh=1.0d-5
    N=SIZE(A,1)

    !Initialize bml matricies
    call bml_zero_matrix(BML_MATRIX_DENSE, BML_ELEMENT_COMPLEX,dp,N,N, A1)
    call bml_zero_matrix(BML_MATRIX_DENSE, BML_ELEMENT_COMPLEX,dp,N,N, B1)
    call bml_zero_matrix(BML_MATRIX_DENSE, BML_ELEMENT_COMPLEX,dp,N,N, C1)

    !Transform matricies into bml and do computation
    call bml_convert_from_dense(BML_MATRIX_DENSE,A,A1,thresh,N)
    call bml_convert_from_dense(BML_MATRIX_DENSE,B,B1,thresh,N)

    call bml_multiply(A1,B1,C1)
    call bml_multiply(B1,A1,C1,-1.0d0,1.0d0)
    call bml_convert_to_dense(C1,C)

    call bml_deallocate(A1)
    call bml_deallocate(B1)
    call bml_deallocate(C1)

  end subroutine prg_commutate


  !> This computes the sparsity of a complex matrix given a threshold value
  !! This routine does:
  !! \f$ f = \frac{N_0}{N_{tot}}\f$ where \f$f\f$ is the sparsity, \f$N_0\f$ is
  !! the number of values less than the threshold, and \f$N_{tot}\f$ is the total
  !! number of values. The sparsity and threshold are printed to the screen.
  !! \param matrix_type the bml matrix type
  !! \param element_type the bml element type
  !! \param sparsity_threshold the threshold for sparsity evaluation
  !! \param a_dense the dense complex matrix to be evaluated for sparsity
  !!
  subroutine prg_get_sparsity_cplxmat(matrix_type,element_type,sparsity_threshold,a_dense)
    implicit none
    character(len=*), intent(in)         :: matrix_type ,element_type
    complex(dp),intent(in)               :: a_dense(:,:)
    type(bml_matrix_t)                   :: a
    real(dp),intent(in)                  :: sparsity_threshold
    real(dp)                             :: convert_threshold
    integer                              :: asize
    character(len=20)                    :: test_type

    asize=SIZE(a_dense,1)
    convert_threshold=1.0d-5
    call bml_zero_matrix(matrix_type, element_type,dp,asize,asize, a)
    call bml_convert_from_dense(matrix_type, a_dense, a, convert_threshold,asize)
    !write(*,*)"thr,sparsity ",sparsity_threshold, bml_get_sparsity(a, sparsity_threshold)

    call bml_deallocate(a)

  end subroutine prg_get_sparsity_cplxmat

  !> This computes the sparsity of a real matrix given a threshold value
  !! This routine does:
  !! \f$ f = \frac{N_0}{N_{tot}}\f$ where \f$f\f$ is the sparsity, \f$N_0\f$ is
  !! the number of values less than the threshold, and \f$N_{tot}\f$ is the total
  !! number of values. The sparsity and threshold are printed to the screen.
  !! \param matrix_type the bml matrix type
  !! \param element_type the bml element type
  !! \param sparsity_threshold the threshold for sparsity evaluation
  !! \param a_dense the dense real matrix to be evaluated for sparsity
  !!
  subroutine prg_get_sparsity_realmat(matrix_type,element_type,sparsity_threshold,a_dense)
    implicit none
    character(len=*), intent(in)         :: matrix_type ,element_type
    real(dp),intent(in)                  :: a_dense(:,:)
    type(bml_matrix_t)                   :: a
    real(dp),intent(in)                  :: sparsity_threshold
    real(dp)                             :: convert_threshold
    integer                              :: asize

    asize=SIZE(a_dense,1)
    convert_threshold=1.0d-5
    call bml_zero_matrix(matrix_type, element_type,dp,asize,asize, a)
    call bml_convert_from_dense(matrix_type,a_dense,a,convert_threshold,asize)
    !write(*,*)"thr,sparsity ", sparsity_threshold,bml_get_sparsity(a,sparsity_threshold)

    call bml_deallocate(a)

  end subroutine prg_get_sparsity_realmat

  !> Performs Liouville-von Neumann integration using leap-frog method.
  !! This routine does:
  !! \f$\hat{\rho}(t+\Delta t)=\hat{\rho}(t-\Delta t) +2\Delta t\frac{\partial
  !! \hat{\rho}(t)}{\partial t}\f$ where the time derivative of the density matrix
  !! is defined as follows:
  !! \f$\frac{\partial\hat{\rho}(t)}{\partial
  !! t}=\frac{-i}{\hbar}\left(S^{-1}\hat{H}(t)\hat{\rho}(t)-\hat{\rho}(t)\hat{H}(t)S^{-1}\right)\f$
  !! \param H the Hamiltonian matrix at time t
  !! \param dt the timestep for integration
  !! \param hbar the Dirac constant (generally taken to be 1 in simulation units)
  !! \param rho_old the density matrix at previous time-step
  !! \param rho_cur the density matrix at current time-step
  !! \param rho_new the density matrix at next time-step
  !!
  subroutine prg_lvni(H, dt, hbar, rho_old, rho_cur, rho_new)
    implicit none
    real(dp)                                    :: thresh, hbar, dt
    complex(dp), allocatable, intent(in)        :: H(:,:), rho_old(:,:), rho_cur(:,:)
    complex(dp), allocatable, intent(out)       :: rho_new(:,:)
    complex(dp), allocatable                    :: partial_rho(:,:), com(:,:)
    complex(dp)                                 :: dR
    call prg_commutate(H,rho_cur,com)

    dR = cmplx(0.0_dp,2/hbar)
    partial_rho = -dR*dt*com
    rho_new = rho_old + partial_rho

  end subroutine prg_lvni

  !> Produce an excitation in the initially calculated density matrix to simulate
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
  !! \param which_atom vector containing atom identification
  !! \param r direction vector for kick based on atom and kick_direc
  !! \param matrix_type the type of BML format
  !! \param thresh the threshold for the BML matrix
  !!
  subroutine prg_kick_density_bml(kick_direc,kick_mag,rho_bml,s_bml,sinv_bml,which_atom,&
             r,matrix_type,thresh)
    implicit none
    integer                             :: i, N
    integer, intent(in)                 :: kick_direc
    integer, allocatable, intent(in)    :: which_atom(:)
    real(dp), allocatable               :: r(:,:)
    real(dp)                            :: kick_mag, thresh
    complex(dp), allocatable            :: tmat1(:), tmat2(:)
    complex(dp)                         :: telem
    type(bml_matrix_t)                  :: T1, T2, rho_bml,s_bml,sinv_bml
    character(len=*), intent(in)        :: matrix_type

    N=bml_get_N(rho_bml)
    allocate (tmat1(N))
    allocate (tmat2(N))

    do i=1,N
       telem = cmplx(0.0_dp,-kick_mag*r(which_atom(i),kick_direc))
       tmat1(i) = exp(telem)
       tmat2(i) = exp(-telem)
    enddo

    call bml_zero_matrix(matrix_type, BML_ELEMENT_COMPLEX,dp,N,N, T1)
    call bml_zero_matrix(matrix_type, BML_ELEMENT_COMPLEX,dp,N,N, T2)
    call bml_add_diagonal(tmat1,T1)
    !call bml_convert_from_dense(matrix_type,tmat1,T1,thresh,N)
    deallocate (tmat1)
    
    call bml_multiply(T1,rho_bml,T2)
    call bml_multiply(T2,s_bml,rho_bml)
    call bml_add_diagonal(tmat2,T1)
    !call bml_convert_from_dense(matrix_type,tmat2,T1,thresh,N)
    deallocate (tmat2)
    call bml_multiply(rho_bml,T1,T2)
    call bml_deallocate(T1)
    
    call bml_multiply(T2,sinv_bml,rho_bml)

    !write(*,*)"#######Density matrix Sparsity######"
    !write(*,*)"thr,sparsity ",thresh,bml_get_sparsity(RHO,thresh)

    call bml_deallocate(T2)

  end subroutine prg_kick_density_bml


  !> This computes the commutator of two matricies.
  !! This routine does:
  !! \f$C = [A,B] = AB - BA\f$
  !! \param A matrix in first index of the commutator
  !! \param B matrix in second index of the commutator
  !! \param C output matrix that is the commutator of A and B
  subroutine prg_commutate_bml(A,B,C)
    implicit none
    type(bml_matrix_t)  :: A, B, C

    call bml_multiply(A,B,C)
    call bml_multiply(B,A,C,-1.0d0,1.0d0)

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
  !! \param rho_cur the density matrix at current time-step
  !! \param rho_new the density matrix at next time-step
  !!
  subroutine prg_lvni_bml(H, sinv_bml, dt, hbar, rho_old, rho_cur,aux_bml,matrix_type)
    implicit none
    real(dp)                            :: hbar, dt
    type(bml_matrix_t)                  :: H, sinv_bml
    type(bml_matrix_t)                  :: rho_cur, rho_old, aux_bml
    complex(dp)                         :: dR
    character(len=*), intent(in)        :: matrix_type

    dR = 2.0_dp*dt*cmplx(0.0_dp,-1.0_dp/hbar)
    call prg_commutate_bml(H,rho_cur,aux_bml)
    call bml_scale_cmplx(dR,aux_bml)
    call bml_add(aux_bml,rho_old,1.0d0,1.0d0)
    call bml_copy(rho_cur,rho_old)
    call bml_copy(aux_bml,rho_cur)

  end subroutine prg_lvni_bml

  !> Constructs the charges from the density matrix.
  !! \param rho_bml Density matrix in bml format.
  !! \param over_bml Overlap matrix in bml format.
  !! \param charges the array of charges.
  !! \param aux_bml the auxiliary matrix in bml format.
  !! \param diag array of diagonal elements.
  !! \param spindex Start and end index for every atom in the system.
  !! \param z
  !! \param nats the number of atoms
  !! \param N
  !! \param matrix_type type of bml matrix.
  !!
  subroutine prg_getcharge(rho_bml,over_bml,charges,aux_bml,diag,z,spindex,N,nats,matrix_type)
    implicit none
    integer                               ::  i, j, k, nats
    integer, allocatable,intent(in)       ::  N(:), spindex(:)
    complex(dp),allocatable               ::  diag(:)
    real(dp), allocatable                 ::  charges(:)
    real(dp), intent(in)                  ::  z(:)
    type(bml_matrix_t), intent(in)        ::  over_bml, rho_bml
    type(bml_matrix_t)                    ::  aux_bml
    character(len=*), intent(in)          ::  matrix_type

    if(.not.allocated(charges)) allocate(charges(nats))
    call bml_multiply(rho_bml,over_bml,aux_bml,2.0_dp)
    call bml_get_diagonal(aux_bml,diag)
    k=0
    do i = 1,nats
       charges(i)=0.0_dp
       do j = 1, N(spindex(i))
          k = k+1
          charges(i) = charges(i) + diag(k)
       enddo
       charges(i) = z(spindex(i)) - charges(i)
    enddo

  end subroutine prg_getcharge

  !> Allocates a number of matricies in BML format required for quantum dynamics
  !! calculations.
  !! \param rho Density matrix in standard format.
  !! \param rho_bml a bml copy of the aforementioned density matrix
  !! \param rhoold Density matrix at past timestep in std format.
  !! \param rhoold_bml a bml copy of the aforementioned rhoold matrix
  !! \param h1 Hamiltonian matrix in std format.
  !! \param h1_bml a bml copy of the aforementioned Hamiltonian matrix
  !! \param aux the auxiliary matrix in bml format
  !! \param diag array of diagonal elements.
  !! \param N the size of the square matricies listed above
  !! \param bmltype the type of bml matricies desired
  !! \param thresh the threshold for the bml matricies
  !!
  subroutine prg_allocate_bml(rho,rho_bml,rhoold,rhoold_bml,h1,&
             h1_bml,aux,diag,N,bmltype,thresh)
    implicit none
    type(bml_matrix_t),intent(inout) :: rho_bml, rhoold_bml, h1_bml,aux
    complex(dp),intent(in)           :: rho(:,:), rhoold(:,:), h1(:,:)
    complex(dp),allocatable          :: diag(:)
    integer                          :: N
    real(dp),intent(in)              :: thresh
    character(len=*), intent(in)     :: bmltype

    allocate(diag(N))
    call bml_zero_matrix(bmltype,BML_ELEMENT_COMPLEX,dp,N,N,rho_bml)
    call bml_convert_from_dense(bmltype,rho,rho_bml,thresh,N)
    call bml_zero_matrix(bmltype,BML_ELEMENT_COMPLEX,dp,N,N,rhoold_bml)
    call bml_convert_from_dense(bmltype,rhoold,rhoold_bml,thresh,N)
    call bml_zero_matrix(bmltype,BML_ELEMENT_COMPLEX,dp,N,N,h1_bml)
    call bml_convert_from_dense(bmltype,h1,h1_bml,thresh,N)
    call bml_zero_matrix(bmltype,BML_ELEMENT_COMPLEX,dp,N,N,aux)

  end subroutine prg_allocate_bml

  !> Deallocates aux matrix and diagonal array used in several quantum dynamics
  !! calculations.
  !! \param aux bml matrix used for storage of values in several routines
  !! \param d complex array storing the diagonal of the charge matrix
  !!
  subroutine prg_deallocate_bml(aux,d)
    implicit none
    type(bml_matrix_t) :: aux
    complex(dp),allocatable :: d(:)

    deallocate(d)
    call bml_deallocate(aux)

  end subroutine prg_deallocate_bml

end module prg_quantumdynamics_mod
