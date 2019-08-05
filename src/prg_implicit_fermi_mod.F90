! The Implicit Recursive Fermi O(N) module.
!! \ingroup PROGRESS
!! \brief This subroutine implements Niklasson's implicit recursive fermi dirac exact 
!! density matrix purification algorithm.
!!
module prg_implicit_fermi_mod

  use bml   
  use prg_normalize_mod
  use prg_timer_mod
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_implicit_fermi
  public :: prg_test_density_matrix

contains

  !> Recursive Implicit Fermi Dirac.
  !! \param h_bml Input Hamiltonian matrix.
  !! \param p_bml Output density matrix.
  !! \param nsteps Number of recursion steps.
  !! \param nocc Number of occupied states.
  !! \param mu Shifted chemical potential
  !! \param beta Input inverse temperature.
  !! \param osteps Outer loop steps to converge chemical potential
  !! \param occErrLimit Occupation error limit.
  !! \param threshold Threshold for multiplication.
  !!
  subroutine prg_implicit_fermi(h_bml, p_bml, nsteps, nocc, &
       mu, beta, osteps, occErrLimit, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: p_bml
    integer, intent(in) :: osteps, nsteps
    real(dp), intent(in) :: nocc, threshold
    real(dp), intent(in) :: occErrLimit, beta
    real(dp), intent(inout) :: mu

    type(bml_matrix_t) :: wtmp_bml, y_bml, d_bml, p2_bml, I_bml, ai_bml
    real(dp) :: trdPdmu, trP0, occErr
    real(dp) :: cnst, ofactor
    real(dp) :: cg_tol
    real(dp), allocatable :: trace(:), gbnd(:)
    character(20) :: bml_type
    integer :: N, M, i, iter, exp_order, method
    logical :: direct = .false.

    method = 0

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(trace(2))   

    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, p2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, d_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, wtmp_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, y_bml)
    if (method .eq. 1) then
       call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, I_bml)
       call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, ai_bml)
    endif
    cg_tol = 0.0001 !Squared Frobenius norm of residual matrix 
    occErr = 1.0_dp
    iter = 0
    exp_order = 2**nsteps
    cnst = beta/(4*exp_order)

    do while ((osteps .eq. 0 .and. occErr .gt. occErrLimit) .or. &
         (osteps .gt. 0 .and. iter .lt. osteps))
       iter = iter + 1

       ! Normalization 
       ! P0 = 0.5*II - cnst*(H0-mu0*II)
       call bml_copy(h_bml, p_bml)
       call prg_normalize_implicit_fermi(p_bml, cnst, mu)

       ! Do direct expansion
       if (direct .eqv. .true.) then 
          call prg_direct_exp(p_bml, y_bml, p2_bml, d_bml, wtmp_bml, exp_order, cg_tol, threshold)
          ! Do recursive expansion
       else if (method .eq. 0) then   
          write(*,*) "Doing CG"
          do i = 1, nsteps

             call bml_multiply_x2(p_bml, p2_bml, threshold, trace)

             ! Y = 2*(P2-P) + II
             call bml_copy(p2_bml, y_bml)
             call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
             call bml_scale_add_identity(y_bml, 2.0_dp, 1.0_dp, threshold)

             ! Solve Y*P = P2
                call prg_conjgrad(y_bml, p_bml, p2_bml, d_bml, wtmp_bml, cg_tol, threshold)
          enddo
       else 
          write(*,*) "Doing NS"
          do i = 1, nsteps

             call bml_multiply_x2(p_bml, p2_bml, threshold, trace)

             ! Y = 2*(P2-P) + II
             call bml_copy(p2_bml, y_bml)
             call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
             call bml_scale_add_identity(y_bml, 2.0_dp, 1.0_dp, threshold)
             if (i .eq. 1) then 
                call prg_conjgrad(y_bml, ai_bml, I_bml, d_bml, wtmp_bml, 0.9_dp, threshold)
             end if
             call prg_newtonschulz(y_bml, ai_bml, d_bml, wtmp_bml, threshold)
             call bml_multiply(ai_bml, p2_bml, p_bml, 1.0_dp, 0.0_dp, threshold)
          enddo

       end if
       !call bml_print_matrix("density matrix",p_bml,0,10,0,10)
       trdPdmu = bml_trace(p_bml)
       trP0 = trdPdmu
       trdPdmu = trdPdmu - bml_sum_squares(p_bml) ! sum p(i,j)**2

       trdPdmu = beta * trdPdmu
       mu = mu + (nocc - trP0)/trdPdmu
       occErr = abs(trP0 - nocc)
    enddo

    ! Adjust occupation
    ! X = II-P0
    call bml_copy(p_bml, d_bml)
    call bml_scale_add_identity(d_bml, -1.0_dp, 1.0_dp, threshold)

    call bml_multiply(p_bml, d_bml, wtmp_bml, 1.0_dp, 0.0_dp, threshold)
    ofactor = ((nocc - trP0)/trdPdmu) * beta
    call bml_add(p_bml, wtmp_bml, 1.0_dp, ofactor, threshold)
    call bml_print_matrix("density matrix",p_bml,0,10,0,10)

    deallocate(trace)

    call bml_deallocate(p2_bml)
    call bml_deallocate(wtmp_bml)
    call bml_deallocate(d_bml)
    call bml_deallocate(y_bml)
    if (method .eq. 1) then
       call bml_deallocate(ai_bml)
       call bml_deallocate(I_bml)
    endif

  end subroutine prg_implicit_fermi

  subroutine prg_direct_exp(p_bml, y_bml, p2_bml, d_bml, wtmp_bml, exp_order, cg_tol, threshold) 

    implicit none 

    type(bml_matrix_t), intent(inout) :: p_bml, y_bml, p2_bml, d_bml, wtmp_bml
    real(dp), intent(in) :: cg_tol, threshold 
    integer, intent(in) :: exp_order
    character(20) :: bml_type
    integer :: N, M 
    real(dp), allocatable :: trace(:)
    type(bml_matrix_t) :: x_bml

    allocate(trace(2))

    bml_type = bml_get_type(p_bml)
    N = bml_get_N(p_bml)
    M = bml_get_M(p_bml) 

    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x_bml)

    call bml_copy(p_bml, y_bml)
    call bml_copy(p_bml, x_bml)

    call bml_scale_add_identity(y_bml, -1.0_dp, 1.0_dp, threshold)

    call prg_exp_by_squaring(y_bml, p2_bml, exp_order, threshold, trace)
    call prg_exp_by_squaring(x_bml, p2_bml, exp_order, threshold, trace)  

    call bml_add(y_bml, x_bml, 1.0_dp, 1.0_dp, threshold)

    call prg_conjgrad(y_bml, p_bml, p2_bml, d_bml, wtmp_bml, cg_tol, threshold)        

    call bml_deallocate(x_bml)
    deallocate(trace)

  end subroutine prg_direct_exp

  subroutine prg_exp_by_squaring(x_bml, x2_bml, exp_order, threshold, trace) 

    implicit none

    type(bml_matrix_t), intent(inout) :: x_bml, x2_bml
    real(dp), intent(in) :: threshold
    integer, intent(in) :: exp_order
    real(dp), allocatable, intent(inout) :: trace(:)
    integer :: k 

    k = exp_order

    do while (k.gt.1)
       call bml_multiply_x2(x_bml, x2_bml, threshold, trace)
       k = k/2 
       call bml_copy(x2_bml, x_bml)
    enddo

  end subroutine prg_exp_by_squaring

  subroutine prg_newtonschulz(a_bml, ai_bml, r_bml, tmp_bml, threshold)

    implicit none

    type(bml_matrix_t), intent(inout) :: ai_bml, r_bml, tmp_bml
    type(bml_matrix_t), intent(in) :: a_bml
    real(dp), intent(in) :: threshold
    real(dp) :: norm
    integer :: i     

    norm = 1.0
    do while(norm > 0.1)
       call bml_copy(ai_bml, tmp_bml)
       call bml_multiply(a_bml, ai_bml, r_bml, 1.0_dp, 0.0_dp, threshold)
       call bml_scale_add_identity(r_bml, -1.0_dp, 1.0_dp, threshold)
       norm = bml_fnorm(r_bml)
       call bml_multiply(tmp_bml, r_bml, ai_bml, 1.0_dp, 1.0_dp, threshold)
       !  write(*,*) "norm = ", norm
    enddo

  end subroutine prg_newtonschulz

  subroutine prg_pcg(A_bml, p_bml, p2_bml, d_bml, wtmp_bml, cg_tol, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: A_bml
    type(bml_matrix_t), intent(inout) :: p_bml, p2_bml, d_bml, wtmp_bml
    real(dp), intent(in) :: cg_tol, threshold

    type(bml_matrix_t) :: M_bml, z_bml
    real(dp), allocatable :: diagonal(:)
    real(dp) :: alpha, beta
    character(20) :: bml_type
    integer :: k,N,M
    real(dp) :: r_norm_old, r_norm_new

    bml_type = bml_get_type(p_bml)
    N = bml_get_N(p_bml)
    M = bml_get_M(p_bml)

    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, z_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, M_bml)

    allocate(diagonal(N))
    call bml_get_diagonal(A_bml, diagonal)
    do k = 1,N
       diagonal(k) = 1.0_dp/diagonal(k)
    enddo
    call bml_set_diagonal(M_bml, diagonal)

    call bml_multiply(A_bml, p_bml, p2_bml, -1.0_dp, 1.0_dp, threshold)
    call bml_multiply(M_bml, p2_bml, z_bml, 1.0_dp, 0.0_dp, threshold) 
    r_norm_new = bml_traceMult(z_bml, p2_bml)
    call bml_copy(z_bml, d_bml)
    k = 0

    do while (bml_sum_squares(p2_bml) .gt. cg_tol)

       write(*,*) "r_norm", bml_sum_squares(p2_bml)
       k = k + 1
       if (k .ne. 1) then 
          beta = r_norm_new/r_norm_old
          call bml_add(d_bml, z_bml, beta, 1.0_dp, threshold)
       endif

       call bml_multiply(A_bml, d_bml, wtmp_bml, 1.0_dp, 0.0_dp, threshold)
       alpha = bml_traceMult(p2_bml,z_bml)/bml_traceMult(d_bml, wtmp_bml)
       call bml_add(p_bml, d_bml, 1.0_dp, alpha, threshold)
       call bml_add(p2_bml, wtmp_bml, 1.0_dp, -alpha, threshold)
       call bml_multiply(M_bml, p2_bml, z_bml, 1.0_dp, 0.0_dp, threshold)
       r_norm_old = r_norm_new
       r_norm_new = bml_traceMult(p2_bml,z_bml)
       if (k .gt. 100) then
          write(*,*) "PCG is not converging"
          stop
       endif
    enddo
    write(*,*) "Number of iterations:", k

    call bml_deallocate(z_bml)
    call bml_deallocate(M_bml)
    deallocate(diagonal)

  end subroutine prg_pcg

  subroutine prg_conjgrad(A_bml, p_bml, p2_bml, d_bml, wtmp_bml, cg_tol, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: A_bml
    type(bml_matrix_t), intent(inout) :: p_bml, p2_bml, d_bml, wtmp_bml
    real(dp), intent(in) :: cg_tol, threshold

    real(dp) :: alpha, beta
    integer :: k
    real(dp) :: r_norm_old, r_norm_new

    call bml_multiply(A_bml, p_bml, p2_bml, -1.0_dp, 1.0_dp, threshold)
    r_norm_new = bml_sum_squares(p2_bml)
    k = 0

    do while (r_norm_new .gt. cg_tol) 

       !   write(*,*) r_norm_new
       k = k + 1
       if (k .eq. 1) then 
          call bml_copy(p2_bml, d_bml)
       else 
          beta = r_norm_new/r_norm_old
          call bml_add(d_bml, p2_bml, beta, 1.0_dp, threshold)
       endif

       call bml_multiply(A_bml, d_bml, wtmp_bml, 1.0_dp, 0.0_dp, threshold)
       alpha = r_norm_new/bml_traceMult(d_bml, wtmp_bml)
       call bml_add(p_bml, d_bml, 1.0_dp, alpha, threshold)
       call bml_add(p2_bml, wtmp_bml, 1.0_dp, -alpha, threshold)
       r_norm_old = r_norm_new
       r_norm_new = bml_sum_squares(p2_bml)
       if (k .gt. 100) then
          write(*,*) "Conjugate gradient is not converging"
          stop
       endif
    enddo
    write(*,*) "Number of iterations:", k
  end subroutine prg_conjgrad

  subroutine prg_test_density_matrix(ham_bml, p_bml, beta, mu, nocc, osteps, occErrLimit, threshold)

    implicit none 

    type(bml_matrix_t), intent(in) :: ham_bml
    type(bml_matrix_t), intent(inout) :: p_bml
    real(dp), intent(in) ::  beta, nocc, occErrLimit, threshold
    real(dp), intent(inout) :: mu
    integer, intent(in) :: osteps
    character(20) :: bml_type
    integer :: N, i, iter
    real(dp) :: trdPdmu, trP0, ofactor, occErr
    real(dp), allocatable ::  eigenvalues(:)
    type(bml_matrix_t) :: eigenvectors_bml,occupation_bml,aux_bml,aux1_bml,i_bml

    bml_type = bml_get_type(p_bml)
    N = bml_get_N(p_bml)

    allocate(eigenvalues(N))

    call bml_zero_matrix(bml_type,bml_element_real,dp,N,N,eigenvectors_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,N,occupation_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,N,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,N,aux1_bml)
    call bml_identity_matrix(bml_type,bml_element_real,dp,N,N,i_bml)

    occErr = 1.0_dp
    iter = 0

    do while ((osteps .eq. 0 .and. occErr .gt. occErrLimit) .or. &
         (osteps .gt. 0 .and. iter .lt. osteps))
       iter = iter + 1

       call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

       do i=1,N   
          eigenvalues(i) = fermi(eigenvalues(i),mu,beta)
       enddo

       call bml_set_diagonal(occupation_bml, eigenvalues) 
       call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp,threshold)
       call bml_transpose(eigenvectors_bml, aux1_bml)
       call bml_multiply(aux_bml, aux1_bml, p_bml, 1.0_dp, 0.0_dp, threshold)
       !call bml_print_matrix("density matrix",p_bml,0,10,0,10)    
          trdPdmu = bml_trace(p_bml)
          trP0 = trdPdmu
         trdPdmu = trdPdmu - bml_sum_squares(p_bml) ! sum p(i,j)**2

        trdPdmu = beta * trdPdmu
        mu = mu + (nocc - trP0)/trdPdmu
        occErr = abs(trP0 - nocc)
    enddo

    ! Adjust occupation
    !! X = II-P0
    call bml_copy(i_bml, aux_bml)
    call bml_add(aux_bml, p_bml, 1.0_dp, -1.0_dp, threshold)

    call bml_multiply(p_bml, aux_bml, aux1_bml, 1.0_dp, 0.0_dp, threshold)
    ofactor = ((nocc - trP0)/trdPdmu) * beta
    call bml_add(p_bml, aux1_bml, 1.0_dp, ofactor, threshold)
   ! call bml_print_matrix("density matrix",p_bml,0,10,0,10)

    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(occupation_bml)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)
    call bml_deallocate(i_bml)

    deallocate(eigenvalues)

  end subroutine prg_test_density_matrix

  !> Gives the Fermi distribution value for energy e.
  !! \param e Energy.
  !! \param mu Fermi energy.
  !!
  real(dp) function fermi(e,mu,beta)

    real(dp), intent(in) :: e, mu, beta

    fermi = 1.0_dp/(1.0_dp+exp(beta*(e-mu)))

  end function fermi

end module prg_implicit_fermi_mod
