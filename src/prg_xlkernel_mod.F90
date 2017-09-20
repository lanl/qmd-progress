!> Add name.
!! \ingroup PROGRESS
!! XL kernel ....
module prg_xlkernel_mod

  use bml
  use prg_kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  type, public :: xlk_type

     !> Kernel type
     character(20) :: kerneltype

     !> Verbosity level
     integer       :: verbose, nrank

     !> Coefficient for mixing
     real(dp)      :: scalecoeff

  end type xlk_type

  public :: prg_parse_xlkernel, prg_rank1, prg_Fermi
  public :: prg_v_kernel_Fermi, prg_kernel_Fermi_full
  private :: prg_MMult, prg_Eig, prg_get_deriv_finite_temp, prg_inv

contains

  !> The parser for the mixer routines.
  !!
  subroutine prg_parse_xlkernel(input,filename)

    type(xlk_type), intent(inout) :: input
    integer, parameter :: nkey_char = 1, nkey_int = 2, nkey_re = 1, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'KernelType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'Full']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Verbose=','Nrank=']
    integer :: valvector_int(nkey_int) = (/ &
         0, 1 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'ScaleCoeff=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.1 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'MixerON=']
    logical :: valvector_log(nkey_log) = (/&
         .true./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'XLKERNEL{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    input%kerneltype = valvector_char(1)

    !Integers
    input%verbose = valvector_int(1)
    input%nrank = valvector_int(2)

    !Logicals

    !Reals
    input%scalecoeff = valvector_re(1)

  end subroutine prg_parse_xlkernel

  subroutine prg_Fermi(D0,QQ,ee,gap,Fe_vec,mu0,H,Z,Nocc,T,OccErrLim,MaxIt,HDIM)

    integer, parameter             :: PREC = 8
    integer(PREC), intent(in)      :: HDIM, Nocc
    real(PREC), parameter          :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
    real(PREC), parameter          :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
    real(PREC), intent(in)         :: H(HDIM,HDIM), Z(HDIM,HDIM)
    real(PREC)                     :: H0(HDIM,HDIM), Fe(HDIM,HDIM)
    real(PREC), intent(out)        :: D0(HDIM,HDIM), QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM), gap
    real(PREC)                     :: X(HDIM,HDIM)
    real(PREC), intent(in)         :: T, OccErrLim
    real(PREC), intent(inout)      :: mu0
    integer(PREC), intent(in)      :: MaxIt
    real(PREC)                     :: OccErr, Occ, dOcc, mu_0, beta, Occ_I
    integer(PREC)                  :: I,ITER

    !! Orthogonalize H with inverse overlap factor Z, where Z'SZ = I
    call MMult(ONE,Z,H,ZERO,X,'T','N',HDIM)  !
    call MMult(ONE,X,Z,ZERO,H0,'N','N',HDIM)   ! H0 = Z'*H*Z
    call Eig(H0,QQ,ee,'V',HDIM)
    !! QQ eigenvectors of H0 and ee the eigenvalues

    OccErr = ONE
    Fe = ZERO
    beta = 1.D0/(kB*T)    ! Temp in Kelvin
    gap = ee(Nocc+1)-ee(Nocc)
    if (mu0 == 0.D0) then
       mu0 = 0.5D0*(ee(Nocc)+ee(Nocc+1)) ! Initial estimate of chemical potential
    endif
    ITER = 0
    do while (OccErr > 1e-9)
       ITER = ITER + 1
       Occ = ZERO
       dOcc = ZERO
       do I = 1,HDIM
          Occ_I = ONE/(exp(beta*(ee(I)-mu0))+ONE)
          Fe(I,I) = Occ_I
          Fe_vec(I) = Occ_I
          Occ = Occ + Occ_I
          dOcc = dOcc + beta*Occ_I*(ONE-Occ_I)
       enddo
       OccErr = abs(Nocc-Occ)
       if (abs(OccErr) > 1e-9) then
          mu0 = mu0 + (Nocc-Occ)/dOcc  ! Adjust occupation by shifting mu
       endif
       if(ITER.gt.MaxIt) then
          OccErr = ZERO
       endif
    enddo
    call MMult(ONE,QQ,Fe,ZERO,X,'N','N',HDIM)
    call MMult(ONE,X,QQ,ZERO,D0,'N','T',HDIM)

  end subroutine prg_Fermi

  subroutine prg_kernel_Fermi_full(KK,JJ,D0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
       MaxIt,eps, m,HDIM, Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
       H_INDEX_END,H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

    integer, parameter             :: PREC = 8
    integer(PREC), intent(in)      :: Nr_atoms, HDIM, Nocc,Max_Nr_Neigh
    real(PREC), parameter          :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
    real(PREC), parameter          :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
    real(PREC), intent(in)         :: Coulomb_acc, TIMERATIO
    real(PREC)                     :: v(Nr_atoms)
    real(PREC), intent(in)         :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
    integer(PREC), intent(in)      :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
    real(PREC), intent(in)         :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
    real(PREC), intent(in)         :: H(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM)
    real(PREC)                     :: H0(HDIM,HDIM), H1(HDIM,HDIM)
    real(PREC), intent(inout)      :: D0(HDIM,HDIM)
    real(PREC)                     :: X(HDIM,HDIM), XX(Nr_atoms,Nr_atoms)
    real(PREC), intent(in)         :: T, eps
    real(PREC), intent(inout)      :: mu0, mu1
    integer(PREC), intent(in)      :: m, MaxIt
    character(10), intent(in)      :: Element_Type(Nr_atoms)
    integer(PREC), intent(in)      :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)         :: nnRx(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)         :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
    real(PREC)                     :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Pot(Nr_atoms)
    real(PREC)                     :: Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms)
    real(PREC)                     :: Coulomb_Pot_Real_dq_v(Nr_atoms), Coulomb_Pot_dq_v(Nr_atoms)
    real(PREC)                     :: Coulomb_Force_Real_I(3),Coulomb_Force_k(3,Nr_atoms)
    real(PREC)                     :: D_dq_v(HDIM,HDIM), H_dq_v(HDIM,HDIM), dq_v(Nr_atoms)
    real(PREC)                     :: dq_dv(Nr_atoms)
    real(PREC), intent(out)        :: KK(Nr_atoms,Nr_atoms),JJ(Nr_atoms,Nr_atoms)
    real(PREC), intent(in)         :: QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM)
    integer(PREC)                  :: I,J,K, ITER, mm

    XX = ZERO
    X = ZERO
    H1 = ZERO
    Coulomb_Pot_dq_v = ZERO
    Coulomb_Pot_k = ZERO
    dq_v = ZERO
    H_dq_v = ZERO

    do J = 1,Nr_atoms  ! TRIVIAL OPEN_MP OR MPI PARALLELISM
       dq_v(J) = ONE
       do I = 1,Nr_atoms
          call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
               dq_v,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
          Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
       enddo
       call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,dq_v,Nr_atoms,Coulomb_acc,TIMERATIO,HDIM,Max_Nr_Neigh)
       Coulomb_Pot_dq_v = Coulomb_Pot_Real+Coulomb_Pot_k

       H_dq_v = ZERO
       do I = 1,Nr_atoms
          do K = H_INDEX_START(I),H_INDEX_END(I)
             H_dq_v(K,K) = Hubbard_U(I)*dq_v(I) + Coulomb_Pot_dq_v(I)
          enddo
       enddo
       call MMult(ONE,S,H_dq_v,ZERO,X,'N','N',HDIM)
       call MMult(ONE,H_dq_v,S,ONE,X,'N','T',HDIM)
       H_dq_v =  (ONE/TWO)*X

       call MMult(ONE,Z,H,ZERO,X,'T','N',HDIM)  !
       call MMult(ONE,X,Z,ZERO,H0,'N','N',HDIM)   ! H0 = Z'*H*Z

       call MMult(ONE,Z,H_dq_v,ZERO,X,'T','N',HDIM)  !
       call MMult(ONE,X,Z,ZERO,H1,'N','N',HDIM)   ! H1 = Z'*H_dq_v*Z

       call get_deriv_finite_temp(D_dq_v,H0,H1,Nocc,T,QQ,ee,Fe_vec,mu0,eps,HDIM)

       call MMult(TWO,Z,D_dq_v,ZERO,X,'N','N',HDIM)
       call MMult(ONE,X,Z,ZERO,D_dq_v,'N','T',HDIM)

       call MMult(ONE,D_dq_v,S,ZERO,X,'N','N',HDIM)
       dq_dv = ZERO
       do I = 1, Nr_atoms
          do K = H_INDEX_START(I), H_INDEX_END(I)
             dq_dv(I) = dq_dv(I) + X(K,K)
          enddo
          JJ(I,J) = dq_dv(I)
       enddo
       dq_v = ZERO
    enddo
    XX = JJ
    do i = 1,Nr_atoms
       XX(i,i) = XX(i,i) - ONE
    enddo
    call Inv(XX,KK,Nr_atoms)

  end subroutine prg_kernel_Fermi_full

  subroutine prg_v_kernel_Fermi(D0,dq_dv,v,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
       MaxIt,eps, m,HDIM, Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
       H_INDEX_END,H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

    integer, parameter             :: PREC = 8
    integer(PREC), intent(in)      :: Nr_atoms, HDIM, Nocc,Max_Nr_Neigh
    real(PREC), parameter          :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
    real(PREC), parameter          :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
    real(PREC), intent(in)         :: Coulomb_acc, TIMERATIO, v(Nr_atoms)
    real(PREC), intent(in)         :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
    integer(PREC), intent(in)      :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
    real(PREC), intent(in)         :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
    real(PREC), intent(in)         :: H(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM)
    real(PREC)                     :: H0(HDIM,HDIM), H1(HDIM,HDIM)
    real(PREC), intent(inout)      :: D0(HDIM,HDIM)
    real(PREC)                     :: X(HDIM,HDIM)
    real(PREC), intent(in)         :: T, eps
    real(PREC), intent(inout)      :: mu0, mu1
    integer(PREC), intent(in)      :: m, MaxIt
    character(10), intent(in)      :: Element_Type(Nr_atoms)
    integer(PREC), intent(in)      :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)         :: nnRx(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)         :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
    real(PREC)                     :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Pot(Nr_atoms)
    real(PREC)                     :: Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms)
    real(PREC)                     :: Coulomb_Pot_Real_dq_v(Nr_atoms), Coulomb_Pot_dq_v(Nr_atoms)
    real(PREC)                     :: Coulomb_Force_Real_I(3),Coulomb_Force_k(3,Nr_atoms)
    real(PREC)                     :: D_dq_v(HDIM,HDIM), H_dq_v(HDIM,HDIM), dq_v(Nr_atoms)
    real(PREC), intent(out)        :: dq_dv(Nr_atoms)
    real(PREC), intent(in)         :: QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM)
    integer(PREC)                  :: I,J,K, ITER, mm

    dq_dv = ZERO
    dq_v = v/norm2(v)

    do I = 1,Nr_atoms
       call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
            dq_v,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
       Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    enddo
    call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,dq_v,Nr_atoms,Coulomb_acc,TIMERATIO,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_dq_v = Coulomb_Pot_Real+Coulomb_Pot_k

    H_dq_v = ZERO
    do I = 1,Nr_atoms
       do K = H_INDEX_START(I),H_INDEX_END(I)
          H_dq_v(K,K) = Hubbard_U(I)*dq_v(I) + Coulomb_Pot_dq_v(I)
       enddo
    enddo
    call MMult(ONE,S,H_dq_v,ZERO,X,'N','N',HDIM)
    call MMult(ONE,H_dq_v,S,ONE,X,'N','T',HDIM)
    H_dq_v =  (ONE/TWO)*X

    call MMult(ONE,Z,H,ZERO,X,'T','N',HDIM)  !
    call MMult(ONE,X,Z,ZERO,H0,'N','N',HDIM)   ! H0 = Z'*H*Z

    call MMult(ONE,Z,H_dq_v,ZERO,X,'T','N',HDIM)  !
    call MMult(ONE,X,Z,ZERO,H1,'N','N',HDIM)   ! H1 = Z'*H_dq_v*Z

    call get_deriv_finite_temp(D_dq_v,H0,H1,Nocc,T,QQ,ee,Fe_vec,mu0,eps,HDIM)

    call MMult(TWO,Z,D_dq_v,ZERO,X,'N','N',HDIM)
    call MMult(ONE,X,Z,ZERO,D_dq_v,'N','T',HDIM)

    call MMult(ONE,D_dq_v,S,ZERO,X,'N','N',HDIM)
    do I = 1, Nr_atoms
       do K = H_INDEX_START(I), H_INDEX_END(I)
          dq_dv(I) = dq_dv(I) + X(K,K)
       enddo
    enddo

  end subroutine prg_v_kernel_Fermi

  subroutine  prg_get_deriv_finite_temp(P1,H0,H1,Nocc,T,Q,ev,fe,mu0,eps,HDIM)

    integer, parameter      :: PREC = 8
    integer(PREC), intent(in)       :: HDIM, Nocc
    real(PREC), parameter   :: ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0
    real(PREC), intent(in)  :: H0(HDIM,HDIM), H1(HDIM,HDIM), Q(HDIM,HDIM), ev(HDIM), fe(HDIM)
    real(PREC)              :: X(HDIM,HDIM), YY(HDIM,HDIM), fermifun
    real(PREC)              :: P02(HDIM,HDIM), II(HDIM,HDIM), XI(HDIM,HDIM)
    real(PREC)              :: DX1(HDIM,HDIM)
    real(PREC)              :: ID0(HDIM,HDIM), T12(HDIM,HDIM), DDT(HDIM,HDIM)
    real(PREC), intent(out) :: P1(HDIM,HDIM) ! DM derivative
    real(PREC), intent(in)  :: T, eps
    real(PREC), intent(inout)  :: mu0
    real(PREC)              ::  OccErr, TrdPdmu, TrP0, TrP1
    real(PREC)              :: beta, cnst, kB, mu1, dh, dm, y, dy, trX, trDDT
    integer(PREC)           :: N, Cnt, i, j, k,l, ind, ind2, fact

    N = HDIM
    kB = 8.61739e-5
    beta = 1.D0/(kB*T)     ! Temp in Kelvin

    ! T12 = F0_ort_Vecs'*F1_ort*F0_ort_Vecs;
    call MMult(ONE,Q,H1,ZERO,X,'T','N',HDIM)
    call MMult(ONE,X,Q,ZERO,T12,'N','N',HDIM)

    ! Divided differences matrix
    DDT = ZERO
    X = ZERO
    do i = 1,N
       do j = 1,N
          if (abs(ev(i)-ev(j)) < 1e-4) then
             ! divided difference from Taylor expansion
             y = (ev(i)+ev(j))/2.D0
             if (abs(beta*(y-mu0)) > 500.D0) then
                DDT(i,j) = 0.D0
             else
                DDT(i,j) = -beta*exp(beta*(y-mu0))/(exp(beta*(y-mu0))+1.D0)**2
             endif
          else
             DDT(i,j) = (fe(i)-fe(j))/(ev(i)-ev(j))
          endif
          X(i,j) = DDT(i,j)*T12(i,j)
       end do
    end do
    mu1 = ZERO
    trX = ZERO
    trDDT = ZERO
    do i = 1,N
       trX = trX + X(i,i)
       trDDT = trDDT + DDT(i,i)
    end do
    mu1 = trX/trDDT
    do i = 1,N
       X(i,i) = X(i,i)-DDT(i,i)*mu1
    end do
    call MMult(ONE,Q,X,ZERO,YY,'N','N',HDIM)
    call MMult(ONE,YY,Q,ZERO,P1,'N','T',HDIM)

  end subroutine prg_get_deriv_finite_temp

  subroutine prg_MMult(alpha,A,B,beta,C,TA,TB,HDIM)

    integer, parameter          :: PREC = 8
    integer(PREC), intent(in)   :: HDIM
    real(PREC), intent(in)      :: A(HDIM,HDIM), B(HDIM,HDIM), alpha, beta
    real(PREC), intent(inout)   :: C(HDIM,HDIM)
    character(1), intent(in)    :: TA, TB

    if (PREC.eq.4) then
       call SGEMM(TA, TB, HDIM, HDIM, HDIM, alpha, &
            A, HDIM, B, HDIM, beta, C, HDIM)
    else
       call DGEMM(TA, TB, HDIM, HDIM, HDIM, alpha, &
            A, HDIM, B, HDIM, beta, C, HDIM)
    endif

  end subroutine prg_MMult

  subroutine prg_Eig(A,Q,ee,type,HDIM)

    integer, parameter             :: PREC = 8
    integer(PREC), intent(in)      :: HDIM
    real(PREC), intent(in)         :: A(HDIM,HDIM)
    real(PREC), intent(out)        :: Q(HDIM,HDIM), ee(HDIM)
    character(1), intent(in)       :: type  ! 'N'(for eigenvaules only) or 'V' (otherwise)

    real(PREC)    :: NONO_EVALS(HDIM), NONO_WORK(1 + 6*HDIM + 2*HDIM*HDIM)
    real(PREC)    :: NONOTMP(HDIM,HDIM), Z(HDIM,HDIM)
    integer(PREC) :: NONO_IWORK(3+5*HDIM)
    integer(PREC) :: I, J, INFO, NONO_LWORK
    real(PREC)    :: INVSQRT,ERR_CHECK,NUMTHRESH

    NONO_LWORK = 1 + 6*HDIM + 2*HDIM*HDIM

    Q = A
    if (PREC.eq.4) then
       call SSYEV(type, 'U', HDIM, Q, HDIM, ee, NONO_WORK, &
            NONO_LWORK, INFO)
    else
       call DSYEV('V', 'U', HDIM, Q, HDIM, ee, NONO_WORK, &
            NONO_LWORK, INFO)
    endif

  end subroutine prg_Eig

  subroutine prg_inv(X,XI,HDIM)

    integer, parameter        :: PREC = 8
    integer(PREC), intent(in) :: HDIM
    integer(PREC)             :: LDA, LWORK, M, N, INFO, IPIV(HDIM)
    real(PREC), intent(in)    :: X(HDIM,HDIM)
    real(PREC), intent(out)   :: XI(HDIM,HDIM)
    real(PREC)                :: WORK(HDIM*HDIM)

    LDA = HDIM
    M = HDIM
    N = HDIM
    LWORK = HDIM*HDIM
    XI = X

    if (PREC.eq.4) then
       call SGETRF(M, N, XI, LDA, IPIV, INFO)
       call SGETRI(N, XI, N, IPIV, WORK, LWORK, INFO)
    else
       call DGETRF(M, N, XI, LDA, IPIV, INFO)
       call DGETRI(N, XI, N, IPIV, WORK, LWORK, INFO)
    endif

  end subroutine prg_inv

  !> Rank1 kernel ....
  !! \param param1 ..
  !! \param verbose Different levels of verbosity.
  subroutine prg_rank1(verbose)

    integer :: i,j,info,s,k,n
    integer :: nint, na, nb
    integer, intent(in) :: verbose
    real(dp), allocatable :: d(:),dl(:),dnewin(:)
    real(dp), allocatable :: dnewout(:)
    real(dp), allocatable :: coef(:,:),b(:),ipiv(:)

    write(*,*)"prt_rank1 Verbosity",verbose

  end subroutine prg_rank1

end module prg_xlkernel_mod
