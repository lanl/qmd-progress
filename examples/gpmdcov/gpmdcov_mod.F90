!> Auxiliary modules.
!! \brief This module will be used to have auxiliary routines.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  integer :: mycount
  logical :: err

  public :: gpmdcov_musearch, gpmdcov_musearch_bisec
  public :: gpmdcov_muDyn, gpmdcov_muFromParts

contains

  subroutine gpmdcov_muDyn(nguess,Nr_SCF)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    real(dp), allocatable :: nguess(:)
    integer, intent(in) :: Nr_SCF

     !Actualize the Fermi level dynamically

        if(.not.allocated(acceprat))then
          allocate(acceprat(2))
          acceprat = 0
          Efstep = 0.1_dp
        endif

        acceprat(2) = acceprat(1)
        acceprat(1) = sign(1.0_dp,tch)

        if(acceprat(2)*acceprat(1) < 0)then
          Efstep = Efstep*0.8_dp
        else
          Efstep = Efstep*1.01_dp
        endif

        if(Nr_SCF.gt.10)then
          if(iscf.gt.10)Ef = Ef -  sign(1.0_dp,tch)*min(tch**2,Efstep)
        else
          Ef = Ef -  sign(1.0_dp,tch)*min(tch**2,Efstep)
        endif

        !Normalize charges to tch
        nguess(:) = nguess(:) - tch/real(sy%nats)

       call gpmdcov_msI("gpmdcov_muDyn","Chemical potential (Mu or Ef) ="//to_string(Ef),lt%verbose,myRank)

  end subroutine gpmdcov_muDyn

  subroutine gpmdcov_muFromParts()
    
    use gpmdcov_vars
    use gpmdcov_writeout_mod

    call gpmdcov_msI("gpmdcov_getmu","In gpmdcov_getmu ...",lt%verbose,myRank)

    norbsInRank = sum(norbsInEachCHAtRank)

    if(allocated(evalsInRank))then
      deallocate(evalsInRank)
      deallocate(dvalsInRank)
    endif

    allocate(evalsInRank(norbsInRank))
    allocate(dvalsInRank(norbsInRank))

    evalsInRank = 0.0_dp
    dvalsInRank = 0.0_dp
    shift = 0
    mycount = 0
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
      do i = 1, norbsInEachCHAtRank(iptt)
        evalsInRank(i+shift) = syprt(ipt)%estr%aux(1,i)
        dvalsInRank(i+shift) = syprt(ipt)%estr%aux(2,i)
        mycount = mycount + 1
      enddo
      shift = shift + norbsInEachCHAtRank(iptt)
    enddo

    if(.not.allocated(norbsInEachCH)) allocate(norbsInEachCH(gpat%TotalParts))
    norbsInEachCH = 0
    nRanks = getNRanks()
    if(.not.allocated(npartsVect))then
      allocate(npartsVect(nRanks))
      allocate(displ(nRanks))
      allocate(norbsInEachRank(nRanks))
    endif

    do i = 1,nRanks
      npartsVect(i) = partsInEachRank(i)
      if(i == 1)then
        shift = 0
      else
        shift = shift + partsInEachRank(i-1)
      endif
      displ(i) = shift
    enddo

    call allGatherVIntParallel(norbsInEachCHAtRank, partsInEachRank(myRank),norbsInEachCH ,npartsVect, displ)

    kk = 0

    norbsInEachRank = 0.0_dp
    do i = 1,nRanks
      do j = 1,partsInEachRank(i)
        kk = kk + 1
        norbsInEachRank(i) = norbsInEachRank(i) + norbsInEachCH(kk)
      enddo
      if(i == 1)then
        shift = 0
      else
        shift = shift + norbsInEachRank(i-1)
      endif
      displ(i) = shift
    enddo
    totalNorbs = sum(norbsInEachRank)

    if(allocated(evalsAll))then
      deallocate(evalsAll)
      deallocate(dvalsAll)
    endif
    allocate(evalsAll(totalNorbs))
    allocate(dvalsAll(totalNorbs))
    evalsAll = 0.0_dp
    dvalsAll = 0.0_dp

    call allGatherVRealParallel(evalsInRank, norbsInRank, evalsAll ,norbsInEachRank, displ)
    call allGatherVRealParallel(dvalsInRank, norbsInRank, dvalsAll ,norbsInEachRank, displ)
    
    if(minval(dvalsAll) < 0)then
      write(*,*)minval(dvalsAll)
      write(223,*)dvalsAll
      stop
    endif

    nocc = bndfilTotal*real(sy%estr%norbs,dp)

    err = .false.
    call gpmdcov_musearch(evalsAll,dvalsAll,beta,nocc,10000,10d-10,Ef,err,myRank,lt%verbose)
    if(err .eqv. .true.)then
      call gpmdcov_musearch_bisec(evalsAll,dvalsAll,beta,nocc,10d-10,Ef,myRank,lt%verbose)
    endif

    call gpmdcov_msI("gpmdcov_muFromParts","Chemical potential (Mu or Ef) ="//to_string(Ef),lt%verbose,myRank)

  end subroutine gpmdcov_muFromParts

  !> Perform musearch.
  !! \param evals Eigenvalues of the system.
  !! \param dvals Partial traces.
  !! \prarm beta Electronic temperature.
  !! \param maxMuIter Maximum numner of iterations.
  !! \param occTol Occupation tolerance.
  !! \param mu Chemical potential.
  !!
  subroutine gpmdcov_musearch(evals, dvals, beta, nocc, maxMuIter, occTol, mu, err, rank,verbose)
    use gpmdcov_writeout_mod
    implicit none
    integer                ::  i, j, norbs
    integer, intent(in)    ::  maxMuIter, rank
    real(dp)               ::  den, occErr, mu0, occ, muMax,muMin
    real(dp), allocatable  ::  fvals(:)
    real(dp), intent(in)   ::  evals(:), dvals(:)
    real(dp), intent(in)   ::  occTol, beta, nocc
    real(dp), intent(inout)   ::  mu
    logical, intent(inout) :: err
    integer, optional, intent(in) :: verbose


    if (present(verbose)) then
      call gpmdcov_msI("gpmdcov_musearch","In gpmdcov_musearch ...",verbose,rank)
    endif

    muMin=minval(evals)
    muMax = maxval(evals)

    norbs = size(evals, dim = 1)
    allocate(fvals(norbs))

    mu0 = mu

    do j = 1, norbs
      fvals(j) = 1.0_dp/(exp(beta*(evals(j)-mu0))+1.0_dp)
    end do

    occ = 0.0_dp
    do j = 1, norbs
      occ = occ + fvals(j)*dvals(j)
    end do

    occErr = abs(occ-nocc)

    do i = 1, maxMuIter

      do j = 1, norbs
        fvals(j) = 1.0_dp/(exp(beta*(evals(j)-mu0))+1.0_dp)
      end do

      occ = 0.0_dp
      den = 0.0_dp
      do j = 1, norbs
        occ = occ + fvals(j)*dvals(j)
        den = den + beta*fvals(j)*(1.0_dp-fvals(j))*dvals(j)
      end do

      occErr = abs(occ - nocc)
      mu0 = mu0 + (nocc - occ)/den
      mu = mu0

      if(mu > muMax .or. mu < muMin)then
        err = .true.
        mu = 0.0_dp
        return
      endif

      if (i .eq. maxMuIter) then
        write(*,*) "Max mu iters reached for NR. Trying bisection ..."
        err = .true.
        mu = 0.0_dp
        return
      else
        if(occErr .lt. occTol) exit
      end if

    end do

    deallocate(fvals)

  end subroutine gpmdcov_musearch


  !> Routine to compute the Fermi level given a set of eigenvalues and a temperature.
  !! It applies the Bisection method over the function:
  !! \f$ g(\mu) = \sum_k 2 f(\epsilon_k - \mu) - N = 0 \f$
  !! Where \f$ f(\epsilon_k - \mu) = \frac{1}{1+\exp{(\epsilon_k - \mu)/(k_bT)}}\f$.
  !! \param evals Eigenvalues of the system (\f$ \{ \epsilon_k \} \f$).
  !! \param kbt Temperature times the Boltzmann's constant  (\f$ k_bT  \f$).
  !! \param bndfil Filing factor (\f$ N_{el}/(2*N_{orbs})\f$).
  !! \param tol Tolerance for the bisection method.
  !! \param Ef Fermi level (\f$ \mu \f$).
  !!
  subroutine gpmdcov_musearch_bisec(evals,dvals,beta,noc,tol,mu,rank,verbose)
    use gpmdcov_writeout_mod
    integer                  ::  i, j, k, m
    integer                  ::  norb
    integer, intent(in)      ::  rank
    integer, optional, intent(in) :: verbose
    real(dp)                 ::  Ft1, Ft2, Kb, Prod, divergTol
    real(dp)                 ::  T, step, tol, nel, fermi, mumax, muMin
    real(dp), intent(in)     ::  noc, evals(:), dvals(:), beta
    real(dp)  ::  mu

    if (present(verbose)) then
      call gpmdcov_msI("gpmdcov_musearch_bisec","In gpmdcov_musearch_bisec ...",verbose,rank)
    endif

    norb = size(evals,dim=1)
    muMin=minval(evals)
    muMax = maxval(evals)
    mu = muMin
    step=abs(muMax-muMin)
    Ft1=0.0_dp
    Ft2=0.0_dp
    prod=0.0_dp

    !Sum of the occupations
    do i=1,norb
      fermi = 1.0_dp/(exp(beta*(evals(i)-mu))+1.0_dp)
      ft1 = ft1 + 1.0_dp*fermi*dvals(i)
    enddo
    ft1=ft1-noc

    do m=1,1000001

      if(m.gt.1000000)then
        write(*,*)"Bisection method in gpmdcov_musearch_bisec not converging ..."
        stop
      endif
      if(mu > muMax .or. mu < muMin)then
        write(*,*)"Bisection method is diverging"
        stop
      endif

      if(abs(ft1).lt.tol)then !tolerance control
        exit
      endif

      mu = mu + step

      ft2=0.0_dp

      !New sum of the occupations
      do i=1,norb
        fermi = 1.0_dp/(exp(beta*(evals(i)-mu))+1.0_dp)
        ft2 = ft2 + 1.0_dp*fermi*dvals(i)
      enddo

      ft2=ft2-noc

      !Product to see the change in sign.
      prod = ft2*ft1
      if(prod.lt.0)then
        mu=mu-step
        step=step/2.0_dp !If the root is inside we shorten the step.
      else
        ft1=ft2  !If not, Ef moves forward.
      endif

    enddo

  end subroutine gpmdcov_musearch_bisec

end module gpmdcov_mod
