!> Auxiliary modules.
!! \brief This module will be used to have auxiliary routines.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_mod

#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif

  private

  integer, parameter :: dp = kind(1.0d0)
  logical :: err

  public :: gpmdcov_musearch, gpmdcov_musearch_bisec
  public :: gpmdcov_muDyn, gpmdcov_muFromParts, algo
  public :: gpmdcov_fermifunction

contains

  subroutine gpmdcov_fermifunction(beta,energies,mu,fermi)
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: energies(:)
    real(dp), intent(in) :: mu
    real(dp), intent(out) :: fermi(:)
    real(dp) :: biggest
    real(dp), allocatable :: exparg(:)

    if(.not.allocated(exparg))then
      allocate(exparg(size(energies)))
    endif

    biggest = log(HUGE(biggest))-1.0_dp
    exparg = beta*(energies(:)-mu)
    exparg = merge(exparg,biggest,exparg.le.biggest)
    fermi = 1.0_dp/(exp(exparg) + 1.0_dp)

  end subroutine gpmdcov_fermifunction

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
    use gpmdcov_dos_mod
    integer :: mycount
    real(dp) :: mls_mu

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

#ifdef DO_MPI
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
      do iptt = 1,gpat%TotalParts
        ipt = iptt
#endif
        do i = 1, norbsInEachCHAtRank(iptt)
          evalsInRank(i+shift) = syprt(ipt)%estr%evals(i)
          dvalsInRank(i+shift) = syprt(ipt)%estr%dvals(i)
          mycount = mycount + 1
        enddo
        shift = shift + norbsInEachCHAtRank(iptt)
      enddo

      !    write(*,*)"norbsInRank",norbsInRank
      !    write(*,*)"evalsInEachCHAtRank",norbsInEachCHAtRank

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

#ifdef DO_MPI
      call allGatherVIntParallel(norbsInEachCHAtRank, partsInEachRank(myRank),norbsInEachCH ,npartsVect, displ)
#else
      norbsInEachCH = norbsInEachCHAtRank
#endif

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

#ifdef DO_MPI
      call allGatherVRealParallel(evalsInRank, norbsInRank, evalsAll ,norbsInEachRank, displ)
      call allGatherVRealParallel(dvalsInRank, norbsInRank, dvalsAll ,norbsInEachRank, displ)
#else
      evalsAll = evalsInRank
      dvalsAll = dvalsInRank
#endif


      nocc = bndfilTotal*real(sy%estr%norbs,dp)

      err = .false.
#ifdef USE_NVTX
      call nvtxStartRange("gpmdcov_musearch",2)
#endif

      mls_mu = mls()
      !call gpmdcov_musearch(evalsAll,dvalsAll,beta,nocc,10000,10d-10,&
      !     &Ef,HOMO,LUMO,err,myRank,egap_glob,lt%verbose)
      

      !if(err .eqv. .true.)then
        call gpmdcov_musearch_bisec(evalsAll,dvalsAll,beta,nocc,10d-10,Ef,myRank,egap_glob,lib_mode,err_status,lt%verbose)
      !endif

#ifdef USE_NVTX
      call nvtxEndRange
#endif
      
      if(err_status) return

      call gpmdcov_msI("gpmdcov_muFromParts","Chemical potential (Mu or Ef) ="//to_string(Ef),lt%verbose,myRank)
      call gpmdcov_msI("gpmdcov_muFromParts","Time for mu search="//to_string(mls() - mls_mu),lt%verbose,myRank)

      if (estrout%write_tdos .and. MDstep .eq. 1) then
       if (gpat%TotalParts .gt. 1) then
           call gpmdcov_msI("gpmdcov_muFromParts","For computing DOS, TotalParts should be 1 ",lt%verbose,myRank)
           stop
       endif
       call compute_dos(estrout%tdos_num_points, estrout%tdos_sigma, estrout%tdos_emin,&
              &estrout%tdos_emax, evalsAll, Ef,  estrout%tdos_output_filename)
       if (estrout%compute_pdos) then
           call gpmdcov_msI("gpmdcov_muFromParts","Computing PDOS ",lt%verbose,myRank)
           !write(*,*) "Hindex ",syprt(1)%estr%hindex
           call compute_local_dos(estrout%tdos_num_points, estrout%pdos_atoms, syprt(1)%estr%hindex, estrout%tdos_sigma,&
                   &estrout%tdos_emin, estrout%tdos_emax,syprt(1)%estr%evects,evalsAll,syprt(1)%estr%over, Ef,&
                   &syprt(1)%estr%zmat, syprt(1)%symbol, estrout%pdos_output_filename)
        endif
        write(*,*) "gpmdcov_mod: called compute_dos to compute TDOS"
        stop
      endif

    end subroutine gpmdcov_muFromParts

    !> Perform musearch.
    !! \param evals Eigenvalues of the system.
    !! \param dvals Partial traces.
    !! \prarm beta Electronic temperature.
    !! \param maxMuIter Maximum numner of iterations.
    !! \param occTol Occupation tolerance.
    !! \param mu Chemical potential.
    !!
    subroutine gpmdcov_musearch(evals, dvals, beta, nocc, maxMuIter, occTol, mu, &
         & HOMO, LUMO, err, rank, egap, verbose)
      use gpmdcov_writeout_mod
      implicit none
      integer                ::  i, j, norbs
      integer, intent(in)    ::  maxMuIter, rank
      real(dp)               ::  den, occErr, mu0, occ, muMax,muMin
      real(dp), allocatable  ::  fvals(:)
      real(dp), intent(in)   ::  evals(:), dvals(:)
      real(dp), intent(in)   ::  occTol, beta, nocc
      real(dp), intent(inout) :: egap
      real(dp), intent(inout)   ::  mu, HOMO, LUMO
      logical, intent(inout) :: err
      integer, optional, intent(in) :: verbose

      mu0 = 0.5d0*(evals(int(floor(Nocc))) + evals(int(floor(Nocc))+1))
      muMin = minval(evals)
      muMax = maxval(evals)
      if (present(verbose)) then
        call gpmdcov_msI("gpmdcov_musearch","In gpmdcov_musearch with muMin,muMax,mu0 &
             &= "//to_string(muMin)//","//to_string(muMax)//","//to_string(mu0),verbose,rank)
      endif
      norbs = size(evals, dim = 1)

      if(.not.allocated(fvals))then
        allocate(fvals(norbs))
      endif

      !  if(mu > (muMin - muMax)/2.0_dp)then
      !      mu = (muMin - muMax)/2.0_dp
      !  else
      !      mu0 = mu !(muMin - muMax)/2.0_dp
      !  endif

      call gpmdcov_fermifunction(beta,evals,mu0,fvals)

      !do j = 1, norbs
      !  fvals(j) = 1.0_dp/(exp(beta*(evals(j)-mu0))+1.0_dp)
      !end do

      occ = 0.0_dp
      do j = 1, norbs
        occ = occ + fvals(j)*dvals(j)
      end do

      occErr = abs(occ-nocc)

      do i = 1, maxMuIter
         
#ifdef USE_NVTX
      call nvtxStartRange("gpmdcov_fermifunction",3)
#endif

        call gpmdcov_fermifunction(beta,evals,mu0,fvals)
        !do j = 1, norbs
        !  fvals(j) = 1.0_dp/(exp(beta*(evals(j)-mu0))+1.0_dp)
        !end do
#ifdef USE_NVTX
      call nvtxEndRange
#endif

        occ = 0.0_dp
        den = 0.0_dp
        do j = 1, norbs
          occ = occ + fvals(j)*dvals(j)
          den = den + beta*fvals(j)*(1.0_dp-fvals(j))*dvals(j)
        end do

        occErr = abs(occ - nocc)
        if (occErr < occTol) then
          mu = mu0
          exit
        elseif (den > 1.0d-12) then
          mu0 = mu0 + (nocc - occ)/den
          mu = mu0
        endif

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

      HOMO = muMin
      LUMO = muMax
      do i = 1,norbs
        if(evals(i) .lt. mu)then
          HOMO = max(evals(i),HOMO)
        elseif(evals(i) .gt. mu)then
          LUMO = min(evals(i),LUMO)
        endif
      enddo

      egap = LUMO-HOMO

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
    subroutine gpmdcov_musearch_bisec(evals,dvals,beta,noc,tol,mu,rank,egap,lib_mode,err_status,verbose)
      use gpmdcov_writeout_mod
      integer                  ::  i, m
      integer                  ::  norb
      integer, intent(in)      ::  rank
      integer, optional, intent(in) :: verbose
      real(dp)                 ::  Ft1, Ft2, Prod, egap
      real(dp)                 ::  step, tol, fermi, mumax, muMin, HOMO, LUMO
      real(dp), intent(in)     ::  noc, evals(:), dvals(:), beta
      real(dp)  ::  mu
      logical :: err_status,lib_mode
      real(dp), allocatable :: fvals(:)

      if(.not.allocated(fvals))then
        allocate(fvals(size(evals)))
      endif

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
      call gpmdcov_fermifunction(beta,evals,mu,fvals)

      do i=1,norb
        fermi = fvals(i)
        ft1 = ft1 + 1.0_dp*fermi*dvals(i)
      enddo
      ft1=ft1-noc

      do m=1,1000001
        if(m.gt.1000000)then
          write(*,*)"Bisection method in gpmdcov_musearch_bisec not converging ..."
          if(lib_mode)then
            err_status = .true.
            return
          else
            stop
          endif
        endif
        if(mu > muMax + 1.0_dp .or. mu < muMin - 1.0_dp)then
          write(*,*)"Bisection method is diverging"
          write(*,*)"muMin=",muMin,"muMax=",muMax
          write(*,*)evals
          stop
        endif

        if(abs(ft1).lt.tol)then !tolerance control
          exit
        endif

        mu = mu + step

        ft2=0.0_dp

        !New sum of the occupations
        call gpmdcov_fermifunction(beta,evals,mu,fvals)
        do i=1,norb
          fermi = fvals(i)
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

      HOMO = muMin
      LUMO = muMax
      do i = 1,norb
        if(evals(i) .lt. mu)then
          HOMO = max(evals(i),HOMO)
        elseif(evals(i) .gt. mu)then
          LUMO = min(evals(i),LUMO)
        endif
      enddo

      !write(*,*)"EGAP BISEC",LUMO-HOMO,mu
      egap = LUMO-HOMO

    end subroutine gpmdcov_musearch_bisec

  end module gpmdcov_mod
