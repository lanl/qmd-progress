
  !> Routine to compute the Fermi level given a set of eigenvalues and a temperature.
  !! It applies the Bisection method over the function:
  !! \f$ g(\mu) = \sum_k 2 f(\epsilon_k - \mu) - N = 0 \f$
  !! Where \f$ f(\epsilon_k - \mu) = \frac{1}{1+\exp{(\epsilon_k - \mu)/(k_bT)}}\f$.
  !! \param eigenvalues Eigenvalues of the system (\f$ \{ \epsilon_k \} \f$).
  !! \param kbt Temperature times the Boltzmann's constant  (\f$ k_bT  \f$).
  !! \param bndfil Filing factor (\f$ N_{el}/(2*N_{orbs})\f$).
  !! \param tol Tolerance for the bisection method.
  !! \param Ef Fermi level (\f$ \mu \f$).
  !!
  subroutine prg_get_flevel(eigenvalues,kbt,bndfil,tol,Ef)

    integer                  ::  i, j, k, m
    integer                  ::  norb
    real(dp)                 ::  Ft1, Ft2, Kb, Prod
    real(dp)                 ::  T, step, tol, nel
    real(dp), intent(in)     ::  bndfil, eigenvalues(:), kbt
    real(dp), intent(inout)  ::  Ef

    norb = size(eigenvalues,dim=1)
    nel = bndfil*2.0_dp*norb
    Ef=eigenvalues(1)
    step=abs((eigenvalues(norb)-eigenvalues(1)))
    Ft1=0.0_dp
    Ft2=0.0_dp
    prod=0.0_dp

    !Sum of the occupations
    do i=1,norb
      ft1 = ft1 + 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo
    ft1=ft1-nel

    do m=1,1000001

      if(m.gt.1000000)then
        stop "Bisection method in prg_get_flevel not converging ..."
      endif

      if(abs(ft1).lt.tol)then !tolerance control
        return
      endif

      ef = ef + step

      ft2=0.0_dp

      !New sum of the occupations
      do i=1,norb
        ft2 = ft2 + 2.0_dp*fermi(eigenvalues(i),ef,kbt)
      enddo

      ft2=ft2-nel

      !Product to see the change in sign.
      prod = ft2*ft1

      if(prod.lt.0)then
        ef=ef-step
        step=step/2.0_dp !If the root is inside we shorten the step.
      else
        ft1=ft2  !If not, Ef moves forward.
      endif

    enddo

  end subroutine prg_get_flevel

