!> Apply constraints to MD between two atoms
!!  Constraints could be potentials, Langrange multipliers, etc.
!! \author Rae Corrigan Grove

module gpmdcov_constraints_mod

contains

  !> Apply a switched harmonic potential to constrain two atoms
  !! 
  !! Harmonic to Linear Energy Function
  !! \param atom1r Location of atom 1 {x,y,z}
  !! \param atom2r Location of atom 2 {x,y,z}
  !! \return force_total The collected force for the two atoms being steered will be updated
  !! \return energy_total The total energy of steering the two atoms together - added to full potential energy total
  subroutine gpmdcov_constraint_harmonicToLinear(atom1r, atom2r, force_total, energy_total, verbose)

    use gpmdcov_writeout_mod
    use gpmdcov_vars
    implicit none

    real(dp), intent(in) :: atom1r(:), atom2r(:)
    real(dp), intent(out) :: energy_total
    real(dp), intent(out) :: force_total(:)
    real(dp), allocatable :: coords(:,:)
    real(dp) :: dcoords(3)
    real(dp), allocatable :: Lbox(:)
    real(dp) :: r_separation
    real(dp) :: switch_low, switch_high, dswitch_low, dswitch_high
    real(dp) :: energy_start, energy_end
    real(dp) :: denergy_start, denergy_end, force_radial
    ! real(dp) :: Lx, Ly, Lz
    integer :: k, natms
    integer, optional, intent(in) :: verbose

    allocate(coords(2,3))
    allocate(Lbox(3))

    do k = 1,3
      coords(1,k) = atom1r(k)
      coords(2,k) = atom2r(k)
      Lbox(k) = sy%lattice_vector(k,k)
    enddo
  
    !> Define separation distance between atoms 1 and 2
    !! Convert to fractional coordinates by dividing by box dimensions
    !! Compute fractional distances
    do k = 1,3
      dcoords(k) = modulo(((atom1r(k) - atom2r(k)) + 0.5_dp*Lbox(k)),Lbox(k)) - 0.5_dp * Lbox(k)
    enddo 
    r_separation = norm2(dcoords)

    call gpmdcov_msI("gpmdcov_constraints","Separation distance between &
            & steered atoms, r_separation "//to_string(r_separation),lt%verbose,myRank)
     
    !> Call switching function
    call switch_fermi(r_separation, gpmdt%smdr0, switch_low, switch_high, dswitch_low, dswitch_high)
     
    ! Energy at the start of the switching function (currently harmonic)
    energy_start = 0.5d0 * gpmdt%smdforceconststart * (r_separation - gpmdt%smdr0)**2
    ! Energy at the end of the switching function (currently linear)
    energy_end = gpmdt%smdforceconstend * r_separation
    ! Smoothly switched energy
    energy_total = switch_low * energy_start + switch_high * energy_end
     
    ! Force - negative derivative of energy_start
    ! Derivative of energy_end will be the input force constant since energy_end is linear in R
    denergy_start = gpmdt%smdforceconststart * (r_separation - gpmdt%smdr0)

    denergy_end = gpmdt%smdforceconstend

    ! Total force across switching function, use chain rule
    force_radial = - dswitch_low * energy_start - switch_low * denergy_start - dswitch_high * energy_end &
                 - switch_high * denergy_end

    !> Determine derivatives wrt x, y, and z directions - multiply by force magnitude to get
    !! forces in x, y, and z directions 

    force_total(:) = dcoords(:)/r_separation*force_radial


    deallocate(coords)

  end subroutine gpmdcov_constraint_harmonicToLinear

    !> Apply a switched harmonic potential to constrain two atoms
  !! 
  !! Harmonic to Linear Energy Function
  !! \param atom1r Location of atom 1 {x,y,z}
  !! \param atom2r Location of atom 2 {x,y,z}
  !! \return force_total The collected force for the two atoms being steered will be updated
  !! \return energy_total The total energy of steering the two atoms together - added to full potential energy total
  subroutine gpmdcov_smd_logistic(atom1r, atom2r, force_total, energy_total, verbose)

    use gpmdcov_writeout_mod
    use gpmdcov_vars
    implicit none

    real(dp), intent(in) :: atom1r(:), atom2r(:)
    real(dp), intent(out) :: energy_total
    real(dp), intent(out) :: force_total(:)
    real(dp), allocatable :: coords(:,:)
    real(dp) :: dcoords(3)
    real(dp), allocatable :: Lbox(:)
    real(dp) :: r_separation, r0, r1, lf
    real(dp) :: switch_low, switch_high, dswitch_low, dswitch_high
    real(dp) :: energy_start, energy_end
    real(dp) :: denergy_start, denergy_end, force_radial
    ! real(dp) :: Lx, Ly, Lz
    integer :: k, natms
    integer, optional, intent(in) :: verbose

    allocate(coords(2,3))
    allocate(Lbox(3))

    do k = 1,3
      coords(1,k) = atom1r(k)
      coords(2,k) = atom2r(k)
      Lbox(k) = sy%lattice_vector(k,k)
    enddo
  
    !> Define separation distance between atoms 1 and 2
    !! Convert to fractional coordinates by dividing by box dimensions
    !! Compute fractional distances
    do k = 1,3
      dcoords(k) = modulo(((atom1r(k) - atom2r(k)) + 0.5_dp*Lbox(k)),Lbox(k)) - 0.5_dp * Lbox(k)
    enddo 
    r_separation = norm2(dcoords)

    call gpmdcov_msI("gpmdcov_constraints","Separation distance between &
            & steered atoms, r_separation "//to_string(r_separation),lt%verbose,myRank)
     
    !> Call switching function
    r0 = gpmdt%smdr0
    r1 = minval(Lbox)*0.2_dp
    lf = 3.0_dp
    
    energy_total = exp(-lf*r0)/(lf*(exp(-lf*r1)-exp(-lf*r0)))*log(abs(exp(-lf*r1)&
            &+exp(-lf*r_separation))/abs(exp(-lf*r0)+exp(-lf*r_separation)))*gpmdt%smdforceconststart
    force_radial = -(1.0_dp - 1.0_dp/(1.0_dp+exp(lf*(r_separation-r0))))/(1.+exp(lf*(r_separation-r1)))*gpmdt%smdforceconststart

    force_total(:) = dcoords(:)/r_separation*force_radial

    deallocate(coords)

  end subroutine gpmdcov_smd_logistic

  !> Polynomial switching function
  !!
  !! \param r_separation Separation distance between steered atoms
  !! \param minsep       Minimum separation distance - steering force applied above this distance
  !! \param maxsep       Maximum separation distance - steering force applied below this distance
  !! \return sl          Low side of switching function
  !! \return sh          High side of switching function
  !! \return dsl         Derivative of low side of switching function
  !! \return dsh         Derivative of high side of switching function
  subroutine switch_polynomial_1(r_separation, minsep, maxsep, sl, sh, dsl, dsh)
    use gpmdcov_vars
    implicit none
    real(dp), intent(in) :: r_separation, minsep, maxsep
    real(dp), intent(out) :: sl, sh, dsl, dsh 
    real(dp) :: x, t
     
    !> If the separation distance, r_separation, is within the constrained distance range,
    !! x gets set to a normalized value between 0 and 1
    !! If the separation distance is outside of the constrained distance range, x is set
    !! to either 0 (when the atoms are closer together than the minimum separation distance)
    !! or 1 (when the atoms are further apart than the maximum separation distance)
    if (r_separation .gt. minsep .AND. r_separation .lt. maxsep) then
        x = (r_separation - minsep)/(maxsep - minsep)
    else if (r_separation .lt. minsep) then
        x = 0.0d0
    else if (r_separation .gt. maxsep) then
        x = 1.0d0
    endif

    t = 1.0d0 - x
    sl = 3.0d0*t**2.0d0 - 2.0d0*t**3.0d0
    dsl = (6.0d0*t - 6.0d0*t**2.0d0) * (-1.0d0)/(maxsep - minsep)
    t = x
    sh = 3.0d0*t**2.0d0 - 2.0d0*t**3.0d0
    dsh = (6.0d0*t - 6.0d0*t**2.0d0) * 1.0d0/(maxsep - minsep)

  end subroutine switch_polynomial_1

  !> Fermi switching function
  !!
  !! \param r     Separation distance between steered atoms
  !! \parma r0    Optimal separation distance between steered atoms - center of switching function
  !! \return sl   Low side of switching function
  !! \return sh   High side of switching function
  !! \return dsl  Derivative of low side of switching function
  !! \return dsh  Derivative of high side of switching function
  subroutine switch_fermi(r, r0, sl, sh, dsl, dsh)
    use gpmdcov_vars
    implicit none
    real(dp), intent(in) :: r, r0
    real(dp), intent(out) :: sl, sh, dsl, dsh
    real(dp) :: x

    x = r - r0
    sl = 1.0d0 - exp(x)/(1.0d0 + exp(x))
    dsl = exp(2.0d0*x)/(1.0d0 + exp(x))**2 - exp(x)/(1.0d0+exp(x))
    sh = exp(x)/(1.0d0 + exp(x))
    dsh = -1.0d0 * (exp(2.0d0*x)/(1.0d0+exp(x))**2) + exp(x)/(1.0d0+exp(x))

  end subroutine switch_fermi

end module gpmdcov_constraints_mod
