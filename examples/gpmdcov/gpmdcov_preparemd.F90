  !>  Preparing for MD
  !!
  subroutine gpmdcov_prepareMD()
    use gpmdcov_vars

    Implicit none

    !> Initialize velocities
    if(.not.allocated(sy%velocity))then
      allocate(sy%velocity(3,sy%nats))
      sy%velocity(1,:) = 0.0_dp; sy%velocity(2,:) = 0.0_dp; sy%velocity(3,:) = 0.0_dp;    ! Initialize velocities
    endif

    !> Kinetic energy in eV (MVV2KE: unit conversion)
    MVV2KE = 166.0538782_dp/1.602176487_dp

    KE2T = 1.0_dp/0.000086173435_dp

  end subroutine gpmdcov_prepareMD

