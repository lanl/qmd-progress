module gpmdcov_field

  public :: gpmdcov_apply_field

  contains    

  !> Apply a constant field perturbation by changing the Coulombic potentials.
  !! Units can be transformed by using the parameter \f$ \lambda \f$.
  !! \param field Direction of the applied field (\f$ \hat{\textbf{E}} \f$).
  !! \param intensity Intensity of the field (\f$ ||\textbf{E}|| \f$)..
  !! \param coordinate Coordinates of the system (\f$ \textbf{r}\f$).
  !! \param lambda Constant to premultiply the perturbation (\f$ \lambda\f$).
  !! \param threshold Threshold value for bml format matrices.
  !! \param verbose Different levels of verbosity.
  !!
  subroutine gpmdcov_apply_field(field,intensity,coordinates,lambda,final_coul_pot_r,verbose)
    use gpmdcov_vars

    integer                                  ::  cont, l, nats
    integer, intent(in)                      ::  verbose
    real(dp)                                 ::  fieldnorm, fieldx, fieldy, fieldz
    real(dp)                                 ::  intensity, lambda, threshold
    real(dp), intent(in)                     ::  coordinates(:,:), field(3)
    real(dp), allocatable, intent(inout)        ::  final_coul_pot_r(:)


    write(*,*)'Applying a constant field perturbation ...'

    fieldx = field(1)
    fieldy = field(2)
    fieldz = field(3)

    fieldnorm=sqrt(fieldx**2+fieldy**2+fieldz**2) !Get the norm

    fieldx=intensity*fieldx/fieldnorm   !Normalize and add intensity
    fieldy=intensity*fieldy/fieldnorm
    fieldz=intensity*fieldz/fieldnorm

    if(.not. allocated(final_coul_pot_r)) allocate(final_coul_pot_r(sy%nats))

    do l=1,nats
       final_coul_pot_r(l)  = final_coul_pot_r(l) + lambda*(fieldx*coordinates(1,l) + &
             & fieldy*coordinates(2,l) + fieldz*coordinates(3,l))
    enddo

  end subroutine gpmdcov_apply_field

end module gpmdcov_field
