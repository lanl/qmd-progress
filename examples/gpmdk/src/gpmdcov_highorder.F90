!> Compute higher order electronic properties, such as polarizability.
!!
module gpmdcov_highorder_mod

  ! Local modules
  use gpmdcov_vars
  use gpmdcov_dm_min_mod
  use gpmdcov_energandforces_mod
  use gpmdcov_prepareMD_mod
  use gpmdcov_mdloop_mod
  use gpmdcov_init_mod
  use gpmdcov_part_mod
  use gpmdcov_assert_mod
  use gpmdcov_mod
  use gpmdcov_diagonalize_mod
  use gpmdcov_kernel_mod
  ! Progress modules
  use prg_ptable_mod
  use prg_response_mod


  !  public :: gpmdcov_get_born_charges_v1
  !public :: gpmdcov_get_born_charges_v2
  public :: gpmdcov_field_forces

contains
#ifdef FULL
  !>  Computes the Born charges
  !! \brief Based on the derivative of the dipole moment with respect
  !! to the coordinates.
  !! \param dipoleMoment Dipole moment of the molecule.
  !! \param bornCharges for each atom
  subroutine gpmdcov_get_born_charges_v1(dipoleMoment,bornCharges)

    implicit none
    real(dp),allocatable :: dipoleMoment(:),savedCharges(:),bornCharges(:,:)
    real(dp),allocatable :: savedDipoleMoment(:),savedCoords(:,:)
    real(dp) :: factor,dr,dfx,dfy,dfz
    integer :: atomI

    if(.not. allocated(dipoleMoment)) allocate(dipoleMoment(3))
    allocate(savedDipoleMoment(3))
    factor = 1.0_dp
    call prg_compute_dipole(sy%net_charge,sy%coordinate,dipoleMoment,factor,lt%verbose)
    write(*,*)"Dipole Moment (e * Ang)=",dipoleMoment
    savedDipoleMoment = dipoleMoment

    !Born Charges
    allocate(savedCoords(3,sy%nats))
    allocate(savedCharges(sy%nats))
    if(.not. allocated(bornCharges)) allocate(bornCharges(9,sy%nats))
    savedCoords = sy%coordinate
    savedCharges = sy%net_charge

    dr = 0.000001
    do atomI = 1,sy%nats

      !write(*,*)"BornCharges",atomI,"----------------------------------"
      sy%coordinate(1,atomI) = sy%coordinate(1,atomI) + dr
      call gpmdcov_InitParts()
      call gpmdcov_Diagonalize_H0()
      call gpmdcov_muFromParts()
      call gpmdcov_FirstCharges(eig)
      call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
      call prg_compute_dipole(sy%net_charge,sy%coordinate,dipoleMoment,factor,lt%verbose)

      dfx = (dipoleMoment(1) - savedDipoleMoment(1))/dr
      dfy = (dipoleMoment(2) - savedDipoleMoment(2))/dr
      dfz = (dipoleMoment(3) - savedDipoleMoment(3))/dr
      sy%coordinate(1,atomI) = savedCoords(1,atomI)
      !write(*,*)"BornCharges",dfx,dfy,dfz,savedCharges(atomI)
      bornCharges(1,atomI) = dfx
      bornCharges(2,atomI) = dfy
      bornCharges(3,atomI) = dfz

      sy%coordinate(2,atomI) = sy%coordinate(2,atomI) + dr
      call gpmdcov_InitParts()
      call gpmdcov_Diagonalize_H0()
      call gpmdcov_muFromParts()
      call gpmdcov_FirstCharges(eig)
      call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
      call prg_compute_dipole(sy%net_charge,sy%coordinate,dipoleMoment,factor,lt%verbose)
      dfx = (dipoleMoment(1) - savedDipoleMoment(1))/dr
      dfy = (dipoleMoment(2) - savedDipoleMoment(2))/dr
      dfz = (dipoleMoment(3) - savedDipoleMoment(3))/dr
      sy%coordinate(2,atomI) = savedCoords(2,atomI)
      !write(*,*)"BornCharges",dfx,dfy,dfz,savedCharges(atomI)
      bornCharges(4,atomI) = dfx
      bornCharges(5,atomI) = dfy
      bornCharges(6,atomI) = dfz

      sy%coordinate(3,atomI) = sy%coordinate(3,atomI) + dr
      call gpmdcov_InitParts()
      call gpmdcov_Diagonalize_H0()
      call gpmdcov_muFromParts()
      call gpmdcov_FirstCharges(eig)
      call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
      call prg_compute_dipole(sy%net_charge,sy%coordinate,dipoleMoment,factor,lt%verbose)
      dfx = (dipoleMoment(1) - savedDipoleMoment(1))/dr
      dfy = (dipoleMoment(2) - savedDipoleMoment(2))/dr
      dfz = (dipoleMoment(3) - savedDipoleMoment(3))/dr
      sy%coordinate(3,atomI) = savedCoords(3,atomI)
      !write(*,*)"BornCharges",dfx,dfy,dfz,savedCharges(atomI)
      bornCharges(7,atomI) = dfx
      bornCharges(8,atomI) = dfy
      bornCharges(9,atomI) = dfz
      !write(*,*)"BornCharges-------------------------------------------"

    enddo

    deallocate(savedCoords)
    deallocate(savedDipoleMoment)


  end subroutine gpmdcov_get_born_charges_v1

  !> Born charges
  !! \brief Will compute the Born charges based on the derivatives
  !! of the forces with respect to the field.
  subroutine gpmdcov_get_born_charges_v2(polarizability,bornCharges,norbi)
    implicit none
    real(dp),allocatable :: dipoleMoment(:),savedCharges(:),bornCharges(:,:)
    real(dp),allocatable :: savedDipoleMoment(:),savedCoords(:,:)
    real(dp) :: lambda,intensity,field(3),threshold
    real(dp), allocatable :: polarizability(:), forcesSaved(:,:)
    integer :: atomI, verbose
    integer, allocatable, intent(in) :: norbi(:)
    type(bml_matrix_t) :: prt_bml
    real(dp) :: relperm, keconst, tfact
    real(dp) :: intensity0,dintensity

    lambda = 1.0_dp !2.0_dp/14.3996437701414_dp
    threshold = 0.0_dp
    verbose = 3

    if(.not. allocated(bornCharges)) allocate(bornCharges(9,sy%nats))

    !This will assume only one part
    call gpmdcov_InitParts()
    allocate(forcesSaved(3,sy%nats))

    intensity0 = 0.0000_dp !eV/Ang
    dintensity = 0.001_dp
    intensity = 1.0D-4
    bornCharges = 0.0_dp

    !We will compute dF/dE in x direction
    field = 0.0_dp ; field(1) = -1.0_dp
    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
         ,prt_bml,threshold,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    !call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
    !   ,prt_bml,threshold,sy%spindex,norbi,verbose)
    sy%net_charge = 0.0_dp
    call gpmdcov_InitParts()
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,1.0_dp,threshold) !H = H + V
    call bml_print_matrix("H",syprt(1)%estr%ham0,0,3,0,3)
    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    if(err_status)return
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)
    forcesSaved = sy%force
    field = 0.0_dp ; field(1) = 1.0_dp
    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
         ,prt_bml,threshold,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    sy%net_charge = 0.0_dp
    call gpmdcov_InitParts()
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,1.0_dp,threshold) !H = H + V
    call bml_print_matrix("H",syprt(1)%estr%ham0,0,3,0,3)
    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)
    bornCharges(1,:) = (sy%force(1,:) - forcesSaved(1,:))/(2.0_dp*intensity)
    bornCharges(2,:) = (sy%force(2,:) - forcesSaved(2,:))/(2.0_dp*intensity)
    bornCharges(3,:) = (sy%force(3,:) - forcesSaved(3,:))/(2.0_dp*intensity)

    !We will compute dF/dE in y direction
    field = 0.0_dp ; field(2) = -1.0_dp
    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
         ,prt_bml,threshold,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    call gpmdcov_InitParts()
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,1.0_dp,lt%threshold) !H = H + V
    call bml_print_matrix("H",syprt(1)%estr%ham0,0,3,0,3)
    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)
    forcesSaved = sy%force
    field = 0.0_dp ; field(2) = 1.0_dp
    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
         ,prt_bml,threshold,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    call gpmdcov_InitParts()
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,1.0_dp,lt%threshold) !H = H + V
    call bml_print_matrix("H",syprt(1)%estr%ham0,0,3,0,3)
    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)
    bornCharges(4,:) = (sy%force(1,:) - forcesSaved(1,:))/(2.0_dp*intensity)
    bornCharges(5,:) = (sy%force(2,:) - forcesSaved(2,:))/(2.0_dp*intensity)
    bornCharges(6,:) = (sy%force(3,:) - forcesSaved(3,:))/(2.0_dp*intensity)

    !We will compute dF/dE in z direction
    field = 0.0_dp ; field(3) = -1.0_dp
    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
         ,prt_bml,threshold,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    call gpmdcov_InitParts()
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,1.0_dp,lt%threshold) !H = H + V
    call bml_print_matrix("H",syprt(1)%estr%ham0,0,3,0,3)
    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)
    forcesSaved = sy%force
    field = 0.0_dp ; field(3) = 1.0_dp
    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,intensity,sy%coordinate,lambda&
         ,prt_bml,threshold,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    call gpmdcov_InitParts()
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,1.0_dp,lt%threshold) !H = H + V
    call bml_print_matrix("H",syprt(1)%estr%ham0,0,3,0,3)
    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)
    bornCharges(7,:) = (sy%force(1,:) - forcesSaved(1,:))/(2.0_dp*intensity)
    bornCharges(8,:) = (sy%force(2,:) - forcesSaved(2,:))/(2.0_dp*intensity)
    bornCharges(9,:) = (sy%force(3,:) - forcesSaved(3,:))/(2.0_dp*intensity)

    relperm = 1.0_dp;
    !! The 14.399 factor corresponds to 1/(4*pi*epsilon0) in eV*Ang
    keconst = 14.3996437701414_dp*relperm;

    deallocate(forcesSaved)

  end subroutine gpmdcov_get_born_charges_v2
#endif

  !> Field forces
  subroutine gpmdcov_get_field_forces(norbi,znucs,field_in)
    implicit none
    real(dp),allocatable :: dipoleMoment(:),savedCharges(:)
    real(dp),allocatable :: savedDipoleMoment(:),savedCoords(:,:)
    real(dp) :: lambda,intensity,field(3),threshold
    real(dp), allocatable ::  forcesSaved(:,:)
    real(dp), allocatable, intent(in) :: field_in(:)
    integer :: atomI, verbose, k
    integer, allocatable, intent(in) :: norbi(:), znucs(:)
    type(bml_matrix_t) :: prt_bml, oham
    real(dp) :: relperm, keconst, tfact, kappa
    real(dp) :: intensity0,dintensity,KE

    threshold = 0.0_dp
    verbose = 3

    kappa = 1.0_dp

    write(*,*)field_in
    field = field_in/norm2(field_in)
    intensity = norm2(field_in)
    lambda = kappa*intensity
    write(*,*)(sy%coordinate)

    call bml_zero_matrix("dense",bml_element_real,dp,norb,norb,prt_bml)
    call prg_pert_constant_field(field,1.0_dp,sy%coordinate,1.0_dp&
         ,prt_bml,0.0_dp,sy%spindex,norbi,verbose,syprt(1)%estr%over)
    call bml_print_matrix("Pert",prt_bml,0,6,0,6)
    call bml_add(syprt(1)%estr%ham0,prt_bml,1.0_dp,-lambda,threshold) !H = H + lambda*V

    call gpmdcov_Diagonalize_H0()
    call gpmdcov_muFromParts()
    if(err_status)return
    call gpmdcov_FirstCharges(eig)
    call gpmdcov_DM_Min_Eig(1000,sy%net_charge,.true.,.false.)
    call gpmdcov_EnergAndForces(sy%net_charge)

    do k = 1,sy%nats
      sy%force(:,k) = sy%force(:,k) + field_in(:)*sy%net_charge(k)
    enddo

  end subroutine gpmdcov_get_field_forces

end module gpmdcov_highorder_mod
