module gpmdcov_lib_mod

  public :: gpmd_compute

contains

  !> High-level program to perform GRAPH-BASED QMD. LIB MODULE
  !!
  !! \ingroup PROGRAMS
  !! \brief  High-level program to perform GRAPH-BASED QMD with a DFTB Hamiltonian
  !!  using a full parallel Graph-based approach.
  !! \note To test this program with the \verbatim run_test \endverbatim script.
  !!
  !! \author C. F. A. Negre
  !! (cnegre@lanl.gov)
  !!
  !! \author S. Mniszewski
  !! (smn@lanl.gov)
  !!
  !! \author A. M. N. Niklasson
  !! (amn@lanl.gov)
  !!
  !! \author M. E. Wall
  !! (mewall@lanl.gov)
  !!
  !! Verbose levels:
  !!   >= 0- Only print tags.
  !!   >= 1- Print useful scalar physical data (e.g., Total Energy)
  !!   >= 2- Print single file output data (e.g., 0-SCF Charges)
  !!   >= 3- Write out trajectory data.
  !!   >= 4- Write out physically meaningful 1D arrays (e.g., Charges for the species)
  !!   >= 5- Write out physically meaningful 2D arrays (e.g., H matrix)
  !!
  !subroutine gpmd_compute(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in,&
  !    &forces_out,charges_out,verb_in)
  !subroutine gpmd_compute(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in,&
  !    &forces_out,charges_out,verb_in)
  !subroutine gpmd_compute(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in,&
  !    &charges_out,forces_out,verb_in)
  subroutine gpmd_compute(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in,&
       &field_in,charges_out,forces_out,dipole_out,energy_out,verb_in)


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
    use gpmdcov_highorder_mod
    use prg_ptable_mod
    use prg_response_mod

    ! Latte modules
#ifdef USE_LATTE
    use gpmdcov_latte_mod
    use setuparray
    use prg_system_mod
    use latteparser
    use constants_mod
    use bml
    use nonoarray
    use NEBLISTARRAY
#endif

    implicit none
    real(dp), allocatable, intent(in) :: coords_in(:,:), lattice_vectors_in(:,:)
    real(dp), allocatable, intent(in) :: field_in(:)
    real(dp), allocatable, intent(inout) :: charges_out(:)
    real(dp), allocatable, intent(inout) :: forces_out(:,:)
    !real(dp), allocatable, intent(inout) :: bornch_out(:,:)
    real(dp), allocatable, intent(inout) :: dipole_out(:)
    integer, allocatable, intent(in) :: atomTypes_in(:)
    real(dp) :: energy_out(1)
    integer, intent(in) :: verb_in
    integer :: k
    integer, allocatable, intent(in) :: atomic_numbers_in(:)

    real(dp),allocatable :: dipoleMoment(:),bornCharges(:,:)
    real(dp),allocatable :: savedDipoleMoment(:),savedCoords(:,:)
    real(dp) :: factor,dr,dfx,dfy,dfz
    integer :: atomI

    !integer :: ipreMD
!!!!!!!!!!!!!!!!!!!!!!!!
    !> Main program driver
!!!!!!!!!!!!!!!!!!!!!!!!

    lib_mode = .true.
    err_status = .false.

    if(verb_in < 0)then
      open(unit=6, file="/dev/null", form="formatted")
    endif

    !> Initialize the program variables and parse input files.
    sy%nats = size(atomTypes_in)
    allocate(sy%symbol(sy%nats))
    allocate(sy%atomic_number(sy%nats))
    allocate(sy%coordinate(3,sy%nats))
    allocate(sy%mass(sy%nats))
    allocate(sy%lattice_vector(3,3))

    sy%coordinate = coords_in
    sy%lattice_vector = lattice_vectors_in

    !In this code symbols, mass, and atomic number are large
    !arrays running through all the atoms, not only atom types.
    !In gpmd, atomTypes are stored in sp%spindex
    allocate(sy%spindex(sy%nats))
    write(*,*)"atom types",atomTypes_in
    do k = 1,sy%nats
      kk = atomTypes_in(k)
      sy%spindex(k) = kk
      sy%mass(k) = element_mass(atomic_numbers_in(kk))
      sy%symbol(k) = element_symbol(atomic_numbers_in(kk))
      sy%atomic_number(k) = atomic_numbers_in(kk)
    enddo

    sy%nsp = size(atomic_numbers_in)
    allocate(sy%splist(sy%nsp))
    allocate(sy%spatnum(sy%nsp))
    allocate(sy%spmass(sy%nsp))

    do k = 1,sy%nsp
      sy%splist(k) = element_symbol(atomic_numbers_in(k))
      sy%spatnum(k) = atomic_numbers_in(k)
      sy%spmass(k) = element_mass(atomic_numbers_in(k))
    enddo
    write(*,*)"spindex",sy%spindex,"splist",sy%splist

    !deallocate(coords_in, attice_vectors_in)
    !deallocate(atomTypes_in)

    write(*,*)"Bef init"
    call gpmdcov_Init(.true.)

    lt%stopat = "gpmdcov_Energ"
    lt%verbose = verb_in


#ifdef USE_LATTE
    call init_latte()
#endif
    write(*,*)"Bef asset"
    call gpmdcov_assert_input(myRank)

    !We give a first guess of the Fermi level.
    Ef =  lt%efermi

    !The inverse of the electronic temperature.
    beta = 1.0_dp/lt%kbt

    !> Initial partition of the system based on the covalency graph.
    call gpmdcov_Part(1)
    if(lt%stopAt == "gpmdcov_Part") stop

    !> Initialize partitions.

    if(lt%method == "DiagEfFull") eig = .false.

    call gpmdcov_InitParts()
    if(err_status)then
            call gpmdcov_Finalize()
            call prg_shutdownParallel 
            call prg_timer_shutdown
            return
    endif

    if(norm2(field_in) < 1.0D-10)then
      if(.not. eig) then
        call gpmdcov_Diagonalize_H0()
        if(lt%MuCalcType == "FromParts" .or. lt%MuCalcType == "Combined")then
          call gpmdcov_muFromParts()
        endif
      endif
      call gpmdcov_FirstCharges(eig)

      !> First SCF loop up to maxscf.
      if(eig)then
#ifdef USE_NVTX
      call nvtxStartRange("gpmdcov_DM_Min",4)
#endif
        call gpmdcov_DM_Min(lt%maxscf,sy%net_charge,.true.)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      else
#ifdef USE_NVTX
      call nvtxStartRange("gpmdcov_DM_Min_Eig",4)
#endif
        call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      endif
#ifdef USE_NVTX
      call nvtxStartRange("gpmdcov_EnergAndForces",1)
#endif
      !> First calculation of energies and forces.
      call gpmdcov_EnergAndForces(sy%net_charge)
      !stop
#ifdef USE_NVTX
      call nvtxEndRange
#endif

    else
#ifdef USE_NVTX
      call nvtxStartRange("gpmdcov_get_field_forces",2)
#endif
      call gpmdcov_get_field_forces(tb%norbi,atomic_numbers_in,field_in)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
    endif

    charges_out(:) = sy%net_charge(:)
    forces_out(:,:) = sy%force(:,:)
    energy_out(1) = Epot

    if(allocated(dipoleMoment)) deallocate(dipoleMoment)
    allocate(dipoleMoment(3))
    factor = 1.0_dp
#ifdef USE_NVTX
      call nvtxStartRange("prg_compute_dipole",3)
#endif
    call prg_compute_dipole(sy%net_charge,sy%coordinate,dipoleMoment,factor,lt%verbose)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
    dipole_out(:) = dipoleMoment(:)

    if(lt%stopAt == "gpmdcov_Energ") then
      lib_init = .true.
      call gpmdcov_Finalize()
      return
    endif

    !> Setup the Molecular Dynamics (MD) calculation.
    call gpmdcov_PrepareMD(gpmdt%temp0)
    if(lt%stopAt == "gpmdcov_PrepareMD") return

    !> Perform the MD simulation.
    call gpmdcov_MDloop()
    if(lt%stopAt == "gpmdcov_MDloop")then
      call gpmdcov_Finalize()
      return
    endif

  end subroutine gpmd_compute

end module gpmdcov_lib_mod
