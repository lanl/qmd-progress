!Small driver to call the gpmd main module used 
!for compiling the code as a library.

program gpmd
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
  use gpmdcov_lib_mod
  ! Progress modules
  use prg_ptable_mod
  use prg_system_mod

  implicit none
  real(dp), allocatable :: coords_in(:,:)
  integer, allocatable :: atomTypes_in(:)
  integer, allocatable :: atomic_numbers_in(:)
  real(dp), allocatable :: lattice_vectors_in(:,:)
  real(dp), allocatable :: forces_out(:,:)
  real(dp), allocatable :: charges_out(:)
  !real(dp), allocatable :: borncharges_out(:,:)
  real(dp), allocatable :: field_in(:)
  real(dp), allocatable :: dipole_out(:)
  type(system_type) :: mysys
  integer :: verb_in

  call prg_parse_system(mysys,"benzene-CN.pdb")
  allocate(coords_in(3,mysys%nats))
  allocate(atomTypes_in(mysys%nsp))
  allocate(atomic_numbers_in(mysys%nsp))
  allocate(lattice_vectors_in(3,3))
  allocate(charges_out(mysys%nats))
  allocate(forces_out(3,mysys%nats))
  coords_in = mysys%coordinate
  atomTypes_in = mysys%spindex
  atomic_numbers_in = mysys%spatnum
  lattice_vectors_in = mysys%lattice_vector
  verb_in = 3
  allocate(field_in(3))
  field_in = 0.0_dp
  call gpmd_compute(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in,&
     &field_in,charges_out,forces_out,dipole_out,verb_in)

  call prg_parse_system(mysys,"benzene-CN.pdb")
  deallocate(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in)
  allocate(coords_in(3,mysys%nats))
  allocate(atomTypes_in(mysys%nsp))
  allocate(atomic_numbers_in(mysys%nsp))
  allocate(lattice_vectors_in(3,3))
  deallocate(charges_out)
  deallocate(forces_out)
  allocate(charges_out(mysys%nats))
  allocate(forces_out(3,mysys%nats))
  coords_in = mysys%coordinate
  atomTypes_in = mysys%spindex
  atomic_numbers_in = mysys%spatnum
  lattice_vectors_in = mysys%lattice_vector
  verb_in = 1
  call gpmd_compute(coords_in,atomTypes_in,atomic_numbers_in,lattice_vectors_in,&
     &field_in,charges_out,forces_out,dipole_out,verb_in)


end program gpmd

