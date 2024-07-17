!> Library interface
!! \brief This file is used to interface to python via iso_c_binding 
!! library. 
!! \param nats Number of total atoms in the system
!! \param nTypes Number of atom types 
!! \param coords_in Coordinates of every atom in the system.
!! Allocation:
!! \verbatim coordinate(3,nats) \endverbatim
!! \param latticeVectors_in Flattened lattice vectors/box 
!! Allocation:
!! \verbatim  lattice_vector(3*3) \endverbatim
!! \verbatim  v1 = lattice_vector(1:3) \endverbatim
!! \param atomTypes_in Atom type index for every atom
!! It gives the species index of a particulat atom. Indexing starts from 0! 
!! Allocation:
!! \verbatim  atomTypes(nats) \endverbatim
!! If we need the index of atom 30 then:
!! \verbatim  atomTypes(30) \endverbatim
!! \param atomicNumbers_in Atomic number for every species.
!! A list with the atomic numbers for every species.
!! Allocation:
!! \verbatim  atomicNumbers(nTypes) \endverbatim
!! \return forces_out Computed forces. Flattened 2D array.
!! \verbatim forces_out(2) \endverbatim : y component of the force for atom 1.
!! \param charges_out Vector of computed charges.
!! \param verb_in Verbosity level.
!!
function gpmd_fortran(nats,nTypes,coords_in,latticeVectors_in,atomTypes_in,atomicNumbers_in,&
    &field_in,charges_out,forces_out,dipole_out,verb_in) result(err) bind(c, name='gpmd_compute')
  
    use iso_c_binding, only: c_char, c_double, c_int, c_bool
    use gpmdcov_vars
    use gpmdcov_lib_mod

    implicit none
    integer(c_int), intent(in), value  :: nats
    integer(c_int), intent(in), value  :: nTypes
    real(c_double), intent(inout)  :: coords_in(3*nats)
    real(c_double), intent(inout)  :: field_in(3)
    real(c_double), intent(inout)  :: forces_out(3*nats)
    real(c_double), intent(inout)  :: charges_out(nats)
    real(c_double), intent(inout)  :: dipole_out(3)
    !real(c_double), intent(inout)  :: bornch_out(9*nats)
    integer(c_int), intent(inout)  :: atomTypes_in(nats)
    integer(c_int) ,intent(inout) :: atomicNumbers_in(nTypes)
    real(c_double) ,intent(inout) :: latticeVectors_in(9)
    integer(c_int), intent(in), value :: verb_in
    logical(c_bool) :: err

    real(dp), allocatable :: coords(:,:)
    real(dp), allocatable :: forces(:,:), charges(:), dipole(:)
    !real(dp), allocatable :: bornch(:,:)
    real(dp), allocatable :: field(:)
    real(dp), allocatable :: latticeVectors(:,:)
    integer, allocatable :: atomTypes(:), atomicNumbers(:)
    integer :: k
    integer :: verb
   
    err = .true.
    allocate(coords(3,nats))
    allocate(atomTypes(nats))
    allocate(atomicNumbers(nTypes))
    allocate(latticeVectors(3,3))
    allocate(charges(nats))
    allocate(forces(3,nats))
    allocate(dipole(3))
    !allocate(bornch(9,nats))
    allocate(field(3))
    
    !Note that arrays appear in another order. We need to rearange 
    !the data. This is because of the column mayor (in python) vs. 
    !row mayor in fortran. 
    do k = 1, nats
      coords(1,k) = coords_in((k-1)*3 + 1)
      coords(2,k) = coords_in((k-1)*3 + 2)
      coords(3,k) = coords_in((k-1)*3 + 3)
    enddo

    latticeVectors(1,1) = latticeVectors_in(1)
    latticeVectors(1,2) = latticeVectors_in(2)
    latticeVectors(1,3) = latticeVectors_in(3)

    latticeVectors(2,1) = latticeVectors_in(4)
    latticeVectors(2,2) = latticeVectors_in(5)
    latticeVectors(2,3) = latticeVectors_in(6)

    latticeVectors(3,1) = latticeVectors_in(7)
    latticeVectors(3,2) = latticeVectors_in(8)
    latticeVectors(3,3) = latticeVectors_in(9)

    atomicNumbers = atomicNumbers_in

    do k = 1,nats !We correct for the indexing
      atomTypes(k) = atomTypes_in(k) + 1
    enddo

    field = field_in

    verb = verb_in
    call gpmd_compute(coords,atomTypes,atomicNumbers,latticeVectors,&
      &field,charges,forces,dipole,verb)

    err = .false.

    !We vectorize/flatten the forces to send back to python
    do k = 1, nats
      forces_out((k-1)*3 + 1) = forces(1,k)
      forces_out((k-1)*3 + 2) = forces(2,k)
      forces_out((k-1)*3 + 3) = forces(3,k)
    enddo

!    do k = 1, nats
!      bornch_out((k-1)*9 + 1) = bornch(1,k)
!      bornch_out((k-1)*9 + 2) = bornch(2,k)
!      bornch_out((k-1)*9 + 3) = bornch(3,k)
!      
!      bornch_out((k-1)*9 + 4) = bornch(4,k)
!      bornch_out((k-1)*9 + 5) = bornch(5,k)
!      bornch_out((k-1)*9 + 6) = bornch(6,k)
!      
!      bornch_out((k-1)*9 + 7) = bornch(7,k)
!      bornch_out((k-1)*9 + 8) = bornch(8,k)
!      bornch_out((k-1)*9 + 9) = bornch(9,k)
!    enddo

    charges_out(:) = charges(:)
    dipole_out(:) = dipole(:)
   
    deallocate(coords)
    deallocate(forces)
    deallocate(charges)
    deallocate(latticeVectors)
    deallocate(atomTypes)
    deallocate(atomicNumbers)
    deallocate(dipole)
    !deallocate(bornch)

    return 

end function gpmd_fortran



