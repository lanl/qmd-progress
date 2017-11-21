!> A module to read and handle chemical systems.
!! \brief This module will be used to build and handle a molecular system.
!! @ingroup PROGRESS
!!
!! \author C. F. A. Negre
!! (cnegre@lanl.gov)
!!
module prg_system_mod

  use bml
  use prg_openfiles_mod
  use prg_ptable_mod
  use prg_graph_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> Electronic structure type
  type, public :: estruct_type  !< The electronic structure type.

     !> Number of orbitals of the system.
     integer :: norbs

     !> Number of electrons.
     integer :: nel

     !> Hindex.
     integer, allocatable :: hindex(:,:)

     !> SCC-Hamiltonian of the system.
     type(bml_matrix_t)  ::  ham

     !> Hamiltonian of the system.
     type(bml_matrix_t)  ::  ham0

     !> Orthogonalized Hamiltonian.
     type(bml_matrix_t)  ::  oham

     !> Overlap matrix of the system.
     type(bml_matrix_t)  ::  over

     !> Density matrix of the system.
     type(bml_matrix_t)  ::  rho

     !> Orthogonalized density matrix.
     type(bml_matrix_t)  ::  orho

     !> Congruence transformation.
     type(bml_matrix_t)  ::  zmat

     !> Real Coulombic contribution.
     real(dp), allocatable  ::  coul_pot_r(:)

     !> Reciprocal Coulombic contribution.
     real(dp), allocatable  ::  coul_pot_k(:)

     !> Slater Koster force.
     real(dp), allocatable  ::  skforce(:,:)

     !> Pulay force.
     real(dp), allocatable  ::  fpul(:,:)

     !> Nonorthogonal Coulombic force.
     real(dp), allocatable  ::  fscoul(:,:)

     !> Band energy.
     real(dp) ::  eband

  end type estruct_type

  !> System type
  type, public :: system_type  !< The molecular system type.

     !> Number of atoms of the system.
     integer :: nats

     !> Chemical Symbols for every atom of the system.
     !! Symbol can be recovered using ptable module and calling the
     !! following routine:
     !! \verbatim system%symbol(i) = element_symbol(system%atomic_number(i)) \endverbatim
     !! Allocation:
     !! \verbatim symbol(nats) \endverbatim
     character(2), allocatable :: symbol(:)

     !> Atomic number for every atom in the system.
     integer, allocatable :: atomic_number(:)

     !> Coordinates of every atom in the system.
     !! Allocation:
     !! \verbatim coordinate(3,nats)
     real(dp), allocatable :: coordinate(:,:)

     !> Velocities for every atom in the system.
     !! Allocation:
     !! \verbatim velocity(3,nats)
     real(dp), allocatable :: velocity(:,:)

     !> Forces acting on every atom in the system.
     !! Allocation:
     !! \verbatim  force(3,nats)
     real(dp), allocatable :: force(:,:)

     !> Charges of every atom in the system.
     !! Allocation:
     !! \verbatim  net_charge(nats)
     real(dp), allocatable :: net_charge(:)

     !> Mass of every atom in the system.
     !! These can be automatically loaded by using the structures of the ptable mod:
     !!  \verbatim system%mass(i) = mass(mystem%atomic_number(i)) \endverbatim
     !! Allocation:
     !! \verbatim  mass(nats)
     real(dp), allocatable :: mass(:)

     !> Lattice vectors of the system.
     !! Use the prg_vectors_to_parameters and parameters_to_vector
     !! to transform from lattice vector to lattice parameters.
     !! Allocation:
     !! \verbatim  lattice_vector(3,3) \endverbatim
     !! \verbatim  v1 = lattice_vector(1,:) \endverbatim
     !! \verbatim  v2 = lattice_vector(2,:) \endverbatim
     !! \verbatim  v3 = lattice_vector(3,:) \endverbatim
     real(dp), allocatable :: lattice_vector(:,:)

     !> Reciprocal vectors of the system.
     !! Allocation:
     !! \verbatim  recip_vector(3,3) \endverbatim
     !! \verbatim  v1 = recip_vector(1,:) \endverbatim
     !! \verbatim  v2 = recip_vector(2,:) \endverbatim
     !! \verbatim  v3 = recip_vector(3,:) \endverbatim
     real(dp), allocatable :: recip_vector(:,:)

     !> Volume of the system (direct space).
     !! \note use prg_get_recip_vects in coulomb_latte_mod to compute this.
     real(dp) :: volr

     !> Volume of the system (direct space).
     !! \note use prg_get_recip_vects in coulomb_latte_mod to compute this.
     real(dp) :: volk

     !> Number of different species.
     !> Number of species or number of differet antom types (symbols) in the system.
     !! This integer is alwas less or equal than the total number of atoms (nsp <= nats).
     !! This information can also be found in tbparams structure and the following equality holds:
     !! \verbatim system%nsp = tbparams%nsp \endverbatim
     integer :: nsp

     !> Species index.
     !! It gives the species index of a particulat atom.
     !! Allocation:
     !! \verbatim  spindex(nats) \endverbatim
     !! If we need the index of atom 30 then:
     !! \verbatim  system%spindex(30) \endverbatim
     integer, allocatable :: spindex(:)

     !> Species symbol list.
     !! A list with the different species e.g. H, C, N, etc with the order corresponding
     !! to the appearence in system%symbol.
     !! Allocation:
     !! \verbatim  splist(nsp) \endverbatim
     !!
     character(2), allocatable :: splist(:)

     !> Species atomic number list.
     !! A list with the atomic numbers for every species
     !! Allocation:
     !! \verbatim  spatnum(nsp) \endverbatim
     !!
     integer, allocatable :: spatnum(:)

     !> Species mass list.
     !! A list with the atomic mass for every species
     !! Allocation:
     !! \verbatim  spmass(nsp) \endverbatim
     real(dp), allocatable :: spmass(:)

     !> User define field.
     real(dp), allocatable :: userdef(:)

     !> Residue index
     integer, allocatable :: resindex(:)

     !> Electronic structure
     type(estruct_type)   ::  estr

  end type system_type

  public :: prg_parse_system, prg_get_nameandext, prg_make_random_system, prg_write_system, prg_write_trajectory
  public :: prg_get_origin, prg_get_covgraph, prg_get_subsystem, prg_translateandfoldtobox, prg_molpartition, prg_get_partial_atomgraph
  public :: prg_destroy_subsystems, prg_get_covgraph_h, prg_collect_graph_p, prg_merge_graph, prg_merge_graph_adj, prg_adj2bml, prg_graph2bml
  public :: prg_graph2vector, prg_vector2graph, prg_sortadj, prg_get_recip_vects, prg_translatetogeomcandfoldtobox
  public :: prg_write_trajectoryandproperty, prg_get_distancematrix
  public :: prg_get_dihedral, prg_wraparound

contains

  !> Get the name and extension of a file.
  !! \param fullfilename Full filename.
  !! \param filename Filename of the system.
  !! \param extension Extension of the file.
  !!
  subroutine prg_get_nameandext(fullfilename,filename,ext)
    implicit none
    character(30), intent(in)         ::  fullfilename
    character(30), intent(inout)      ::  filename
    character(3), intent(inout)       ::  ext
    character(1), allocatable         ::  tempc(:)
    character(len=30)                 ::  tempcflex
    integer                           ::  lenc

    lenc=len(adjustl(trim(fullfilename)))
    if(.not.allocated(tempc))allocate(tempc(lenc))
    tempcflex = adjustl(trim(fullfilename))
    filename = adjustl(trim(tempcflex(1:lenc-4)))
    ext = adjustl(trim(tempcflex(lenc-2:lenc+1)))

  end subroutine prg_get_nameandext


  !> The parser for the chemical system.
  !! \param system System to be constructed.
  !! \param filename Filename of the system.
  !! \param extin Extension of the file.
  !!
  subroutine prg_parse_system(system,filename,extin)
    implicit none
    character(1)                    ::  onechar
    character(10)                   ::  dummyc(10)
    character(2)                    ::  possibleNewSymbol, twochar
    character(2), allocatable       ::  spTempSymbols(:)
    character(3), optional, intent(in)  ::  extin
    character(3)                    ::  extension
    character(30)                   ::  dummy, io_name, nametmp
    character(60)                   ::  pdbformat
    character(len=*)                ::  filename
    integer                         ::  dummyi(10), header_lines, i, io_unit
    integer                         ::  ios, itemp, j, lines_to_lattice
    integer                         ::  max_lines, nats, nsp, verbose
    integer                         ::  lines_to_masses, lines_to_atom
    logical                         ::  islattice
    real(dp)                        ::  abc_angles(2,3), dummyr(10), lbx1
    real(dp)                        ::  lbx2, lby1, lby2, lbz1
    real(dp)                        ::  lbz2, max_x, max_y, max_z
    real(dp)                        ::  min_x, min_y, min_z, scfactor
    type(system_type), intent(out)  ::  system

    verbose = 1
    islattice = .false.

    max_x = 0.0_dp; max_y = 0.0_dp; max_z = 0.0_dp;
    min_x = 0.0_dp; min_y = 0.0_dp; min_z = 0.0_dp;

    if(allocated(system%symbol)) stop "ERROR: System already allocated"

    if(.not.present(extin))then
       call prg_get_nameandext(filename,nametmp,extension)
       filename = nametmp
    else
       extension = extin
    endif

    select case(extension)

    case("xyz")

       !! For xyz format see http://openbabel.org/wiki/XYZ_%28format%29
       io_name=trim(filename)//".xyz"
       call prg_open_file_to_read(io_unit,io_name)
       read(io_unit,*)nats
       read(io_unit,*)
       system%nats = nats
       allocate(system%symbol(nats))
       allocate(system%atomic_number(nats))
       allocate(system%coordinate(3,nats))
       allocate(system%mass(nats))
       allocate(system%lattice_vector(3,3))

       system%lattice_vector = 0.0_dp

       do i=1,nats
          read(io_unit,*)system%symbol(i),system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i)
          system%atomic_number(i) = element_atomic_number(system%symbol(i))
          system%mass(i) = element_mass(system%atomic_number(i))
          max_x = max(max_x,system%coordinate(1,i))
          min_x = min(min_x,system%coordinate(1,i))
          max_y = max(max_y,system%coordinate(2,i))
          min_y = min(min_y,system%coordinate(2,i))
          max_z = max(max_z,system%coordinate(3,i))
          min_z = min(min_z,system%coordinate(3,i))
          ! write(*,*)system%symbol(i),system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i)
       enddo

       !The following is not part of an xyz format.
       !VMD, babel, xmakemol and pymol can still read this file.
       read(io_unit,*,iostat=ios)dummy
       if(dummy.eq."#lattice")then
          if(verbose.eq.1)write(*,*)"There is a lattice ..."
          read(io_unit,*)system%lattice_vector(1,1),system%lattice_vector(1,2),system%lattice_vector(1,3)
          read(io_unit,*)system%lattice_vector(2,1),system%lattice_vector(2,2),system%lattice_vector(2,3)
          read(io_unit,*)system%lattice_vector(3,1),system%lattice_vector(3,2),system%lattice_vector(3,3)
       else
          !We will create the lattice vectors if they are not pressent.
          !It will add 10.0 Ang to each coordinate.
          system%lattice_vector = 0.0_dp
          system%lattice_vector(1,1) = max_x - min_x + 10.0_dp
          system%lattice_vector(2,2) = max_y - min_y + 10.0_dp
          system%lattice_vector(3,3) = max_z - min_z + 10.0_dp
       endif

       close(io_unit)

    case("pdb")

       !! For PDB format see http://www.wwpdb.org/documentation/file-format
       io_name=trim(filename)//".pdb"
       call prg_open_file_to_read(io_unit,io_name)
       header_lines = 0
       lines_to_lattice = 0
       max_lines = 1000000

       !Counting header lines
       do i=1,max_lines
          read(io_unit,*)dummy
          if(dummy.eq."ATOM".or.dummy.eq."HETATM") then
             exit
          else
             header_lines = header_lines + 1
             if(dummy.eq."CRYST1")then
                lines_to_lattice = header_lines
                islattice = .true.
             endif
          endif
       enddo

       nats = 1
       do i=1,max_lines
          read(io_unit,*)dummy
          if(dummy.eq."ATOM".or.dummy.eq."HETATM") then
             nats = nats + 1
          else
             exit
          endif
       enddo
       close(io_unit)

       call prg_open_file_to_read(io_unit,io_name)
       do i=1,header_lines
          if(i.eq.lines_to_lattice)then
             read(io_unit,*)dummy,abc_angles(1,1),abc_angles(1,2),abc_angles(1,3)&
                  ,abc_angles(2,1),abc_angles(2,2),abc_angles(2,3)
          else
             read(io_unit,*)dummy
          endif
       enddo

       system%nats = nats
       allocate(system%symbol(nats))
       allocate(system%atomic_number(nats))
       allocate(system%coordinate(3,nats))
       allocate(system%mass(nats))
       allocate(system%lattice_vector(3,3))

       system%lattice_vector = 0.0_dp

       pdbformat= '(A4,A2,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'

       do i=1,nats
          read(io_unit,pdbformat)dummyc(1),dummyc(2),dummyi(1), &
               dummyc(4),dummyc(5),dummyc(6),dummyc(7),dummyi(2),dummyc(8),&
               system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
               dummyr(1),dummyr(2),system%symbol(i),dummyc(10)

          ! In case there are no symbols in the last column:
          if(dummyc(4).ne."".and.system%symbol(i).eq."")then
             onechar=adjustl(trim(dummyc(4)))
             if(onechar.ne."H".and. &
                  onechar.ne."B".and. &
                  onechar.ne."C".and. &
                  onechar.ne."N".and. &
                  onechar.ne."O".and. &
                  onechar.ne."F".and. &
                  onechar.ne."P".and. &
                  onechar.ne."S".and. &
                  onechar.ne."I")then
                twochar=adjustl(trim(dummyc(4)))
                system%symbol(i)=twochar
             else
                system%symbol(i)=onechar
             endif
             system%atomic_number(i) = element_atomic_number_upper(system%symbol(i))
             system%symbol(i) = element_symbol(system%atomic_number(i)) !upper to lower char.
          else
             system%atomic_number(i) = element_atomic_number(system%symbol(i))
          endif

          ! Getting the atomic mass
          system%mass(i) = element_mass(system%atomic_number(i))

          ! Getting the system limits.
          max_x = max(max_x,system%coordinate(1,i))
          min_x = min(min_x,system%coordinate(1,i))
          max_y = max(max_y,system%coordinate(2,i))
          min_y = min(min_y,system%coordinate(2,i))
          max_z = max(max_z,system%coordinate(3,i))
          min_z = min(min_z,system%coordinate(3,i))

       enddo

       if(islattice)then
          call prg_parameters_to_vectors(abc_angles,system%lattice_vector)
       else
          !We will create the lattice vectors if they are not pressent
          !It will add 10.0 Ang to each coordinate.
          system%lattice_vector = 0.0_dp
          system%lattice_vector(1,1) = max_x - min_x + 10.0_dp
          system%lattice_vector(2,2) = max_y - min_y + 10.0_dp
          system%lattice_vector(3,3) = max_z - min_z + 10.0_dp
       endif

       close(io_unit)

    case("ltt")

       !! For old inputblock.dat (old dat) format see the LATTE manual.
       io_name=trim(filename)//".ltt"
       call prg_open_file_to_read(io_unit,io_name)
       read(io_unit,*)dummy, nats
       read(io_unit,*)scfactor
       system%nats = nats
       allocate(system%symbol(nats))
       allocate(system%atomic_number(nats))
       allocate(system%coordinate(3,nats))
       allocate(system%mass(nats))
       allocate(system%lattice_vector(3,3))

       !!The gets the llattice vectors from the box boundaries.
       read(io_unit,*)lbx1,lbx2,lby1,lby2,lbz1,lbz2
       system%lattice_vector = 0.0_dp
       system%lattice_vector(1,1) = (lbx2-lbx1)*scfactor
       system%lattice_vector(2,2) = (lby2-lby1)*scfactor
       system%lattice_vector(3,3) = (lbz2-lbz1)*scfactor

       do i=1,nats
          read(io_unit,*)system%symbol(i),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
          system%coordinate(1,i)= system%coordinate(1,i)*scfactor
          system%coordinate(2,i)= system%coordinate(2,i)*scfactor
          system%coordinate(3,i)= system%coordinate(3,i)*scfactor
          system%atomic_number(i) = element_atomic_number(system%symbol(i))
          system%mass(i) = element_mass(system%atomic_number(i))
          ! write(*,*)system%symbol(i),system%coordinate(1,i)&
          ! ,system%coordinate(2,i),system%coordinate(3,i)
       enddo

       close(io_unit)

    case("dat")

       !! For new inputblock.dat (dat) format see the LATTE manual.
       io_name=trim(filename)//".dat"
       call prg_open_file_to_read(io_unit,io_name)
       read(io_unit,*)nats
       system%nats = nats
       allocate(system%symbol(nats))
       allocate(system%atomic_number(nats))
       allocate(system%coordinate(3,nats))
       allocate(system%mass(nats))
       allocate(system%lattice_vector(3,3))

       !!The gets the llattice vectors from the box boundaries.
       system%lattice_vector = 0.0_dp
       read(io_unit,*)system%lattice_vector(1,1),system%lattice_vector(1,2),system%lattice_vector(1,3)
       read(io_unit,*)system%lattice_vector(2,1),system%lattice_vector(2,2),system%lattice_vector(2,3)
       read(io_unit,*)system%lattice_vector(3,1),system%lattice_vector(3,2),system%lattice_vector(3,3)

       do i=1,nats
          read(io_unit,*)system%symbol(i),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
          system%coordinate(1,i)= system%coordinate(1,i)
          system%coordinate(2,i)= system%coordinate(2,i)
          system%coordinate(3,i)= system%coordinate(3,i)
          system%atomic_number(i) = element_atomic_number(system%symbol(i))
          system%mass(i) = element_mass(system%atomic_number(i))
       enddo

       close(io_unit)

    case("gen")

       stop "gen format is not implemented for reading yet"

    case("lmp")

       !! For Lammps data.* input file.
       !! For more information see: http://lammps.sandia.gov/doc/2001/data_format.html
       io_name=adjustl(trim(filename))//".lmp"
       call prg_open_file(io_unit,io_name)
       read(io_unit,*)dummyc(1)
       read(io_unit,*)system%nats

       do i=1,100
          read(io_unit,*)dummyc(1),dummyc(2)!,dummyc(3)

          if(adjustl(trim(dummyc(2))) == "atom")then
             lines_to_atom = i + 2
          endif
          if(adjustl(trim(dummyc(1))) == "Masses")then
             lines_to_masses = i + 2
             exit
          endif
       enddo
       close(io_unit)

       call prg_open_file(io_unit,io_name)

       do i=1,lines_to_atom-1
          read(io_unit,*)dummyc(1)
       enddo

       read(io_unit,*)system%nsp

       do i=1,lines_to_masses - lines_to_atom - 4
          read(io_unit,*)dummyc(1)
       enddo

       read(io_unit,*)min_x, max_x
       read(io_unit,*)min_y, max_y
       read(io_unit,*)min_z, max_z

       allocate(system%lattice_vector(3,3))
       system%lattice_vector = 0.0_dp
       system%lattice_vector(1,1) = max_x - min_x
       system%lattice_vector(2,2) = max_y - min_y
       system%lattice_vector(3,3) = max_z - min_z

       read(io_unit,*)dummyc(1)

       allocate(system%spmass(system%nsp))
       allocate(system%splist(system%nsp))

       do i=1,system%nsp
          read(io_unit,*)dummyi(1),system%spmass(i)
       enddo

       do i=1,system%nsp
          do j=1,size(element_mass,dim=1)
             if(abs(system%spmass(i) - element_mass(j)) < 1.0) then
                system%splist(i) = element_symbol(j)
                exit
             endif
          enddo
       enddo

       do i=1,100000
          read(io_unit,*)dummyc(1)
          if(dummyc(1) == "Atoms")exit
       enddo

       allocate(system%spindex(system%nats))
       allocate(system%coordinate(3,system%nats))
       allocate(system%symbol(system%nats))

       do i=1,system%nats
          read(io_unit,*)dummyi(1),dummyi(2),system%spindex(i),dummyr(1),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
          write(*,*)dummyi(1),dummyi(2),system%spindex(i),dummyr(1),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
          system%symbol(i) = system%splist(system%spindex(i))
       enddo

       close(io_unit)

    case default

       stop "The file extension is not valid. Only xyz, lmp, dat and pdb formats are implemented"

    end select

    !!Computing species and species index
    !!Assigning species index to the system elements.
    nsp=1
    allocate(spTempSymbols(150)) !Temporary vector to accumulate the species.
    spTempSymbols(1)=system%Symbol(1)

    if(.not.allocated(system%spindex)) allocate(system%spindex(nats))

    do i=1,system%nats
       itemp = 0
       do j=1,nsp !If there is a new chemical symbol
          if(trim(system%symbol(i)).ne.trim(spTempSymbols(j)))then
             itemp = itemp + 1
             possibleNewSymbol = system%symbol(i)
          endif
       enddo
       if(itemp == nsp)then
          nsp = nsp+1
          spTempSymbols(nsp)=possibleNewSymbol
       endif
    enddo

    if(allocated(system%splist))deallocate(system%splist)
    if(allocated(system%spmass))deallocate(system%spmass)

    allocate(system%splist(nsp))
    allocate(system%spmass(nsp))
    allocate(system%spatnum(nsp))

    do i=1,nsp
       system%splist(i) = spTempSymbols(i)
       system%spatnum(i) = element_atomic_number(system%splist(i))
       system%spmass(i) = element_mass(system%spatnum(i))
    enddo

    deallocate(spTempSymbols)
    system%nsp = nsp

    !> Assignment of species index for every atom.
    !! \todo Integrate this loop in the loop for building the splist.
    do i = 1,system%nats
       do j=1,nsp
          if(trim(system%symbol(i)).eq.trim(system%splist(j)))then
             system%spindex(i)=j
          endif
       enddo
    enddo

  end subroutine prg_parse_system

  !>  Write system in .xyz, .dat or pdb file.
  !! \param system System to be constructed.
  !! \param filename File name.
  !! \param extension Extension of the file.
  !!
  subroutine prg_write_system(system,filename,extension)
    implicit none
    character(*)                   ::  filename
    character(10)                  ::  dummyc(10)
    character(11)                  ::  xyzformat
    character(3)                   ::  extension
    character(30)                  ::  io_name
    character(60)                  ::  pdbformat
    integer                        ::  dummyi(10), i, io_unit, nats
    integer, allocatable           ::  resindex(:)
    real(dp)                       ::  abc_angles(2,3), dummyr(10), origin(3)
    type(system_type), intent(in)  ::  system

    nats = system%nats

    select case(extension)

    case("xyz")

       io_name=trim(filename)//".xyz"
       call prg_open_file(io_unit,io_name)
       write(io_unit,*)nats
       write(io_unit,*)io_name, "Generated by the PROGRESS library"
       xyzformat = '(A2,3F10.5)'
       do i=1,nats
          write(io_unit,xyzformat)system%symbol(i),system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i)
       enddo

       !The following is not part of an xyz format but
       !VMD, babel, xmakemol and pymol can still read this file.
       write(io_unit,*)"#lattice vectors"
       write(io_unit,"(3F10.5)")system%lattice_vector(1,1),system%lattice_vector(1,2),system%lattice_vector(1,3)
       write(io_unit,"(3F10.5)")system%lattice_vector(2,1),system%lattice_vector(2,2),system%lattice_vector(2,3)
       write(io_unit,"(3F10.5)")system%lattice_vector(3,1),system%lattice_vector(3,2),system%lattice_vector(3,3)

       close(io_unit)

    case("pdb")

       call prg_vectors_to_parameters(system%lattice_vector,abc_angles)

       dummyi = 0
       dummyc = ""
       dummyc(1) = "ATOM"
       dummyr = 0.0_dp

       io_name=trim(filename)//".pdb"
       call prg_open_file(io_unit,io_name)

       write(io_unit,'(A6,1X,A50)')"REMARK","Generated by PROGRESS library"
       write(io_unit,'(A5,1X,A20)')"TITLE",io_name
       write(io_unit,'(A6,3F9.3,3F7.2,A16)')"CRYST1",abc_angles(1,1),abc_angles(1,2),abc_angles(1,3)&
            ,abc_angles(2,1),abc_angles(2,2),abc_angles(2,3)," P 1           1"
       write(io_unit,'(A5,A20)')"MODEL","        1"

       pdbformat= '(A4,A2,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2 )'

       if(.not.allocated(system%resindex))then
          allocate(resindex(nats))
          resindex = 1
       else
          resindex = system%resindex
       endif

       if(allocated(system%net_charge))then
          do i=1,nats
             write(io_unit,pdbformat)dummyc(1),dummyc(2),i, &
                  system%symbol(i),dummyc(5),"MOL",dummyc(7),resindex(i),dummyc(8),&
                  system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                  dummyr(1),system%net_charge(i),system%symbol(i),dummyc(10)
          enddo
       else
          do i=1,nats
             write(io_unit,pdbformat)dummyc(1),dummyc(2),i, &
                  system%symbol(i),dummyc(5),"MOL",dummyc(7),resindex(i),dummyc(8),&
                  system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                  dummyr(1),dummyr(2),system%symbol(i),dummyc(10)
          enddo
       endif

       write(io_unit,'(A3)')"TER"
       write(io_unit,'(A3)')"ENDMDL"
       close(io_unit)

    case("ltt")

       !! For old inputblock.dat (old dat) format see the LATTE manual.
       io_name=trim(filename)//".ltt"
       call prg_open_file(io_unit,io_name)
       write(io_unit,*)"NATS=  ", system%nats

       !Here we will always write the system in 1:1 scale
       write(io_unit,*)"1.0","  Generated by the PROGRESS library"

       !! This gets the lattice vectors from the box boundaries.
       write(io_unit,"(6F10.5)")0.0,system%lattice_vector(1,1),0.0,system%lattice_vector(2,2)&
            ,0.0,system%lattice_vector(3,3)

       do i=1,nats
          write(io_unit,"(A2,3F10.5)")system%symbol(i),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
       enddo

       close(io_unit)

    case("dat")

       !! For new inputblock.dat (nat) format see the LATTE manual.
       io_name=trim(filename)//".dat"
       call prg_open_file(io_unit,io_name)
       write(io_unit,*)system%nats

       !! This gets the lattice vectors from the box boundaries.
       write(io_unit,"(3F10.5)")system%lattice_vector(1,1),system%lattice_vector(1,2),system%lattice_vector(1,3)
       write(io_unit,"(3F10.5)")system%lattice_vector(2,1),system%lattice_vector(2,2),system%lattice_vector(2,3)
       write(io_unit,"(3F10.5)")system%lattice_vector(3,1),system%lattice_vector(3,2),system%lattice_vector(3,3)

       do i=1,nats
          write(io_unit,"(A2,3F10.5)")system%symbol(i),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
       enddo

       close(io_unit)

    case("gen")

       !! For gen format see DFTB+ manual.
       io_name=trim(filename)//".gen"
       call prg_open_file(io_unit,io_name)
       write(io_unit,*)system%nats,"S"
       write(io_unit,'(103A2)')(system%splist(i),i=1,system%nsp)
       do i=1,nats
          write(io_unit,"(I10,I10,3F10.5)")i,system%spindex(i),system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
       enddo
       write(io_unit,*)"0.0 0.0 0.0"
       write(io_unit,'(3F10.5)')system%lattice_vector(1,1),system%lattice_vector(1,2),system%lattice_vector(1,3)
       write(io_unit,'(3F10.5)')system%lattice_vector(2,1),system%lattice_vector(2,2),system%lattice_vector(2,3)
       write(io_unit,'(3F10.5)')system%lattice_vector(3,1),system%lattice_vector(3,2),system%lattice_vector(3,3)

       close(io_unit)

    case("lmp")

       !! For Lammps data.* input file.
       !! For more information see: http://lammps.sandia.gov/doc/2001/data_format.html
       io_name= adjustl(trim(filename))//".lmp"
       call prg_open_file(io_unit,io_name)
       write(io_unit,*)"LAMMPS Description"
       write(io_unit,*)""
       write(io_unit,*)system%nats,"atoms"
       write(io_unit,*)""
       write(io_unit,*)system%nsp,"atom types"
       write(io_unit,*)""
       write(io_unit,*)0.0_dp,system%lattice_vector(1,1),"xlo xhi"
       write(io_unit,*)0.0_dp,system%lattice_vector(2,2),"ylo yhi"
       write(io_unit,*)0.0_dp,system%lattice_vector(3,3),"zlo zhi"
       write(io_unit,*)""
       write(io_unit,*)"Masses"
       write(io_unit,*)""
       do i=1,system%nsp
          write(io_unit,*)"  ",i,system%spmass(i)
       enddo
       write(io_unit,*)""
       write(io_unit,*)"Atoms"
       write(io_unit,*)""
       do i=1,nats
          write(io_unit,'(3I5,A6,3F10.5)')i,1,system%spindex(i),"   0.0",system%coordinate(1,i)&
               ,system%coordinate(2,i),system%coordinate(3,i)
       enddo

       close(io_unit)

    case default

       stop "The file extension is not valid. Only xyz, dat, lmp, pdb or gen formats are implemented"

    end select

  end subroutine prg_write_system

  !>  Write trajectory in .xyz, .dat or pdb file.
  !! \param system System to be appended to the trajectory file.
  !! \param iter Simulation step.
  !! \param each Writing frequency.
  !! \param filename File name for the trajectory.
  !! \param extension Extension of the file.
  !!
  subroutine prg_write_trajectory(system,iter,each,prg_deltat,filename,extension)
    implicit none
    character(*)                   ::  filename
    character(10)                  ::  dummyc(10)
    character(11)                  ::  xyzformat
    character(20)                  ::  io_name
    character(3)                   ::  extension
    character(60)                  ::  pdbformat
    integer                        ::  dummyi(10), i, io_unit, nats
    integer, intent(in)            ::  iter, each
    integer, allocatable           ::  resindex(:)
    real(dp)                       ::  abc_angles(2,3), dummyr(10)
    real(dp), intent(in)           ::  prg_deltat
    type(system_type), intent(in)  ::  system

    if(mod(iter,each).eq.0.or.iter.eq.0.or.iter.eq.1)then
       nats = system%nats

       select case(extension)

       case("xyz")

          io_name=trim(filename)//".xyz"
          io_unit=get_file_unit(-1)

          if(iter.eq.0.or.iter.eq.1)then
             open(unit=io_unit,file=io_name,Status='unknown')
          else
             open(unit=io_unit,file=io_name,Access = 'append',Status='old')
          endif

          write(io_unit,*)system%nats
          write(io_unit,*)"frame", iter
          if(allocated(system%net_charge))then
             do i=1,nats
                write(io_unit,*)system%symbol(i),system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                     system%net_charge(i)
             enddo
          else
             do i=1,nats
                write(io_unit,*)system%symbol(i),system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i)
             enddo
          endif

          close(io_unit)

       case("pdb")

          call prg_vectors_to_parameters(system%lattice_vector,abc_angles)

          dummyi = 0
          dummyc = ""
          dummyc(1) = "ATOM"
          dummyr = 0.0_dp

          io_name=trim(filename)//".pdb"
          io_unit=get_file_unit(-1)

          if(iter.eq.0.or.iter.eq.1)then
             open(unit=io_unit,file=io_name,Status='unknown')
          else
             open(unit=io_unit,file=io_name,Access = 'append',Status='old')
          endif

          write(io_unit,'(A8,A4,F10.5)')"Trajectory","  t=",iter*prg_deltat
          write(io_unit,'(A24)')"THIS IS A SIMULATION BOX"
          write(io_unit,'(A6,3F9.3,3F7.2,A16)')"CRYST1",abc_angles(1,1),abc_angles(1,2),abc_angles(1,3)&
               ,abc_angles(2,1),abc_angles(2,2),abc_angles(2,3)," P 1           1"
          write(io_unit,'(A5,I6)')"MODEL",iter

          pdbformat= '(A4,A2,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2 )'

          if(.not.allocated(system%resindex))then
             allocate(resindex(nats))
             resindex = 1
          else
             resindex = system%resindex
          endif


          if(allocated(system%net_charge))then
             do i=1,nats
                write(io_unit,pdbformat)dummyc(1),dummyc(2),i, &
                     system%symbol(i),dummyc(5),"MOL",dummyc(7),resindex(i),dummyc(8),&
                     system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                     dummyr(1),system%net_charge(i),system%symbol(i),dummyc(10)
             enddo
          else
             do i=1,nats
                write(io_unit,pdbformat)dummyc(1),dummyc(2),i, &
                     system%symbol(i),dummyc(5),"MOL",dummyc(7),resindex(i),dummyc(8),&
                     system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                     dummyr(1),dummyr(2),system%symbol(i),dummyc(10)
             enddo
          endif

          write(io_unit,'(A3)')"TER"
          write(io_unit,'(A3)')"ENDMDL"
          close(io_unit)

       case("dat")

          write(*,*)"write traj not implemented for dat "

       case default

          stop "The file extension is not valid. Only pdb and xyz formats are implemented"

       end select

    endif

  end subroutine prg_write_trajectory

  !>  Write trajectory and atomic properties. Only pdb file.
  !! \param system System to be appended to the trajectory file.
  !! \param iter Simulation step.
  !! \param each Writing frequency.
  !! \param prg_deltat Integration step.
  !! \param scalarprop Scalar property to plot on atoms.
  !! \param filename File name for the trajectory.
  !! \param extension Extension of the file.
  !!
  subroutine prg_write_trajectoryandproperty(system,iter,each,prg_deltat,scalarprop,filename,extension)
    implicit none
    character(*)                   ::  filename
    character(10)                  ::  dummyc(10)
    character(11)                  ::  xyzformat
    character(20)                  ::  io_name
    character(3)                   ::  extension
    character(60)                  ::  pdbformat
    integer                        ::  dummyi(10), i, io_unit, nats
    integer, intent(in)            ::  iter, each
    integer, allocatable           ::  resindex(:)
    real(dp)                       ::  abc_angles(2,3), dummyr(10)
    real(dp)                       ::  maxprop, minprop, realtmp
    real(dp), intent(in)           ::  prg_deltat, scalarprop(:)
    type(system_type), intent(in)  ::  system


    ! maxprop = maxval(scalarprop(:))
    ! minprop = minval(scalarprop(:))

    if(mod(iter,each).eq.0.or.iter.eq.0.or.iter.eq.1)then
       nats = system%nats

       select case(extension)

       case("pdb")

          call prg_vectors_to_parameters(system%lattice_vector,abc_angles)

          dummyi = 0
          dummyc = ""
          dummyc(1) = "ATOM"
          dummyr = 0.0_dp

          io_name=trim(filename)//".pdb"
          io_unit=get_file_unit(-1)

          if(iter.eq.0.or.iter.eq.1)then
             open(unit=io_unit,file=io_name,Status='unknown')
          else
             open(unit=io_unit,file=io_name,Access = 'append',Status='old')
          endif

          write(io_unit,'(A8,A4,F10.5)')"Trajectory","  t=",iter*prg_deltat
          write(io_unit,'(A24)')"THIS IS A SIMULATION BOX"
          write(io_unit,'(A6,3F9.3,3F7.2,A16)')"CRYST1",abc_angles(1,1),abc_angles(1,2),abc_angles(1,3)&
               ,abc_angles(2,1),abc_angles(2,2),abc_angles(2,3)," P 1           1"
          write(io_unit,'(A5,I6)')"MODEL",iter

          pdbformat= '(A4,A2,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2 )'

          if(.not.allocated(system%resindex))then
             allocate(resindex(nats))
             resindex = 1
          else
             resindex = system%resindex
          endif

          do i=1,nats
             ! realtmp = (scalarprop(i) - minprop)/(maxprop-minprop) - 1.0_dp !Normalization
             write(io_unit,pdbformat)dummyc(1),dummyc(2),i, &
                  system%symbol(i),dummyc(5),"MOL",dummyc(7),resindex(i),dummyc(8),&
                  system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                  dummyr(1),scalarprop(i),system%symbol(i),dummyc(10)
          enddo

          write(io_unit,'(A3)')"TER"
          write(io_unit,'(A3)')"ENDMDL"
          close(io_unit)

       case("xyz")

          io_name=trim(filename)//".xyz"
          io_unit=get_file_unit(-1)

          if(iter.eq.0.or.iter.eq.1)then
             open(unit=io_unit,file=io_name,Status='unknown')
          else
             open(unit=io_unit,file=io_name,Access = 'append',Status='old')
          endif

          write(io_unit,*)system%nats
          write(io_unit,*)"frame", iter
          do i=1,nats
             ! realtmp = (scalarprop(i) - minprop)/(maxprop-minprop) - 1.0_dp !Normalization
             write(io_unit,*)system%symbol(i),system%coordinate(1,i),system%coordinate(2,i),system%coordinate(3,i),&
                  scalarprop(i)
             ! write(*,*)scalarprop(i)
          enddo

          close(io_unit)

       case default

          stop "The file extension is not valid. Only pdb format is implemented"

       end select

    endif

  end subroutine prg_write_trajectoryandproperty


  !>  Make random Xx system.
  !! \param system System to be construucted.
  !! \param nats Number of atoms.
  !! \param lx length of the box for the x coordinate.
  !! \param ly length of the box for the y coordinate.
  !! \param lz length of the box for the z coordinate.
  !!
  subroutine prg_make_random_system(system,nats,seed,lx,ly,lz)

    implicit none
    integer                         ::  i, nats, seed, seed1(12)
    real(dp)                        ::  lx, ly, lz, ran
    type(system_type), intent(out)  ::  system

    !Allocating the system
    system%nats = nats
    allocate(system%symbol(nats))
    allocate(system%atomic_number(nats))
    allocate(system%coordinate(3,nats))
    allocate(system%mass(nats))

    seed1(1) = seed
    call random_seed(put=seed1)

    do i=1,nats

       system%atomic_number(i) = 1
       call random_number(ran)
       system%coordinate(i,1) = ran*lx
       call random_number(ran)
       system%coordinate(i,2) = ran*ly
       call random_number(ran)
       system%coordinate(i,3) = ran*lz

       system%symbol(i) = element_symbol(system%atomic_number(i))
       system%mass(i) = element_mass(1)

    enddo

    allocate(system%lattice_vector(3,3))
    system%lattice_vector = 0.0_dp
    system%lattice_vector(1,1) = lx
    system%lattice_vector(2,2) = ly
    system%lattice_vector(3,3) = lz

  end subroutine prg_make_random_system

  !> Transforms the lattice parameters into lattice vectors.
  !! \param abc_angles 2x3 array containing the lattice parameters.
  !! abc_angles(1,1) = a, abc_angles(1,2) = b, and abc_angles(1,3) = c
  !! abc_angles(2,1) = \f$ \alpha \f$ , abc_angles(2,2) = \f$ \beta \f$ and abc_angles(2,3) = \f$ \gamma \f$
  !! \param lattice_vector 3x3 array containing the lattice vectors.
  !! lattice_vector(1,:) = \f$ \overrightarrow{a} \f$
  !!
  subroutine prg_parameters_to_vectors(abc_angles,lattice_vector)
    implicit none
    real(dp)               ::  angle_alpha, angle_beta, angle_gamma, lattice_a
    real(dp)               ::  lattice_b, lattice_c, pi
    real(dp), intent(in)   ::  abc_angles(2,3)
    real(dp), intent(out)  ::  lattice_vector(3,3)

    pi = 3.14159265358979323846264338327950_dp

    lattice_a = abc_angles(1,1)
    lattice_b = abc_angles(1,2)
    lattice_c = abc_angles(1,3)
    angle_beta = abc_angles(2,1)
    angle_alpha = abc_angles(2,2)
    angle_gamma = abc_angles(2,3)

    angle_alpha = 2.0_dp*pi*angle_alpha/360.0_dp
    angle_beta = 2.0_dp*pi*angle_beta/360.0_dp
    angle_gamma = 2.0_dp*pi*angle_gamma/360.0_dp

    lattice_vector(1,1) = lattice_a
    lattice_vector(1,2)=0
    lattice_vector(1,3)=0

    lattice_vector(2,1)=lattice_b*cos(angle_gamma)
    lattice_vector(2,2)=lattice_b*sin(angle_gamma)
    lattice_vector(2,3)=0

    lattice_vector(3,1)=lattice_c*cos(angle_beta)
    lattice_vector(3,2)=lattice_c*( cos(angle_alpha)-cos(angle_gamma)* &
         cos(angle_beta) )/sin(angle_gamma)
    lattice_vector(3,3)=sqrt(lattice_c**2 - lattice_vector(3,2)**2 - lattice_vector(3,3)**2)

  end subroutine prg_parameters_to_vectors

  !> Transforms the lattice vectors into lattice parameters.
  !! \param lattice_vector 3x3 array containing the lattice vectors.
  !! lattice_vector(1,:) = \f$ \overrightarrow{a} \f$
  !! \param abc_angles 2x3 array containing the lattice parameters.
  !! abc_angles(1,1) = a, abc_angles(1,2) = b and abc_angles(1,3) = c
  !! abc_angles(2,1) = \f$ \alpha \f$, abc_angles(2,2) = \f$ \beta \f$, and abc_angles(2,3) = \f$ \gamma \f$.
  !!
  subroutine prg_vectors_to_parameters(lattice_vector,abc_angles)
    implicit none
    real(dp)               ::  angle_alpha, angle_beta, angle_gamma, lattice_a
    real(dp)               ::  lattice_b, lattice_c, pi
    real(dp), intent(in)   ::  lattice_vector(3,3)
    real(dp), intent(out)  ::  abc_angles(2,3)

    pi = 3.14159265358979323846264338327950_dp

    lattice_a = sqrt(lattice_vector(1,1)**2+lattice_vector(1,2)**2+lattice_vector(1,3)**2)
    lattice_b = sqrt(lattice_vector(2,1)**2+lattice_vector(2,2)**2+lattice_vector(2,3)**2)
    lattice_c = sqrt(lattice_vector(3,1)**2+lattice_vector(3,2)**2+lattice_vector(3,3)**2)

    angle_gamma = dot_product(lattice_vector(1,:),lattice_vector(2,:))/(lattice_a*lattice_b)
    angle_beta = dot_product(lattice_vector(1,:),lattice_vector(3,:))/(lattice_a*lattice_c)
    angle_alpha = dot_product(lattice_vector(2,:),lattice_vector(3,:))/(lattice_b*lattice_c)


    angle_alpha = 360.0_dp*acos(angle_alpha)/(2.0_dp*pi)
    angle_beta = 360.0_dp*acos(angle_beta)/(2.0_dp*pi)
    angle_gamma = 360.0_dp*acos(angle_gamma)/(2.0_dp*pi)

    abc_angles(1,1) = lattice_a
    abc_angles(1,2) = lattice_b
    abc_angles(1,3) = lattice_c

    abc_angles(2,1) = angle_alpha
    abc_angles(2,2) = angle_beta
    abc_angles(2,3) = angle_gamma

  end subroutine prg_vectors_to_parameters

  !> Get the origin of the coordinates.
  !! \param coords Coordinates of teh system (see system_type).
  !! \param origin (min(x),min(y),min(z)) set as the origin of the system.
  !!
  subroutine prg_get_origin(coords,origin)
    implicit none
    integer                              ::  i
    real(dp)                             ::  max_x, max_y, max_z, min_x
    real(dp)                             ::  min_y, min_z
    real(dp), allocatable, intent(inout)  ::  origin(:)
    real(dp), intent(in)                 ::  coords(:,:)

    if(.not.allocated(origin)) allocate(origin(3))

    max_x = -1.0d5 ; max_y = -1.0d5 ; max_z = -1.0d5 ;
    min_x =  1.0d5 ; min_y =  1.0d5 ; min_z =  1.0d5 ;

    write(*,*)size(coords,dim=2)

    ! Getting the system limits.
    do i=1,size(coords,dim=2)
       max_x = max(max_x,coords(1,i))
       write(*,*)coords(1,i)
       min_x = min(min_x,coords(1,i))
       max_y = max(max_y,coords(2,i))
       min_y = min(min_y,coords(2,i))
       max_z = max(max_z,coords(3,i))
       min_z = min(min_z,coords(3,i))
    enddo

    origin(1) = min_x
    origin(2) = min_y
    origin(3) = min_z

  end subroutine prg_get_origin

  !> Get the distance matrix.
  !! \param coords Coordinates of the system (see system_type).
  !! \param dmat Distance matrix (nats x nats).
  !!
  subroutine prg_get_distancematrix(coords,dmat)
    implicit none
    integer                              ::  i,j,nats
    real(dp), intent(in)                 ::  coords(:,:)
    real(dp), intent(out), allocatable   ::  dmat(:,:)

    nats = size(coords,dim=2)
    if(.not.allocated(dmat)) allocate(dmat(nats,nats))

    ! Getting the system limits.
    do i=1,nats
      do j=1,nats
        dmat(i,j) = norm2(coords(:,i)-coords(:,j))
      enddo
    enddo

  end subroutine prg_get_distancematrix


  !> Translate and fold to box.
  !! \param coords Coordinates of the system (see system_type).
  !! \param lattice_vectors System lattice vectors.
  !! \param origin (min(x),min(y),min(z)) set as the origin of the system.
  !!
  subroutine prg_translateandfoldtobox(coords,lattice_vectors,origin)
    implicit none
    integer                              ::  i
    real(dp)                             ::  max_x, max_y, max_z, min_x
    real(dp)                             ::  min_y, min_z
    real(dp), allocatable, intent(inout) ::  origin(:),coords(:,:)
    real(dp), intent(in)                 ::  lattice_vectors(:,:)

    if(.not.allocated(origin)) allocate(origin(3))

    max_x = -1.0d5 ; max_y = -1.0d5 ; max_z = -1.0d5 ;
    min_x =  1.0d5 ; min_y =  1.0d5 ; min_z =  1.0d5 ;

    write(*,*)size(coords,dim=2)

    ! Getting the system limits.
    do i=1,size(coords,dim=2)
       max_x = max(max_x,coords(1,i))
       write(*,*)coords(1,i)
       min_x = min(min_x,coords(1,i))
       max_y = max(max_y,coords(2,i))
       min_y = min(min_y,coords(2,i))
       max_z = max(max_z,coords(3,i))
       min_z = min(min_z,coords(3,i))
    enddo

    origin(1) = min_x
    origin(2) = min_y
    origin(3) = min_z

    coords(1,:) = coords(1,:) - origin(1)
    coords(2,:) = coords(2,:) - origin(2)
    coords(3,:) = coords(3,:) - origin(3)

    coords(1,:) = mod(coords(1,:)*1000000,1000000*lattice_vectors(1,1))/1000000;
    coords(2,:) = mod(coords(2,:)*1000000,1000000*lattice_vectors(2,2))/1000000;
    coords(3,:) = mod(coords(3,:)*1000000,1000000*lattice_vectors(3,3))/1000000;

    origin = 0.0_dp

    !Get new origin.
    origin(1) = -1.0d-1 ; origin(2) = -1.0d-1; origin(3) = -1.0d-1

  end subroutine prg_translateandfoldtobox

  !> Wrap around atom i using pbc
  !! \param coords Coordinates of the system (see system_type).
  !! \param lattice_vectors System lattice vectors.
  !! \param index Index atom to wrap around
  !!
  subroutine prg_wraparound(coords,lattice_vectors,index)
    implicit none
    integer                              ::  i, nats
    integer, intent(in)                  ::  index
    real(dp), allocatable, intent(inout) ::  coords(:,:)
    real(dp), allocatable                ::  origin(:)
    real(dp), intent(in)                 ::  lattice_vectors(:,:)

    if(.not.allocated(origin)) allocate(origin(3))

    nats=size(coords,dim=2)

    origin(1) = -coords(1,index) + lattice_vectors(1,1)/2.0_dp
    origin(2) = -coords(2,index) + lattice_vectors(2,2)/2.0_dp
    origin(3) = -coords(3,index) + lattice_vectors(3,3)/2.0_dp

    coords(1,:) = coords(1,:) + origin(1)
    coords(2,:) = coords(2,:) + origin(2)
    coords(3,:) = coords(3,:) + origin(3)

    !$omp parallel do default(none) private(i) &
    !$omp shared(coords,lattice_vectors,nats)
    do i=1,nats
      if(coords(1,i) > lattice_vectors(1,1))coords(1,i)=coords(1,i)-lattice_vectors(1,1)
      if(coords(2,i) > lattice_vectors(2,2))coords(2,i)=coords(2,i)-lattice_vectors(2,2)
      if(coords(3,i) > lattice_vectors(3,3))coords(3,i)=coords(3,i)-lattice_vectors(3,3)
      if(coords(1,i) < 0.0_dp)coords(1,i)=coords(1,i)+lattice_vectors(1,1)
      if(coords(2,i) < 0.0_dp)coords(2,i)=coords(2,i)+lattice_vectors(2,2)
      if(coords(3,i) < 0.0_dp)coords(3,i)=coords(3,i)+lattice_vectors(3,3)
    enddo
    !$end omp parallel do

  end subroutine prg_wraparound


  !> Translate to geometric center.
  !! \param coords Coordinates of the system (see system_type).
  !! \param lattice_vectors System lattice vectors.
  !! \param origin (min(x),min(y),min(z)) set as the origin of the system.
  !!
  subroutine prg_translatetogeomcandfoldtobox(coords,lattice_vectors,origin)
    implicit none
    integer                              ::  i
    real(dp), allocatable, intent(inout) ::  origin(:),coords(:,:)
    real(dp), intent(in)                 ::  lattice_vectors(:,:)
    real(dp)                             ::  geomc(3)

    if(.not.allocated(origin)) allocate(origin(3))

    write(*,*)size(coords,dim=2)

    ! Getting the geometric center.
    do i=1,size(coords,dim=2)
       geomc = geomc + coords(:,i)
    enddo

    geomc = geomc/size(coords,dim=2)

    coords(1,:) = coords(1,:) - geomc(1)
    coords(2,:) = coords(2,:) - geomc(2)
    coords(3,:) = coords(3,:) - geomc(3)

    coords(1,:) = mod(coords(1,:)*1000000,1000000*lattice_vectors(1,1))/1000000;
    coords(2,:) = mod(coords(2,:)*1000000,1000000*lattice_vectors(2,2))/1000000;
    coords(3,:) = mod(coords(3,:)*1000000,1000000*lattice_vectors(3,3))/1000000;

    origin = 0.0_dp

    !Shift origin slightly
    origin(1) = -1.0d-1 ; origin(2) = -1.0d-1; origin(3) = -1.0d-1

  end subroutine prg_translatetogeomcandfoldtobox


  !> Get the volume of the cell and the reciprocal vectors:
  !! This soubroutine computes:
  !! - \f$ b_1 = \frac{1}{V_c} a_1 \times a_2 \f$
  !! - \f$ b_2 = \frac{1}{V_c} a_2 \times a_3 \f$
  !! - \f$ b_3 = \frac{1}{V_c} a_3 \times a_1 \f$
  !! - \f$ V_c = || a_1\cdot (a_2 \times a_3)||\f$
  !! - \f$ V_{BZ} = || b_1\cdot (b_2 \times b_3)||\f$
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param recip_vectors Reciprocal vectors of the system.
  !! \param volr Volume of the cell.
  !! \param volk Volume of the reciprocal cell.
  subroutine prg_get_recip_vects(lattice_vectors,recip_vectors,volr,volk)
    implicit none
    real(dp)                             ::  a1xa2(3), a2xa3(3), a3xa1(3), b2xb3(3)
    real(dp)                             ::  pi
    real(dp), allocatable, intent(inout)  ::  recip_vectors(:,:)
    real(dp), intent(in)                 ::  lattice_vectors(:,:)
    real(dp), intent(inout)              ::  volk, volr

    if(.not.allocated(recip_vectors))allocate(recip_vectors(3,3))
    volr=0.0_dp; volk=0.0_dp

    pi = 3.14159265358979323846264338327950_dp

    a1xa2(1) = lattice_vectors(1,2)*lattice_vectors(2,3) - lattice_vectors(1,3)*lattice_vectors(2,2)
    a1xa2(2) = -lattice_vectors(1,1)*lattice_vectors(2,3) + lattice_vectors(1,3)*lattice_vectors(2,1)
    a1xa2(3) =  lattice_vectors(1,1)*lattice_vectors(2,2) - lattice_vectors(1,2)*lattice_vectors(2,1)

    a2xa3(1) = lattice_vectors(2,2)*lattice_vectors(3,3) - lattice_vectors(2,3)*lattice_vectors(3,2)
    a2xa3(2) = -lattice_vectors(2,1)*lattice_vectors(3,3) + lattice_vectors(2,3)*lattice_vectors(3,1)
    a2xa3(3) =  lattice_vectors(2,1)*lattice_vectors(3,2) - lattice_vectors(2,2)*lattice_vectors(3,1)

    a3xa1(1) = lattice_vectors(3,2)*lattice_vectors(1,3) - lattice_vectors(3,3)*lattice_vectors(1,2)
    a3xa1(2) = -lattice_vectors(3,1)*lattice_vectors(1,3) + lattice_vectors(3,3)*lattice_vectors(1,1)
    a3xa1(3) =  lattice_vectors(3,1)*lattice_vectors(1,2) - lattice_vectors(3,2)*lattice_vectors(1,1)

    !Get the volume of the cell
    volr = lattice_vectors(1,1)*a2xa3(1)+ lattice_vectors(1,2)*a2xa3(2)+lattice_vectors(1,3)*a2xa3(3)

    !Get the reciprocal vectors
    recip_vectors(:,1) = 2.0_dp*pi*a2xa3(:)/volr
    recip_vectors(:,2) = 2.0_dp*pi*a3xa1(:)/volr
    recip_vectors(:,3) = 2.0_dp*pi*a1xa2(:)/volr

    b2xb3(1) = recip_vectors(2,2)*recip_vectors(3,3) - recip_vectors(2,3)*recip_vectors(3,2)
    b2xb3(2) = -recip_vectors(2,1)*recip_vectors(3,3) + recip_vectors(2,3)*recip_vectors(3,1)
    b2xb3(3) =  recip_vectors(2,1)*recip_vectors(3,2) - recip_vectors(2,2)*recip_vectors(3,1)

    volk = recip_vectors(1,1)*b2xb3(1)+ recip_vectors(1,2)*b2xb3(2)+recip_vectors(1,3)*b2xb3(3)

  end subroutine prg_get_recip_vects

  !> Get the dihedral angle given four atomic positions.
  !! \param sy System structure
  !! \param id1 Atom index 1
  !! \param id2 Atom index 1
  !! \param id3 Atom index 1
  !! \param id4 Atom index 1
  !! \param dihedral Output dihedral angle
  !!
  subroutine prg_get_dihedral(coords,id1,id2,id3,id4,dihedral)

    real(dp)                          ::  mv1, mv2, v1(3), v2(3)
    real(dp)                          ::  dotprod, cosdir, v2xv20(3), v1xv10(3)
    real(dp)                          ::  v10(3),v20(3), cprod(3), normcprod, sindir
    real(dp), intent(in)              ::  coords(:,:)
    real(dp), intent(out)             ::  dihedral
    integer                           ::  i
    integer, intent(in)               ::  id1,id2,id3,id4
    character(2)                      ::  index1, index2, index3, index4

    v1=coords(:,id4) - coords(:,id3)
    v10=coords(:,id2) - coords(:,id3)
    v2=coords(:,id1) - coords(:,id2)
    v20=coords(:,id3) - coords(:,id2)

    v1xv10(1)=v1(2)*v10(3)-v1(3)*v10(2)
    v1xv10(2)=-(v1(1)*v10(3)-v1(3)*v10(1))
    v1xv10(3)=v1(1)*v10(2)-v1(2)*v10(1)

    v2xv20(1)=v2(2)*v20(3)-v2(3)*v20(2)
    v2xv20(2)=-(v2(1)*v20(3)-v2(3)*v20(1))
    v2xv20(3)=v2(1)*v20(2)-v2(2)*v20(1)

    dotprod = v1xv10(1)*v2xv20(1) + v1xv10(2)*v2xv20(2) + v1xv10(3)*v2xv20(3)

    cprod(1)=v1xv10(2)*v2xv20(3)-v1xv10(3)*v2xv20(2)
    cprod(2)=-(v1xv10(1)*v2xv20(3)-v1xv10(3)*v2xv20(1))
    cprod(3)=v1xv10(1)*v2xv20(2)-v1xv10(2)*v2xv20(1)

    normcprod=sqrt(cprod(1)*cprod(1) + cprod(2)*cprod(2) + cprod(3)*cprod(3))

    mv1= sqrt(v1xv10(1)*v1xv10(1) + v1xv10(2)*v1xv10(2) + v1xv10(3)*v1xv10(3))
    mv2= sqrt(v2xv20(1)*v2xv20(1) + v2xv20(2)*v2xv20(2) + v2xv20(3)*v2xv20(3))

    cosdir = dotprod/(mv1*mv2)
    sindir = normcprod/(mv1*mv2)

    sindir = cprod(3)
    dihedral=sign(1.0d0,sindir)*acos(-cosdir)
    dihedral=-360*dihedral/(2.0*3.14159265359)
    if(dihedral < 0)dihedral = 360 + dihedral

  end subroutine prg_get_dihedral

  !> Get the covalency graph in bml format.
  !! \brief This is the graph composed by the covalent bonds (edges)
  !! that are determined with the VDW radius.
  !! \param sy System structure (see system_type).
  !! \param nnStructMindist Minimun distance between atoms.
  !! \param nnStruct The neigbors J to I within Rcut that are all within the box.
  !! \param nrnnstruct Number of neigbors to I within Rcut that are all within the box.
  !! \param bml_type The bml type for constructing the graph.
  !! \param gconv_bml Covanlency graph in bml format.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_get_covgraph(sy,nnStructMindist,nnStruct,nrnnstruct,bml_type,factor,gcov_bml,mdimin,verbose)
    implicit none
    character(20), intent(in)          ::  bml_type
    integer                            ::  i, j, jj, mdim
    integer, intent(in)                ::  mdimin
    integer, intent(in)                ::  nnStruct(:,:), nrnnstruct(:)
    integer, optional, intent(in)      ::  verbose
    real(dp)                           ::  d, dvdw, factor, ra(3),rb(3),rab(3)
    real(dp)                           ::  Lx, Ly, Lz
    real(dp), intent(in)               ::  nnStructMindist(:,:)
    type(bml_matrix_t), intent(inout)  ::  gcov_bml
    type(system_type), intent(in)      ::  sy
    logical(1), allocatable            ::  ispresent(:)

    if(bml_get_N(gcov_bml).gt.0) call bml_deallocate(gcov_bml)

    if(mdimin > 0)then
       mdim = mdimin
    else
       mdim = sy%nats
    endif

    call bml_zero_matrix(bml_type,bml_element_real,kind(1.0),sy%nats,mdim,gcov_bml)

    Lx = sy%lattice_vector(1,1)
    Ly = sy%lattice_vector(2,2)
    Lz = sy%lattice_vector(3,3)

    if(present(verbose).and.verbose >= 1)then
       write(*,*)" "
       write(*,*)"Building covalency graph ..."
       write(*,*)" "
    endif

    allocate(ispresent(sy%nats))

    do i=1,sy%nats
       ra = sy%coordinate(:,i)
       ispresent = .false.
       do j=1,nrnnstruct(i)

          if(nrnnstruct(i) <= 1)write(*,*)"WARNING: Atom" ,i,"is desconnected"

          jj = nnStruct(j,i)
          if(.not.ispresent(jj))then
             rb = sy%coordinate(:,jj)

             rab(1) = modulo((Rb(1) - Ra(1) + Lx/2.0_dp),Lx) - Lx/2.0_dp
             rab(2) = modulo((Rb(2) - Ra(2) + Ly/2.0_dp),Ly) - Ly/2.0_dp
             rab(3) = modulo((Rb(3) - Ra(3) + Lz/2.0_dp),Lz) - Lz/2.0_dp

             d = norm2(rab)

             if(d == 0.0_dp .and. i .ne. jj)write(*,*)"WARNING: Atom" ,i,"and atom",jj,&
                  "are on top of each other"
             dvdw = factor*(element_covr(sy%atomic_number(i)) + element_covr(sy%atomic_number(jj)))
             if(d < dvdw .and. d > 0.0_dp)then
                ispresent(jj) = .true.
                !         call bml_set_element_new(gcov_bml,i,jj,1.0_dp)
                !         call bml_set_element_new(gcov_bml,jj,i,1.0_dp)
                call bml_set_element_new(gcov_bml,i,jj,1.0)
                !         call bml_set_element_new(gcov_bml,jj,i,1.0)
             endif
          endif
       enddo

       !     call bml_set_element_new(gcov_bml,i,i,1.0_dp)
    enddo

    deallocate(ispresent)

  end subroutine prg_get_covgraph


  subroutine prg_get_covgraph_int(sy,nnStructMindist,nnStruct,nrnnstruct,bml_type,factor,gcov_bml,mdimin,verbose)
    implicit none
    character(20), intent(in)          ::  bml_type
    integer                            ::  i, j, jj, mdim
    integer, intent(in)                ::  mdimin
    integer, intent(in)                ::  nnStruct(:,:), nrnnstruct(:)
    integer, optional, intent(in)      ::  verbose
    real(dp)                           ::  d, dvdw, factor
    real(dp), intent(in)               ::  nnStructMindist(:,:)
    type(bml_matrix_t), intent(inout)  ::  gcov_bml
    type(system_type), intent(in)      ::  sy

    if(bml_get_N(gcov_bml).gt.0) call bml_deallocate(gcov_bml)

    if(mdimin > 0)then
       mdim = mdimin
    else
       mdim = sy%nats
    endif

    call bml_zero_matrix(bml_type,bml_element_real,dp,sy%nats,mdim,gcov_bml)

    if(present(verbose).and.verbose >= 1)then
       write(*,*)" "
       write(*,*)"Building covalency graph ..."
       write(*,*)" "
    endif

    do i=1,sy%nats
       do j=1,nrnnstruct(i)
          if(nrnnstruct(i) <= 1)write(*,*)"WARNING: Atom" ,i,"is desconnected"
          jj = nnStruct(j,i)
          d = nnStructMindist(j,i)
          if(d == 0.0_dp .and. i .ne. jj)write(*,*)"WARNING: Atom" ,i,"and atom",j,&
               "are on top of each other"
          dvdw = factor*(element_covr(sy%atomic_number(i)) + element_covr(sy%atomic_number(jj)))
          if(d < dvdw .and. d > 0.0_dp)then
             call bml_set_element_new(gcov_bml,i,jj,1.0_dp)
             call bml_set_element_new(gcov_bml,jj,i,1.0_dp)
          endif
       enddo
       call bml_set_element_new(gcov_bml,i,i,1.0_dp)
    enddo

  end subroutine prg_get_covgraph_int


  !> Get the covanlency graph.
  !! \brief This is the graph composed by the covalent bonds (edges)
  !! that are determined with the VDW radius.
  !! \param sy System structure (see system_type).
  !! \param nnStructMindist Minimun distance between atoms.
  !! \param nnStruct The neigbors J to I within Rcut that are all within the box.
  !! \param nrnnstruct Number of neigbors to I within Rcut that are all within the box.
  !! \param bml_type The bml type for constructing the graph.
  !! \param gconv_bml Covanlency graph in bml format.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_get_covgraph_h(sy,nnStructMindist,nnStruct,nrnnstruct,rcut,&
       graph_h,mdimin,verbose)
    implicit none
    integer                            ::  i, j, jj, ncount, mdim
    integer, intent(in)                ::  nnStruct(:,:), nrnnstruct(:), mdimin
    integer, optional, intent(in)      ::  verbose
    real(dp)                           ::  d, dvdw, ra(3),rb(3),rab(3)
    real(dp)                           ::  Lx, Ly, Lz
    real(dp), intent(in)               ::  nnStructMindist(:,:), rcut
    integer, allocatable, intent(inout)  ::  graph_h(:,:)
    type(system_type), intent(in)      ::  sy


    if(mdimin > 0)then
       mdim = mdimin
    else
       mdim = sy%nats
    endif

    if(.not.allocated(graph_h)) allocate(graph_h(mdim,sy%nats))

    if(present(verbose).and.verbose >= 1)then
       write(*,*)" "
       write(*,*)"Building H covalency graph ..."
       write(*,*)" "
    endif

    graph_h = 0

    Lx = sy%lattice_vector(1,1)
    Ly = sy%lattice_vector(2,2)
    Lz = sy%lattice_vector(3,3)

    !$omp parallel do default(none) private(i) &
    !$omp private(ncount,j,jj,ra,rb,rab,d) &
    !$omp shared(sy,nrnnstruct,nnStruct,Lx,Ly,Lz,graph_h,rcut)
    do i=1,sy%nats
       ncount = 0
       do j=1,nrnnstruct(i)
          jj = nnStruct(j,i)
          ra=sy%coordinate(:,i)
          rb=sy%coordinate(:,jj)
          rab(1) = modulo((Rb(1) - Ra(1) + Lx/2.0_dp),Lx) - Lx/2.0_dp
          rab(2) = modulo((Rb(2) - Ra(2) + Ly/2.0_dp),Ly) - Ly/2.0_dp
          rab(3) = modulo((Rb(3) - Ra(3) + Lz/2.0_dp),Lz) - Lz/2.0_dp
          d = norm2(rab)
          if(d.lt.rcut.and.d.gt.0.0_dp)then
             ncount=ncount+1
             graph_h(ncount,i) = jj
          endif
       enddo
    enddo
    !omp end parallel do

  end subroutine prg_get_covgraph_h

  !> Get a subsystem out of the total system.
  !!
  !! \brief This will get a subsystem from the total system
  !! guided by a partition.
  !!
  !! \param sy System structure (see system_type).
  !! \param lsize Core+Hallo subsystem size.
  !! \param indices Partition indices.
  !! \param sbsy Subsystem to be extracted.
  !!
  subroutine prg_get_subsystem(sy,lsize,indices,sbsy,verbose)
    implicit none
    integer                           ::  i, nsptmp, prev
    integer, intent(in)               ::  indices(:), lsize
    integer, optional, intent(in)     ::  verbose
    type(system_type), intent(in)     ::  sy
    type(system_type), intent(inout)  ::  sbsy

    if(present(verbose).and.verbose >= 1)then
       write(*,*)" "
       write(*,*)"Extracting subsystem ..."
       write(*,*)" "
    endif

    sbsy%nats = lsize

    if(allocated(sbsy%symbol))then
       deallocate(sbsy%symbol)
       deallocate(sbsy%coordinate)
       deallocate(sbsy%atomic_number)
       deallocate(sbsy%velocity)
       deallocate(sbsy%force)
       deallocate(sbsy%net_charge)
       deallocate(sbsy%mass)
       deallocate(sbsy%lattice_vector)
       deallocate(sbsy%spindex)
    endif

    allocate(sbsy%symbol(sbsy%nats))
    allocate(sbsy%coordinate(3,sbsy%nats))
    allocate(sbsy%atomic_number(sbsy%nats))
    allocate(sbsy%velocity(3,sbsy%nats))
    allocate(sbsy%force(3,sbsy%nats))
    allocate(sbsy%net_charge(sbsy%nats))
    allocate(sbsy%mass(sbsy%nats))
    allocate(sbsy%lattice_vector(3,3))
    allocate(sbsy%spindex(sbsy%nats))

    ! allocate(sbsy%recip_vector(3,3))

    sbsy%lattice_vector = sy%lattice_vector
    sbsy%volr = sy%volr
    sbsy%volk = sy%volk

    nsptmp = 0
    do i=1,sbsy%nats
       sbsy%symbol(i) = sy%symbol(indices(i)+1) !Indices from the graph partition start from 0
       sbsy%coordinate(:,i) = sy%coordinate(:,indices(i)+1)
       sbsy%atomic_number(i) = sy%atomic_number(indices(i)+1)
       sbsy%mass(i) = sy%mass(indices(i)+1)
       sbsy%spindex(i) = sy%spindex(indices(i)+1)
       if(sbsy%spindex(i).gt.nsptmp)then
          nsptmp = nsptmp + 1
       endif
    enddo

    if(nsptmp.lt.sy%nsp)then
       write(*,*)"WARNING: nsp = ",nsptmp
       write(*,*)"The subsystem contains less species that the system ..."
       write(*,*)"A generalization to parts where subsystem contains"
       write(*,*)"less species that the system needs to be added ..."
       write(*,*)"See prg_get_subsystem in prg_system_mod"
       !     stop
    endif

    sbsy%nsp = sy%nsp

    if(allocated(sbsy%spatnum))then
       deallocate(sbsy%spatnum)
       deallocate(sbsy%spmass)
    endif

    allocate(sbsy%spatnum(sbsy%nsp))
    allocate(sbsy%spmass(sbsy%nsp))

    sbsy%spatnum = sy%spatnum
    sbsy%spmass = sy%spmass

  end subroutine prg_get_subsystem


  !> Destroy allocated subsystem.
  !!
  !! \brief This routine will deallocate all the arrays of the structures.
  !!
  !! \param sy System to de deallocated (see system_type).
  !!
  subroutine prg_destroy_subsystems(sbsy,verbose)
    implicit none
    type(system_type), intent(inout)  ::  sbsy
    integer, optional, intent(in)     ::  verbose
    integer                           ::  nparts

    if(present(verbose).and.verbose >= 1)then
       write(*,*)" "
       write(*,*)"Deallocating the multiple subsystem ..."
       write(*,*)" "
    endif


    if(allocated(sbsy%symbol))then
       deallocate(sbsy%symbol)
       deallocate(sbsy%coordinate)
       deallocate(sbsy%atomic_number)
       deallocate(sbsy%velocity)
       deallocate(sbsy%force)
       deallocate(sbsy%net_charge)
       deallocate(sbsy%mass)
       deallocate(sbsy%spindex)
    endif

    if(allocated(sbsy%lattice_vector))deallocate(sbsy%lattice_vector)
    if(allocated(sbsy%resindex))deallocate(sbsy%resindex)

    if(allocated(sbsy%spatnum))then
       deallocate(sbsy%spatnum)
       deallocate(sbsy%spmass)
    endif

    if(bml_get_N(sbsy%estr%ham) > 0)   call bml_deallocate(sbsy%estr%ham)
    if(bml_get_N(sbsy%estr%ham0) > 0)  call bml_deallocate(sbsy%estr%ham0)
    if(bml_get_N(sbsy%estr%oham) > 0)  call bml_deallocate(sbsy%estr%oham)
    if(bml_get_N(sbsy%estr%over) > 0)  call bml_deallocate(sbsy%estr%over)
    if(bml_get_N(sbsy%estr%rho) > 0)   call bml_deallocate(sbsy%estr%rho)
    if(bml_get_N(sbsy%estr%orho) > 0)  call bml_deallocate(sbsy%estr%orho)
    if(bml_get_N(sbsy%estr%zmat) > 0)  call bml_deallocate(sbsy%estr%zmat)


    if(allocated(sbsy%estr%coul_pot_r)) deallocate(sbsy%estr%coul_pot_r)
    if(allocated(sbsy%estr%coul_pot_k)) deallocate(sbsy%estr%coul_pot_k)
    if(allocated(sbsy%estr%skforce))    deallocate(sbsy%estr%skforce)
    if(allocated(sbsy%estr%fpul))       deallocate(sbsy%estr%fpul)
    if(allocated(sbsy%estr%fscoul))     deallocate(sbsy%estr%fscoul)
    if(allocated(sbsy%estr%hindex))     deallocate(sbsy%estr%hindex)


  end subroutine prg_destroy_subsystems

  !> Partition by molecule.
  !!
  !! \param sy System structure.
  !! \param npart Number of parts.
  !! \param nnStructMindist Minimum distance between neighbors.
  !! \param nnStruct The neighbors J to I within Rcut that are all within the box.
  !! \param nrnnstruct Number of neighbors to I within Rcut that are all within the box.
  !! \param hetatm Atom to be taken as the "center" of the by molecule partition.
  !! \param gp Graph partition structure.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_molpartition(sy,npart,nnStructMindist,nnStruct,nrnnstruct,hetatm,gp,verbose)
    implicit none
    character(2), intent(in)                   ::  hetatm
    integer                                    ::  cnt, i, j, jj
    integer                                    ::  maxparts, nmax
    integer, allocatable                       ::  core_count(:), cores(:,:), partof(:)
    integer, intent(in)                        ::  nnStruct(:,:), nrnnstruct(:)
    integer, intent(inout)                     ::  npart
    integer, optional, intent(inout)           ::  verbose
    logical, allocatable                       ::  inpart(:)
    real(dp)                                   ::  d, dvdw
    real(dp), intent(in)                       ::  nnStructMindist(:,:)
    type(graph_partitioning_t), intent(inout)  ::  gp
    type(system_type), intent(in)              ::  sy

    if(present(verbose).and.verbose >= 1)then
       write(*,*)" "; write(*,*)"Partitioning by molecule ..."; write(*,*)" "
    endif

    if(allocated(gp%nnodesInPart))then
       call prg_destroyGraphPartitioning(gp)
    endif

    nmax = sy%nats
    maxparts = sy%nats

    allocate(core_count(maxparts))

    allocate(cores(nmax,maxparts))
    allocate(inpart(sy%nats))
    allocate(partof(sy%nats))

    cores = -1
    inpart = .false.
    partof=-1
    core_count = 0
    npart = 0
    do i=1,sy%nats
       if(.not.inpart(i).and.trim(sy%symbol(i)) == trim(hetatm))then
          npart = npart + 1
          core_count(npart) = core_count(npart) + 1
          cores(core_count(npart),npart) = i
          partof(i) = npart
          do j=1,nrnnstruct(i)
             jj = nnStruct(j,i)
             d = nnStructMindist(j,i)
             dvdw = 1.3_dp
             if(d.lt.dvdw.and.d.gt.0.001_dp)then
                if(.not.inpart(jj))then
                   core_count(npart) = core_count(npart) + 1
                   cores(core_count(npart),npart) = jj
                   inpart(jj) = .true.
                   partof(jj) = npart
                endif
             endif
          enddo
       endif
    enddo

    call prg_initGraphPartitioning(gp, "molecules", npart, sy%nats, sy% nats)

    ! Initialize and fill up subgraph structure
    ! Assign node ids (mapped to orbitals as rows) to each node in each
    do i = 1, npart
       gp%nnodesInPartAll(i) = core_count(i)
       call prg_initSubgraph(gp%sgraph(i), i, gp%totalNodes2)
       allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
       gp%nnodesInPart(i) = core_count(i)
    enddo

    !Assign node ids to sgraph
    do i=1, npart
       do j=1,core_count(i)
          gp%sgraph(i)%nodeInPart(j) = cores(j,i) - 1
       end do
    enddo

  end subroutine prg_molpartition

  !> Get partial subgraph based on the Density matrix.
  !!
  !! \param rho_bml Density matix in bml format.
  !! \param hindex Start and end index for every atom in the system.
  !! \param gch_bml Atom based graph in bml format.
  !! \param threshold Threshold value for constructing the graph.
  !! \param verbose Verbosity levels.
  !!
  subroutine prg_get_partial_atomgraph(rho_bml,hindex,gch_bml,threshold,verbose)
    implicit none
    character(20)                      ::  bml_type
    integer                            ::  i, ii, j, jj
    integer                            ::  nch, norbs
    integer, intent(in)                ::  hindex(:,:)
    integer, optional, intent(in)      ::  verbose
    logical                            ::  connection
    logical, allocatable               ::  iconnectedtoj(:)
    real(dp), allocatable              ::  row(:), rowat(:)
    real(dp), intent(in)               ::  threshold
    type(bml_matrix_t), intent(in)     ::  rho_bml
    type(bml_matrix_t), intent(inout)  ::  gch_bml

    norbs = bml_get_N(rho_bml)
    bml_type = bml_get_type(rho_bml)
    nch = size(hindex,dim=2)
    allocate(rowat(nch))
    allocate(row(norbs))
    allocate(iconnectedtoj(nch))

    if(bml_get_N(gch_bml) > 0) call bml_deallocate(gch_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,nch,nch,gch_bml)

    do i = 1, nch
       rowat = 0.0_dp
       do ii = hindex(1,i),hindex(2,i)  ! i,j block
          call bml_get_row(rho_bml,ii,row)
          iconnectedtoj = .false.
          do j = 1, nch
             connection = .false.

             if(.not.iconnectedtoj(j).and.rowat(j) == 0)then
                do jj = hindex(1,j),hindex(2,j)
                   if(abs(row(jj)) > threshold)then
                      connection = .true.  !We exit if there is a connection.
                      exit
                   endif
                enddo
                if(connection) iconnectedtoj(j) = .true.
             endif
             if(connection) rowat(j) = 1.0_dp
          enddo

       enddo
       call bml_set_row(gch_bml,i,rowat)
    enddo

    deallocate(rowat)
    deallocate(row)
    deallocate(iconnectedtoj)

  end subroutine prg_get_partial_atomgraph


  !> Collect the small graph to build the full graph.
  !!
  !! \param rho_bml Density matix in bml format.
  !! \param nc Number of core atoms.
  !! \param nats Number of atoms.
  !! \param hindex Hindex for the small part (see haindex)
  !! \param chindex Core-hallo index for the small part.
  !! \param graph_p Graph in an "ellpack" format.
  !! \param threshold Threshold to buil the density based atom projected graph.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_collect_graph_p(rho_bml,nc,nats,hindex,chindex,graph_p,threshold,mdimin,verbose)
    implicit none
    character(20)                       ::  bml_type
    integer                             ::  i, ifull, ii, j
    integer                             ::  jfull, jj, nch, ncounti
    integer                             ::  norbs, mdim
    logical(1), allocatable             ::  rowatfull(:)
    integer, allocatable, intent(inout) ::  graph_p(:,:)
    integer, intent(in)                 ::  chindex(:), hindex(:,:), nats, nc
    integer, intent(in)                 ::  mdimin
    integer, optional, intent(in)       ::  verbose
    logical                             ::  connection
    logical, allocatable                ::  iconnectedtoj(:)
    real(dp), allocatable               ::  row(:)
    real(dp), intent(in)                ::  threshold
    type(bml_matrix_t), intent(in)      ::  rho_bml

    norbs = bml_get_N(rho_bml)
    bml_type = bml_get_type(rho_bml)
    nch = size(hindex,dim=2)
    allocate(rowatfull(nats))
    allocate(row(norbs))
    allocate(iconnectedtoj(nch))

    if(mdimin > 0)then
       mdim = mdimin
    else
       mdim = nats
    endif

    if(.not.allocated(graph_p))then
       allocate(graph_p(mdim,nats))
       graph_p = 0
    endif

    !$omp parallel do default(none) private(i) &
    !$omp private(ncounti,rowatfull,ii,j,jfull,ifull,iconnectedtoj) &
    !$omp private(row,connection) &
    !$omp shared(graph_p,nc,chindex,hindex,rho_bml,nch,threshold)
    do i = 1, nc

       ifull = chindex(i) + 1 !Map it to the full system
       rowatfull = .false.
       ii=1
       ncounti = 0
       do while (graph_p(ii,ifull).ne.0)  !Unfolding the connections of ifull
          rowatfull(graph_p(ii,ifull)) = .true.
          ncounti = ncounti + 1
          ii = ii+1
       enddo

       do ii = hindex(1,i),hindex(2,i)  ! i,j block
          call bml_get_row(rho_bml,ii,row)
          iconnectedtoj = .false.
          !
          do j = 1, nch
             jfull = chindex(j) + 1  !Cause the count starts form 0
             connection = .false.
             if(.not.iconnectedtoj(j).and..not.rowatfull(jfull))then
                do jj = hindex(1,j),hindex(2,j)
                   if(abs(row(jj)) > threshold)then
                      connection = .true.  !We exit if there is a connection.
                      exit
                   endif
                enddo
                if(connection) iconnectedtoj(j) = .true.
             endif
             if(connection)then  !I conected to j >>> ifull connected to jfull
                !Map back to the full system graph.
                rowatfull(jfull) = .true.
                ncounti = ncounti + 1
                graph_p(ncounti,ifull) = jfull
             endif
          enddo
          !
       enddo

    enddo
    !$omp end parallel do

    deallocate(rowatfull)
    deallocate(iconnectedtoj)
    deallocate(row)

  end subroutine prg_collect_graph_p


  !> Get partial subgraph based on the Density matrix.
  !!
  !! \param graph_p Density matix based graph in bml format.
  !! \param graph_h Hamiltonian matix based graph in bml format.
  !!
  subroutine prg_merge_graph(graph_p,graph_h)
    implicit none
    integer                               ::  i, ii, j, nats
    integer                               ::  ncounti, maxnz
    integer,              intent(inout)   ::  graph_h(:,:)
    integer, intent(inout)                ::  graph_p(:,:)
    logical(1), allocatable               ::  rowpatfull(:)

    nats = size(graph_p,dim=2)
    maxnz = size(graph_p,dim=1)

    allocate(rowpatfull(nats))

    !$omp parallel do default(none) private(i) &
    !$omp private(ncounti,rowpatfull,ii,j) &
    !$omp shared(graph_p,graph_h,nats,maxnz)
    do i = 1, nats
       ncounti = 0
       rowpatfull = .false.

       do ii = 1,maxnz
          if(graph_p(ii,i) == 0) exit
          rowpatfull(graph_p(ii,i)) = .true.
          ncounti = ncounti + 1
       enddo

       do ii = 1,maxnz
          j = graph_h(ii,i)
          if(j==0)exit
          if(.not.rowpatfull(j))then
             ncounti = ncounti + 1
             graph_p(ncounti,i) = j
          endif
       enddo

    enddo
    !$omp end parallel do

    deallocate(rowpatfull)


  end subroutine prg_merge_graph


  !> Get partial subgraph based on the Density matrix.
  !!
  !! \param graph_p Density matix based graph in "ellpack type format".
  !! \param graph_h Hamiltonian matix based graph in "ellpack type format".
  !! \param xadj CSR start values for the adjacency matrix.
  !! \param adjncy CSR positions of adjacency matrix.
  !!
  subroutine prg_merge_graph_adj(graph_p,graph_h,xadj,adjncy)
    implicit none
    integer                             ::  i, ii, j, nats
    integer                             ::  ncounti, ncountot
    integer, allocatable, intent(inout)  ::  adjncy(:), graph_h(:,:), graph_p(:,:), xadj(:)
    logical(1), allocatable             ::  rowpatfull(:)

    nats = size(graph_p,dim=2)
    allocate(rowpatfull(nats))
    allocate(xadj(nats+1))
    allocate(adjncy(nats*nats))

    xadj = 0.0
    adjncy = 0.0

    !$omp parallel do default(none) private(i) &
    !$omp private(ncounti,rowpatfull,ii,j,ncountot) &
    !$omp shared(graph_p,graph_h,nats)
    do i = 1, nats
       ncounti = 0
       rowpatfull = .false.
       ii = 1
       do while (graph_p(ii,i).ne.0)  !Unfolding the connections of i
          rowpatfull(graph_p(ii,i)) = .true.
          ncounti = ncounti + 1
          ii = ii+1
       enddo

       ii = 1
       do while (graph_h(ii,i).ne.0)  !Unfolding the connections of i
          j = graph_h(ii,i)
          if(.not.rowpatfull(j))then
             ncounti = ncounti + 1
             ncountot = ncountot + 1
             graph_p(ncounti,i) = j
          endif
          ii = ii+1
       enddo
    enddo
    !$omp end parallel do


    ncountot = 0

    !$omp parallel do default(none) private(i) &
    !$omp private(ncounti,ii,ncountot) &
    !$omp shared(graph_p,nats,xadj,adjncy)
    do i = 1, nats
       ncounti = 0
       ii = 1
       do while (graph_p(ii,i).gt.0)
          ncounti = ncounti + 1
          ncountot = ncountot + 1
          adjncy(ncountot) = graph_p(ii,i)
          ii = ii+1
       enddo
       xadj(i) = ncountot - ncounti + 1
    enddo
    !$omp end parallel do

    xadj(nats+1) = ncountot + 1

    deallocate(rowpatfull)
    deallocate(graph_p)

  end subroutine prg_merge_graph_adj

  !> prg_adj2bml
  !!
  !! \param xadj CSR start values for the adjacency matrix.
  !! \param adjncy CSR positions of adjacency matrix.
  !! \param bml_type bml format.
  !! \param g_bml graph in bml format.
  !!
  subroutine prg_adj2bml(xadj,adjncy,bml_type,g_bml)
    implicit none
    character(20), intent(in)          ::  bml_type
    integer                            ::  i, ii, j, nats
    integer, intent(in)                ::  adjncy(:), xadj(:)
    real(dp), allocatable              ::  row(:)
    type(bml_matrix_t), intent(inout)  ::  g_bml

    nats = size(xadj,dim=1) - 1
    if(bml_get_N(g_bml).gt.0) call bml_deallocate(g_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,nats,nats,g_bml)
    allocate(row(nats))

    !$omp parallel do default(none) private(i) &
    !$omp private(j,row) &
    !$omp shared(nats,g_bml,xadj,adjncy)
    do i = 1, nats
       row = 0.0_dp
       do j = xadj(i), xadj(i+1) - 1
          row(adjncy(j)) = 1.0_dp
       enddo
       call bml_set_row(g_bml,i,row)
    enddo

    deallocate(row)

  end subroutine prg_adj2bml

  !> Graph2bml
  !!
  !! \param graph Atom based graph in "ellpack" like format.
  !! \param bml_type Bml type (usually ellpack for graph starage)
  !! \param g_bml Graph in bml format.
  !!
  subroutine prg_graph2bml(graph,bml_type,g_bml)
    implicit none
    character(20), intent(in)             ::  bml_type
    integer                               ::  i, ii, nats, mdim
    integer, allocatable, intent(inout)   ::  graph(:,:)
    !   real(dp), allocatable                 ::  row(:)
    real(4), allocatable                 ::  row(:)
    type(bml_matrix_t), intent(inout)     ::  g_bml

    nats = size(graph,dim=2)
    !   if(bml_get_N(g_bml).GT.0) call bml_deallocate(g_bml)
    !   call bml_zero_matrix(bml_type,bml_element_real,kind(1.0),nats,nats,g_bml)
    allocate(row(nats))

    mdim = bml_get_m(g_bml)

    !$omp parallel do default(none) private(i) &
    !$omp private(ii,row) &
    !$omp shared(nats,g_bml,graph,mdim)
    do i = 1, nats
       ii = 1
       row = 0.0
       do while (graph(ii,i).gt.0)
          row(graph(ii,i)) = 1.0
          !       call bml_set_element_new(g_bml,i,graph(ii,i),1.0)
          !       call bml_set_element_new(g_bml,graph(ii,i),i,1.0)
          ii = ii+1
       enddo
       call bml_set_row(g_bml,i,row,0.5_dp)
    enddo
    !$omp end parallel do

    deallocate(graph)
    deallocate(row)

  end subroutine prg_graph2bml


  !> Vectorize graph.
  !!
  !! \param graph Ellpack graph.
  !! \param vector Vector to store the graph.
  !!
  subroutine prg_graph2vector(graph,vector,maxnz)
    implicit none
    integer, intent(inout)                ::  graph(:,:)
    integer, allocatable                 ::  vector(:)
    integer :: i, j, maxnz, nats, ncount

    nats = size(graph,dim=2)
    allocate(vector(nats*maxnz))

    ncount = 0

    !$omp parallel do default(none) private(i) &
    !$omp private(ncount,j) &
    !$omp shared(nats,maxnz,vector,graph)
    do i = 1, nats
       do j = 1, maxnz
          ncount = (i-1)*maxnz + j !ncount + 1
          vector(ncount) = graph(j,i)
       enddo
    enddo
    !$omp end parallel do

  end subroutine prg_graph2vector

  !> Back to graph.
  !!
  !! \param vector Vector to store the graph.
  !! \param graph Ellpack graph.
  !!
  subroutine prg_vector2graph(vector,graph,maxnz)
    implicit none
    integer, intent(inout)                ::  graph(:,:)
    integer, allocatable, intent(inout)   ::  vector(:)
    integer :: i, j, maxnz, nats, ncount

    nats = size(graph,dim=2)
    ncount = 0

    !$omp parallel do default(none) private(i) &
    !$omp private(ncount,j) &
    !$omp shared(nats,maxnz,vector,graph)
    do i = 1, nats
       do j = 1, maxnz
          ncount = (i-1)*maxnz + j !ncount + 1
          graph(j,i) = vector(ncount)
       enddo
    enddo
    !$omp end parallel do

    deallocate(vector)

  end subroutine prg_vector2graph

  !> Sort adj
  !! NOTE: this might not be needed anymre since the bml_get_adj routine is sorting
  !! the values.
  !!
  subroutine prg_sortadj(xadj, adjncy)
    implicit none
    integer, intent(inout)   ::  xadj(:)
    integer, allocatable, intent(inout)   ::  adjncy(:)
    integer :: i, j, N, l, first, last, nx, irank, mintmp
    integer, allocatable :: tmpvect(:), newadjncy(:)

    nx = size(xadj,dim=1)
    allocate(tmpvect(nx))
    allocate(newadjncy(xadj(nx) - 1))

    !$omp parallel do default(none) private(i) &
    !$omp private(first,last,j,l,tmpvect,N,mintmp) &
    !$omp shared(nx,adjncy,xadj,newadjncy)
    do i = 1, nx-1
       first = xadj(i)
       last = xadj(i+1) - 1
       write(*,*)first,last
       N = last - first + 1
       tmpvect(1:N) = 0
       tmpvect(1:N) = adjncy(first:last)

       do j = 1,N
          do l = j+1,N
             if(tmpvect(j) >= tmpvect(l).and.tmpvect(l).ne.0)then
                mintmp = tmpvect(l)
                tmpvect(l) = tmpvect(j)
                tmpvect(j) = mintmp
             endif
          enddo
       enddo

       adjncy(first:last) = tmpvect(1:N)
       newadjncy(first:last) = tmpvect(1:N)

    enddo
    !$omp end parallel do
    deallocate(tmpvect)
    deallocate(adjncy)

    allocate(adjncy(xadj(nx) - 1))  !Reallocation for adjusting size

    adjncy = newadjncy

    deallocate(newadjncy)

  end subroutine prg_sortadj

end module prg_system_mod
