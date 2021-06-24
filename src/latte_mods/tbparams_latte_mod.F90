!> A module to read and handle the latte TB params and electronic structure data.
!! @ingroup LATTE
!! \brief This module will be used load the parameter containded in:
!! \verbatim
!! /latteTBparams/electrons.dat
!! /latteTBparams/ppots.ortho
!! /latteTBparams/bondint.ortho
!! /latteTBparams/ppots.nonortho
!! /latteTBparams/bondint.nonortho
!! \endverbatim
!!
module tbparams_latte_mod

  use prg_openfiles_mod
  use bml

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> General TB params type
  !!
  type, public :: tbparams_type

    !> To have different verbosity levels.
    integer :: verbose

    !> Number of species or number of differet atom types (symbols) in the system.
    !! This is less or equal than the total number of atoms.
    integer :: nsp

    !> Element symbol for every species.
    !! A list with the different species e.g. H, C, N, etc with the order corresponding to the appearence in
    !! system@symbol.
    !! Allocation:
    !! \verbatim splist(nsp) \endverbatim
    !! This is information is also containded in system type:
    !! \verbatim tbparams%splist = system%splist \endverbatim
    !! \note In some cases the species list of system an tbparams could differ if we allow
    !! for a particular parameterization.
    character(2), allocatable :: splist(:)

    !> Basis kind for every species (s, sp, spd, spdf ...)
    !! This vector contains the type of basis set for every species.
    !! Allocation:
    !! \verbatim basis(nsp) \endverbatim
    !! Example: O atom will have "s" and "p" orbitals, hence, if species number 3 is an
    !! O atom:
    !! \verbatim tbparams%basis(3) = "sp" \endverbatim
    character(4), allocatable :: basis(:)

    !> Number of electrons for every species.
    !!To get the total number of electrons of the system we can proceed as follows:
    !! \verbatim sum(numel(system%spindex(i)),i=1,system%nats) \endverbatim
    real(dp), allocatable :: numel(:)

    !> Ionization energies (Es, Ep, Ed, Ef).
    !! This is read from the tbparams table.
    !! Example to get Ep of atom 30:
    !! \verbatim Ep = tbparams%onsite_energ(system%spindex(30)) \endverbatim
    real(dp), allocatable :: onsite_energ(:,:)

    !> Number of orbitals of species i.
    !! In a tipical TB model O species have four orbitas, one "s" and 3 "p" (px,py,pz).
    !! Example to get the number of orbitals of atom 30:
    !! \verbatim norbi(system%spindex(30)) \endverbatim
    integer, allocatable :: norbi(:)

    !> Atomic Mass.
    !! This data can also be gather from ptable module.
    !! Example:
    !! \verbatim tbparams%mass(i) = ptable%mass(tbparams%atomic_number(i))
    real(dp), allocatable :: mass(:)

    !> Hubbard parameter U.
    !! This is the onsite e-e potential repulsion.
    !! it runs over the species list.
    real(dp), allocatable :: HubbardU(:)

    !> Spin related parameters.
    real(dp), allocatable :: w(:,:)

  end type tbparams_type

  !> Type to store the Integral paramenters for every pair.
  !! This type will store all the possible pair bond integrals
  !! Allocation: If intPairs is of type intPairs_type, then
  !! \verbatim intPairs(nsp,nsp) \endverbatim
  !! and for any particular pair of species i,j then:
  !! \verbatim intPairs(i,j)%intParams(maxparams,maxint) \endverbatim
  !! Where maxint is the maximum number of diferent integral's kind.
  !! - Up to now there are only a total of four integral kinds: sss,sps,pps and ppp.
  !! - Example: If we want to know the parameter number 10 for the integral pps of the pair i,j we can do:
  !! \verbatim  intPairs(i,j)%intParams(10,3)
  type, public :: intPairs_type
    real(dp),allocatable :: intParams(:,:)
  end type intPairs_type

  !> Type to store the Pair potential paramenters for every pair in the system.
  !! This type will store all the possible pair potential.
  !! Allocation: If ppot is of type ppot_type, then
  !! \verbatim ppot(nsp,nsp) \endverbatim
  !! and for any particular pair of species i,j then:
  !! \verbatim ppot(i,j)%potparams(maxparams) \endverbatim
  !! Where maxparams is the maximum number of diferent parameters.
  !! - Example: If we want to know the parameter number 10 for the pair i,j we can do:
  !! \verbatim  ppot(i,j)%potparams(10)
  type, public :: ppot_type
    real(dp),allocatable :: potparams(:)
  end type ppot_type

  public :: load_latteTBparams, load_bintTBparamsH, write_latteTBparams, write_bintTBparamsH
  public :: load_PairPotTBparams

contains

  !> This routine loads element properties TB parameters.
  !! \param latteTBparams Latte tight-bingind parameters type.
  !! \param splist Element symbol for all the species in the system.
  !! \param parampath Path to the LATTE parameters.
  !! \note The number of species is the total number of different atoms symbols
  !! in the system.
  !! splist can be gather from the system type:
  !! \verbatim spList(:) = system%splist(:) \endverbatim
  subroutine load_latteTBparams(latteTBparams,splist,paramPath)
    implicit none
    character(100)                      ::  io_name, tableformat
    character(2)                        ::  dummy, tablesymbol
    character(2), intent(in)            ::  splist(:)
    character(4)                        ::  tablesbasis
    character(len=*)                    ::  paramPath
    integer                             ::  i, io_unit, j, maxre
    integer                             ::  noelem, nsp
    real(dp), allocatable               ::  tablevector(:)
    type(tbparams_type), intent(inout)  ::  latteTBparams

    io_name=trim(parampath)//"/electrons.dat"
    write(*,*)""; write(*,*)"Reading parameters from: ",io_name
    call prg_open_file_to_read(io_unit,io_name)
    read(io_unit,*)dummy,noelem
    read(io_unit,*)dummy
    maxre = 11
    allocate(tablevector(maxre))

    nsp = size(splist,dim=1)
    latteTBparams%nsp = nsp

    allocate(latteTBparams%splist(nsp))
    allocate(latteTBparams%numel(nsp))
    allocate(latteTBparams%basis(nsp))
    allocate(latteTBparams%onsite_energ(4,nsp))
    allocate(latteTBparams%mass(nsp))
    allocate(latteTBparams%HubbardU(nsp))
    allocate(latteTBparams%w(4,nsp))
    allocate(latteTBparams%norbi(nsp))

    latteTBparams%splist = ""

    do i=1,noelem
      tableformat='(A2,A4,11F10.5)'
      read(io_unit,*)tablesymbol,tablesbasis,(tablevector(j),j=1,maxre)
      write(*,tableformat)tablesymbol,tablesbasis,(tablevector(j),j=1,maxre)
      !Comparing the sybols in splist with the symbols from the table
      do j=1,nsp
        if(adjustl(trim(splist(j))).eq.adjustl(trim(tablesymbol)))then
          latteTBparams%splist(j)=trim(tablesymbol)
          latteTBparams%basis(j)=tablesbasis
          latteTBparams%numel(j)=tablevector(1)
          latteTBparams%onsite_energ(1,j) = tablevector(2) !Es
          latteTBparams%onsite_energ(2,j) = tablevector(3) !Ep
          latteTBparams%onsite_energ(3,j) = tablevector(4) !Ed
          latteTBparams%onsite_energ(4,j) = tablevector(5) !Ef
          latteTBparams%mass(j) = tablevector(6)
          latteTBparams%HubbardU(j) = tablevector(7)
          latteTBparams%w(1,j) = tablevector(8)  !Wss
          latteTBparams%w(2,j) = tablevector(9)   !Wpp
          latteTBparams%w(3,j) = tablevector(10)  !Wdd
          latteTBparams%w(4,j) = tablevector(11)  !Wff
        endif
      enddo
    enddo

    !Check missing species
    do i=1,nsp
      if(latteTBparams%splist(i).EQ."")then
        write(*,*) "WARNING: Species",splist(i)," has no matching parameters."
        write(*,*) "Please provide paramenter for",splist(i)," in",io_name
        write(*,*)i,latteTBparams%splist(i),splist(i)
        stop ""
      endif
    enddo

    !Constructing the number of orbitals per atom information.
    do i=1,nsp
      select case(trim(latteTBparams%basis(i)))
      case("s")
        latteTBparams%norbi(i)=1
      case("sp")
        latteTBparams%norbi(i)=4
      case("spd")
        latteTBparams%norbi(i)=9
      case("spdf")
        latteTBparams%norbi(i)=16
      case default
        write(*,*)"No basis set assigned to", latteTBparams%norbi(i)
        stop ""
      end select
    enddo

  end subroutine load_latteTBparams

  !> This routine loads the bond integral TB parameters for the hamiltonian.
  !! \param splist List of different species.
  !! \param onsite_energ Onsites energies for all the species (tbparams%onsite_energ)
  !! \param typeA First symbol for every pair combination.
  !! \param typeB First symbol for every pair combination.
  !! In particular:
  !! - \verbatim typeA(i,j) = splist(i) \endverbatim
  !! - \verbatim typeB(i,j) = splist(j) \endverbatim
  !! \param intKind Integral kinds (intKind(1)="sss"(\f$ \sigma_{s-s}\f$))
  !! \param onsitesH onsitesS Onsite energies for every orbital of a particular species.
  !! Example2: If pair system%sybol(i) = "O" then:
  !! \verbatim onsites(1,stem%spindex(i)) = Es \endverbatim
  !! \verbatim onsites(2,stem%spindex(i)) = Ep \endverbatim
  !! \verbatim onsites(3,stem%spindex(i)) = Ep \endverbatim
  !! \verbatim onsites(4,stem%spindex(i)) = Ep \endverbatim
  !! \param intPairs Integral parameters for every pair (see intPairs_type)
  !! \param parampath Path to the bond integral hamiltonian/overlap LATTE parameters.
  subroutine load_bintTBparamsH(spList,onsite_energ,typeA,typeB,&
       intKind,onsitesH,onsitesS,intPairsH,intPairsS,paramPath)
    implicit none
    character(100)                                  ::  io_name, tableformat
    character(2)                                    ::  dummy
    character(2), allocatable                       ::  pairTableType(:,:)
    character(2), allocatable, intent(inout)        ::  TypeA(:,:), TypeB(:,:)
    character(2), intent(in)                        ::  spList(:)
    character(3), allocatable                       ::  intKind(:), intTableKind(:)
    character(len=*), intent(in)                    ::  paramPath
    integer                                         ::  i, io_unit, j, k
    integer                                         ::  l, m, maxints, maxre
    integer                                         ::  npoints, nsp
    real(dp), allocatable                           ::  onsitesH(:,:), onsitesS(:,:), tablevector(:,:)
    real(dp), intent(in)                            ::  onsite_energ(:,:)
    type(intpairs_type), allocatable, intent(inout)  ::  intPairsH(:,:), intPairsS(:,:)


    io_name=trim(parampath)//"/bondints.nonortho"
    write(*,*)""; write(*,*)"Reading parameters from ",io_name
    call prg_open_file_to_read(io_unit,io_name)

    read(io_unit,*)dummy,npoints
    read(io_unit,*)dummy
    write(*,*)dummy,npoints
    maxre = 16
    allocate(tablevector(maxre,npoints))
    allocate(intTableKind(npoints))
    allocate(pairTableType(2,npoints))
    !
    tableformat='(A2,A4,A3,16F10.5)'
    do i=1,npoints
      read(io_unit,*)pairTableType(1,i),pairTableType(2,i),intTableKind(i),(tablevector(j,i),j=1,maxre)
      write(*,tableformat)pairTableType(1,i),pairTableType(2,i),intTableKind(i),(tablevector(j,i),j=1,maxre)
    enddo

    maxints = 4
    nsp = size(splist,dim=1)
    allocate(TypeA(nsp,nsp))
    allocate(TypeB(nsp,nsp))
    allocate(intKind(maxints))

    allocate(intPairsH(nsp,nsp))
    allocate(intPairsS(nsp,nsp))
    allocate(onsitesH(maxints,nsp))
    allocate(onsitesS(maxints,nsp))

    !>Integral types
    !! \todo add all the possiible integrals including d and f orbitals.
    intKind(1) = "sss"
    intKind(2) = "sps"
    intKind(3) = "pps"
    intKind(4) = "ppp"

    onsitesH = 0.0_dp
    onsitesS = 0.0_dp

    !     write(*,*)"Get All pairs"
    !Get All pairs
    do i=1,nsp
      !       write(*,*)i
      do j=1,nsp
        TypeA(i,j) = splist(i)
        TypeB(i,j) = splist(j)
        allocate(intPairsH(i,j)%intParams(maxre,maxints))
        allocate(intPairsS(i,j)%intParams(maxre,maxints))
        intPairsH(i,j)%intParams=0.0_dp
        intPairsS(i,j)%intParams=0.0_dp
        do k=1,maxints
          do l=1,npoints
            if(adjustl(trim(TypeA(i,j))).EQ.adjustl(trim(pairTableType(1,l)))&
                 .and.adjustl(trim(TypeB(i,j))).EQ.adjustl(trim(pairTableType(2,l))))then
              if(adjustl(trim(intKind(k))).EQ.adjustl(trim(intTableKind(l))))then
                do m=1,8 !Loop over all the param for the pair type and the kind
                  intPairsH(i,j)%intParams(m,k) = tablevector(m,l)    !H paramenters
                  intPairsS(i,j)%intParams(m,k) = tablevector(m+8,l)  !S paramenters
                enddo
              endif
            endif
            if(adjustl(trim(TypeA(i,j))).eq.adjustl(trim(TypeB(i,j))))then
              select case(adjustl(trim(intKind(k))))
              case("sss")
                onsitesH(k,i) = onsite_energ(1,i)
                onsitesS(k,i) = 1.0_dp
              case("sps")
                onsitesH(k,i) = onsite_energ(2,i)
                onsitesS(k,i) = 1.0_dp
              case("pps")
                onsitesH(k,i) = onsite_energ(2,i)
                onsitesS(k,i) = 1.0_dp
              case("ppp")
                onsitesH(k,i) = onsite_energ(2,i)
                onsitesS(k,i) = 1.0_dp
              case default
                stop "No valid onsite"
              end select
            endif
          enddo
          call scale_tail(intPairsH(i,j)%intParams(:,k))
          call scale_tail(intPairsS(i,j)%intParams(:,k))
        enddo
      enddo
    enddo

    !symmetrization
    do i=1,nsp
      do j=1,nsp
        do k=1,maxints
          if(sum(abs(intPairsH(i,j)%intParams(:,k))) .lt. 0.0001_dp)then
            intPairsH(i,j)%intParams(:,k) = intPairsH(j,i)%intParams(:,k)
          endif
          write(*,'(2I2,A3,100F10.5)')i,j,intKind(k),(intPairsH(i,j)%intParams(l,k),l=1,maxre)
          if(sum(abs(intPairsS(i,j)%intParams(:,k))) .lt. 0.0001_dp)then
            intPairsS(i,j)%intParams(:,k) = intPairsS(j,i)%intParams(:,k)
          endif
          write(*,'(2I2,A3,100F10.5)')i,j,intKind(k),(intPairsS(i,j)%intParams(l,k),l=1,maxre)
        enddo
      enddo
    enddo

  end subroutine load_bintTBparamsH

  !> This routine loads the pair potential latte TB parameters.
  !! \param parampath path to parameters file.
  !! \param splist Species list. Runs through all the species index.
  !! \param ppot Structuee with the pair potential.
  !! Allocation: If ppot is of type ppot_type, then
  !! \verbatim ppot(nsp,nsp) \endverbatim
  !! and for any particular pair of species i,j then:
  !! \verbatim ppot(i,j)%params(maxparams) \endverbatim
  !! Where in this case maxparams = 16.
  !!
  subroutine load_PairPotTBparams(parampath,splist,ppot)
    implicit none
    character(1)                                ::  dummy
    character(100)                              ::  io_name, tableformat
    character(2)                                ::  PairPotType1, PairPotType2
    character(2), intent(in)                    ::  splist(:)
    character(len=*), intent(in)                ::  paramPath
    integer                                     ::  i, io_unit, j, k
    integer                                     ::  npots, nsp
    real(dp)                                    ::  DDPOLY, DELTA, DELTA2, DELTA3
    real(dp)                                    ::  DELTA4, DPOLY, EXPTMP, POLY
    real(dp)                                    ::  PotCoef(16), R1, R1SQ, RCUT
    real(dp)                                    ::  SCL_R1
    type(ppot_type), allocatable, intent(inout)  ::  ppot(:,:)

    io_name=trim(parampath)//"/ppots.nonortho"
    write(*,*)""; write(*,*)"Reading parameters from ",io_name
    call prg_open_file_to_read(io_unit,io_name)

    nsp = size(splist,dim=1)

    if(.not.allocated(ppot))then
      allocate(ppot(nsp,nsp))
      do i=1,nsp
        do j=1,nsp
          allocate(ppot(i,j)%potparams(16))
          ppot(i,j)%potparams = 0.0_dp
        enddo
      enddo
    endif

    read(io_unit,*) dummy, npots
    read(io_unit,*) dummy

    do i=1,npots
      read(io_unit,*)PairPotType1,PairPotType2,(PotCoef(j),j=1,10)
      do j=1,nsp
        do k=1,nsp
          if(PairPotType1.eq.spList(j))then
            if(PairPotType2.eq.spList(k))then
              ppot(j,k)%potparams(:) = PotCoef(:)
            endif
          endif
        enddo
      enddo
    enddo

    write(*,*)"Pair potentials:"

    do i=1,nsp
      do j=1,nsp
        if(ppot(i,j)%potparams(1).eq.0.0_dp)then
          ppot(i,j)%potparams(:) = ppot(j,i)%potparams(:)
        endif
      enddo
    enddo

    do i=1,nsp
      do j=1,nsp
        ! do PPID = 1:10
        PotCoef = ppot(i,j)%potparams(:)

        R1 = PotCoef(9);

        RCUT = PotCoef(10);

        R1SQ = R1*R1;

        POLY = R1*(PotCoef(2) + R1*(PotCoef(3) +  R1*(PotCoef(4) + R1*PotCoef(5))));
        SCL_R1 = PotCoef(1)*exp(POLY);
        EXPTMP = PotCoef(6)*exp(PotCoef(7)*(R1 - PotCoef(8)));
        PotCoef(11) = SCL_R1 + EXPTMP;
        DPOLY = PotCoef(2) + 2*PotCoef(3)*R1 + 3*PotCoef(4)*R1SQ + 4*PotCoef(5)*R1*R1SQ;
        PotCoef(12) = DPOLY*SCL_R1 + PotCoef(7)*EXPTMP;
        DDPOLY = 2*PotCoef(3) + 6*PotCoef(4)*R1 + 12*PotCoef(5)*R1SQ;
        PotCoef(13) = 0.5*((DPOLY*DPOLY + DDPOLY)*SCL_R1 + PotCoef(7)*PotCoef(7)*EXPTMP);

        !At the end of the join function:
        DELTA = RCUT - R1;
        DELTA2 = DELTA*DELTA;
        DELTA3 = DELTA2*DELTA;
        DELTA4 = DELTA3*DELTA;

        PotCoef(14) = (-1/DELTA3)*(3*PotCoef(13)*DELTA2 + 6*PotCoef(12)*DELTA + 10*PotCoef(11));
        PotCoef(15) = (1/DELTA4)*(3*PotCoef(13)*DELTA2 + 8*PotCoef(12)*DELTA + 15*PotCoef(11));

        PotCoef(16) = (-1/(10*DELTA3))*(6*PotCoef(15)*DELTA2 + 3*PotCoef(14)*DELTA + PotCoef(13));

        ppot(i,j)%potparams(:) = PotCoef

        write(*,'(2A2,16F10.5)')spList(i),spList(j),(ppot(i,j)%potparams(k),k=1,16)

      enddo
    enddo

  end subroutine load_PairPotTBparams


  subroutine scale_tail(a)
    implicit none
    integer                  ::  i
    real(dp)                 ::  ddpoly, prg_delta, prg_delta2, prg_delta3
    real(dp)                 ::  prg_delta4, dpoly, polynom, r1
    real(dp)                 ::  r1sq, rcut, rmod, scl_r1
    real(dp), intent(inout)  ::  a(:)

    if(abs(a(1)) < 1e-12)then
      do i=9,14
        a(i) = 0;
      enddo
    else
      r1 = a(7);
      rcut = a(8);
      r1sq = r1*r1;
      rmod = r1 - a(6);
      polynom = rmod*(a(2) + rmod*(a(3) + rmod*(a(4) + a(5)*rmod)));
      scl_r1 = exp(polynom);
      prg_delta = rcut - r1;
      ! now we're using a 6th order polynomial: fitted to value, first,
      ! and second derivatives at r1 and r_cut
      a(9) = scl_r1;
      rmod = r1 - a(6);
      dpoly = a(2) + 2*a(3)*rmod + 3*a(4)*rmod*rmod + 4*a(5)*rmod*rmod*rmod;
      a(10) = dpoly*scl_r1;
      ddpoly = 2*a(3) + 6*a(4)*rmod + 12*a(5)*rmod*rmod;
      a(11) = (dpoly*dpoly + ddpoly)*scl_r1/2;
      prg_delta2 = prg_delta*prg_delta;
      prg_delta3 = prg_delta2*prg_delta;
      prg_delta4 = prg_delta3*prg_delta;
      a(12) = (-1/prg_delta3)*(3*a(11)*prg_delta2 + 6*a(10)*prg_delta + 10*a(9));
      a(13) = (1/prg_delta4)*(3*a(11)*prg_delta2 + 8*a(10)*prg_delta + 15*a(9));
      a(14) = (-1/(10*prg_delta3))*(6*a(13)*prg_delta2 + 3*a(12)*prg_delta + a(11));
    endif
    write(*,'(16f10.5)')(a(i),i=1,16)

  end subroutine scale_tail

  subroutine write_latteTBparams(latteTBparams,filename)
    implicit none
    type(tbparams_type), intent(in) :: latteTBparams
    character(len=*), intent(in) :: filename
    character(100) :: io_name,tableformat
    character(2) :: dummy, tablesymbol
    character(4) :: tablesbasis
    integer :: io_unit,maxre,i,j,nsp

    io_name=trim(filename)
    write(*,*)"Writting TB parameters to ",io_name
    call prg_open_file(io_unit,io_name)
    write(io_unit,*)"Noelem= ",latteTBparams%nsp
    write(io_unit,'(13A10)')"Element","basis","Numel","Es","Ep","Ed",&
         "Ef","Mass","HubbardU","Wss","Wpp","Wdd","Wff"

    do i=1,latteTBparams%nsp
      tableformat='(A2,A4,11F10.5)'
      write(io_unit,tableformat)latteTBparams%splist(i),latteTBparams%basis(i)&
           ,latteTBparams%numel(i)&
           ,latteTBparams%onsite_energ(1,i)&
           ,latteTBparams%onsite_energ(2,i)&
           ,latteTBparams%onsite_energ(3,i)&
           ,latteTBparams%onsite_energ(4,i)&
           ,latteTBparams%mass(i) &
           ,latteTBparams%HubbardU(i)&
           ,latteTBparams%w(1,i)&
           ,latteTBparams%w(2,i)&
           ,latteTBparams%w(3,i)&
           ,latteTBparams%w(4,i)
    enddo

    close(io_unit)

  end subroutine write_latteTBparams

  subroutine write_bintTBparamsH(typeA,typeB,&
       intKind,intPairsH,intPairsS,filename)
    implicit none
    character(100)                                  ::  io_name, tableformat
    character(2), allocatable, intent(inout)        ::  TypeA(:,:), TypeB(:,:)
    character(3), allocatable                       ::  intKind(:)
    character(len=*), intent(in)                    ::  filename
    integer                                         ::  i, io_unit, j, k
    integer                                         ::  l, m, npoints, nsp
    type(intpairs_type), allocatable, intent(inout)  ::  intPairsH(:,:), intPairsS(:,:)

    io_name=trim(filename)
    write(*,*)""; write(*,*)"Writting parameters to ",io_name
    call prg_open_file(io_unit,io_name)
    nsp=size(typeA,dim=1)
    npoints = ((nsp*nsp + nsp)/2)*size(intKind,dim=1)
    write(io_unit,*)"Npoints ",npoints
    write(io_unit,'(19A10)')"Element1","Element2","Kind","H0","B1","B2","B3","B4",&
         "B5","R1","Rcut","S0","B1","B2","B3","B4","B5","R1","Rcut"
    tableformat='(A2,1X,A2,1X,A2,1X,A2,1X,A3,16F10.5)'
    do i=1,nsp
      do j=i,nsp
        do k = 1,size(intKind,dim=1)
          write(io_unit,tableformat)TypeA(i,j),"  ",TypeB(i,j),"  ",intKind(k),&
               (intPairsH(i,j)%intParams(m,k),m=1,8),(intPairsS(i,j)%intParams(m,k),m=1,8)
        enddo
      enddo
    enddo

  end subroutine write_bintTBparamsH

end module tbparams_latte_mod
