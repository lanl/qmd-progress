!> Kernel related routined.
!! \brief This module will be used to compute quantities related to the Kernel.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_kernel_mod

  use prg_kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_getKernel, gpmdcov_parseKernel

  !> General latte input variables type.
  !!
  type, public :: kernel_type

    !> Kernel type (if full or other)
    character(20) :: kernelType

    !> Verbosity level.
    integer :: verbose

    !> Threshold values for matrix elements.
    real(dp) :: threshold

    !> Use in mixing.
    logical :: kernelMixing

    !> Update each
    integer :: updateEach

  end type kernel_type

  type(kernel_type), public :: kernel

contains


  !> The parser for the kernel type.
  !!
  subroutine gpmdcov_parseKernel(kernel,filename)

    implicit none
    type(kernel_type) :: kernel
    integer, parameter :: nkey_char = 1, nkey_int = 2, nkey_re = 1, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'KernelType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'Full']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Verbose=', 'UpdateEach=']
    integer :: valvector_int(nkey_int) = (/ &
         0   ,     0  /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Threshold=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.00001 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'KernelMixing=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'KERNEL{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    kernel%kernelType = valvector_char(1)

    !Integers
    kernel%verbose = valvector_int(1)
    kernel%updateEach = valvector_int(2)

    !Reals
    kernel%threshold = valvector_re(1)

    !Logical
    kernel%kernelmixing = valvector_log(1)

  end subroutine gpmdcov_parseKernel


  subroutine gpmdcov_getKernel(nats)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: my_coul_pot(:),my_coul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_r(:,:), my_coul_pot_k(:), my_coul_pot_r(:)
    real(dp), allocatable :: Jacob(:), lwork(:), work(:)
    real(dp) :: mlsi
    integer :: nats, l, m, lm, info
    integer, allocatable :: ipiv(:)
    type(system_type), allocatable    :: mysyprt(:)

    if(.not. allocated(chargePertVect))allocate(chargePertVect(nats))

    allocate(mysyprt(gpat%TotalParts))

    allocate(my_coul_forces_k(3,nats))
    allocate(my_coul_forces_r(3,nats))
    allocate(my_coul_pot_k(nats))
    allocate(my_coul_pot_r(nats))

    allocate(Jacob(nats*nats))

    chargePertVect=0.0_dp
    do i=1,nats
      call gpmdcov_msI("gpmdcov_get_kernel","Constructing response for atom ="//to_string(i),lt%verbose,myRank)
      mlsi = mls()
      chargePertVect=0.0_dp
      chargePertVect(i)=1.0_dp
      call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,my_coul_forces_r,my_coul_pot_r);

      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,my_coul_forces_k,my_coul_pot_k);

 
        call gpmdcov_msIII("gpmdcov_get_kernel","Time for coulomb "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

#ifdef DO_MPI
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif
        norb = syprt(ipt)%estr%norbs

        if(.not.allocated(mysyprt(ipt)%estr%coul_pot_k))then
          allocate(mysyprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
          allocate(mysyprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))
          allocate(mysyprt(ipt)%net_charge(syprt(ipt)%nats))
        endif

        mysyprt(ipt)%estr%coul_pot_k = 0.0_dp
        mysyprt(ipt)%estr%coul_pot_r = 0.0_dp
        mysyprt(ipt)%net_charge = 0.0_dp

        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          mysyprt(ipt)%estr%coul_pot_k(j) = my_coul_pot_k(jj)
          mysyprt(ipt)%estr%coul_pot_r(j) = my_coul_pot_r(jj)
          mysyprt(ipt)%net_charge(j) = chargePertVect(jj)
        enddo
        
        !Get part(j)\%hamPert given mysyprt(j)\%over, and mysyprt(j)\%pot
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,mysyprt(ipt)%estr%ham)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,mysyprt(ipt)%estr%mat)
        call bml_copy_new(mysyprt(ipt)%estr%mat,mysyprt(ipt)%estr%ham)

        mlsi = mls()
        call gpmdcov_msI("gpmdcov_getKernel","Entering prg_get_hscf to construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(mysyprt(ipt)%estr%mat,syprt(ipt)%estr%over,mysyprt(ipt)%estr%ham,syprt(ipt)%spindex,&
             syprt(ipt)%estr%hindex,tb%hubbardu,mysyprt(ipt)%net_charge,&
             mysyprt(ipt)%estr%coul_pot_r,mysyprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)
        call gpmdcov_msIII("gpmdcov_get_kernel","Time for H construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call prg_orthogonalize(mysyprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,mysyprt(ipt)%estr%oham,&
             & lt%threshold,lt%bml_type,lt%verbose)

        call prg_toEigenspace(mysyprt(ipt)%estr%oham,mysyprt(ipt)%estr%ham,syprt(ipt)%estr%evects,lt%threshold,lt%verbose)
        !call bml_print_matrix("hamPert",mysyprt(ipt)%estr%ham,0,10,0,10)
        mlsi = mls()
        call prg_canon_response(mysyprt(ipt)%estr%mat,mysyprt(ipt)%estr%ham,nocc,beta,syprt(ipt)%estr%evects,&
             &syprt(ipt)%estr%aux(1,:),ef,12,norb)
        call gpmdcov_msIII("gpmdcov_get_kernel","Time for Canonincal Response construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        mlsi = mls()
        call prg_get_charges(mysyprt(ipt)%estr%mat, syprt(ipt)%estr%over, syprt(ipt)%estr%hindex, mysyprt(ipt)%net_charge, tb%numel,&
             syprt(ipt)%spindex, lt%mdim, lt%threshold)
        call gpmdcov_msIII("gpmdcov_get_kernel","Time for getting charges"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        !Constructing 1D matrix J to prepare for all-gather
        do j=1,gpat%sgraph(ipt)%llsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          Jacob(i*(nats - 1) + jj) = mysyprt(ipt)%net_charge(j)
        enddo

      enddo

    enddo

    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
      deallocate(mysyprt(ipt)%estr%coul_pot_k)
      deallocate(mysyprt(ipt)%estr%coul_pot_r)
      deallocate(mysyprt(ipt)%net_charge)
    enddo

    deallocate(mysyprt)
    deallocate(my_coul_pot_k)
    deallocate(my_coul_pot_r)
    deallocate(my_coul_forces_k)
    deallocate(my_coul_forces_r)

#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
      call prg_sumRealReduceN(Jacob, nats*nats)
    endif
#endif

    call prg_wait

    if(.not. allocated(Ker))then
      allocate(Ker(nats,nats))
    endif

    lm = 0
    do l=1,nats
      do m=1,nats
        lm = lm + 1
        Ker(l,m) = Jacob(lm)
        if(l == m) Ker(l,m) = Ker(l,m) - 1.0_dp
      enddo
    enddo

    deallocate(Jacob)

    allocate(work(nats+nats*nats))
    allocate(ipiv(nats))
    call DGETRF(nats, nats, Ker, nats, ipiv, info)
    call DGETRI(nats, Ker, nats, ipiv, work, nats+nats*nats, info)
    deallocate(work)
    deallocate(ipiv)

  end subroutine gpmdcov_getKernel

end module gpmdcov_kernel_mod
