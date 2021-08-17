!> Kernel Preconditioner related routined.
!! \brief This module will be used to compute quantities related to the
!! Kernel preconditioner.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_kernel_mod

  use prg_kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_getKernel, gpmdcov_getKernel_byBlocks, gpmdcov_parseKernel
  public :: gpmdcov_applyKernel, gpmdcov_getKernel_byParts

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


  !> Subroutine for computing the Kernel Preconditioner.
  !! \brief This particular routine builds the full kernel
  !! considering the system as a whole. For every perturbation
  !! on atom "i" partially fills column "i" of a full Jacobian
  !! (inverse of the Kernel). Paralelzation is obtained over
  !! computing the perturbation for the part.
  subroutine gpmdcov_getKernel(nats)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: my_coul_pot(:),my_coul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_r(:,:), my_coul_pot_k(:), my_coul_pot_r(:)
    real(dp), allocatable :: Jacob(:), lwork(:), work(:), ptcoul_pot_k(:),ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:), Q(:,:),H1(:,:),P1(:,:),row2(:)
    real(dp), allocatable :: mynumel(:)
    real(dp) :: mlsi, KSUM
    integer :: nats, l, m, lm, info
    integer, allocatable :: ipiv(:)
    type(system_type), allocatable    :: mysyprt(:)
    type(bml_matrix_t) :: zq_bml, myaux_bml, zqt_bml
    type(bml_matrix_t) :: ptham_bml, ptrho_bml, ptaux_bml

    !Charge perturbation vector. A charge +1 is added to an atomic
    !site to compute the response.
    if(.not. allocated(chargePertVect))allocate(chargePertVect(nats))

    !The coulomb potentials will be computed for every perturbation.
    allocate(my_coul_forces_k(3,nats))
    allocate(my_coul_forces_r(3,nats))
    allocate(my_coul_pot_k(nats))
    allocate(my_coul_pot_r(nats))

    allocate(Jacob(nats*nats))
    Jacob = 0.0_dp

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

      call gpmdcov_msI("gpmdcov_get_kernel","Time for coulomb "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

#ifdef DO_MPI
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif
        norb = syprt(ipt)%estr%norbs

        allocate(ptcoul_pot_k(syprt(ipt)%nats))
        allocate(ptcoul_pot_r(syprt(ipt)%nats))
        allocate(ptnet_charge(syprt(ipt)%nats))

        ptcoul_pot_k = 0.0_dp
        ptcoul_pot_r = 0.0_dp
        ptnet_charge = 0.0_dp

        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          ptcoul_pot_k(j) = my_coul_pot_k(jj)
          ptcoul_pot_r(j) = my_coul_pot_r(jj)
          ptnet_charge(j) = chargePertVect(jj)
        enddo

        !Get part(j)\%hamPert given mysyprt(j)\%over, and mysyprt(j)\%pot
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptham_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptrho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zq_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zqt_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptaux_bml)

        mlsi = mls()
        call gpmdcov_msI("gpmdcov_getKernel","Entering prg_get_hscf to construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(ptaux_bml,syprt(ipt)%estr%over,ptham_bml,syprt(ipt)%spindex,&
             syprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,lt%mdim,lt%threshold)
        call gpmdcov_msI("gpmdcov_get_kernel","Time for H construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)
        deallocate(ptcoul_pot_k)
        deallocate(ptcoul_pot_r)
        ! call prg_orthogonalize(mysyprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,mysyprt(ipt)%estr%oham,&
        !      & lt%threshold,lt%bml_type,lt%verbose)

        ! call prg_toEigenspace(mysyprt(ipt)%estr%oham,mysyprt(ipt)%estr%ham,syprt(ipt)%estr%evects,lt%threshold,lt%verbose)
        ! mlsi = mls()

        call bml_multiply(syprt(ipt)%estr%zmat,syprt(ipt)%estr%evects,zq_bml, 1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)

        !allocate(P1(norb,norb))
        !allocate(H1(norb,norb))
        !allocate(Q(norb,norb))
        !call bml_export_to_dense(ptham_bml,H1)
        !call bml_export_to_dense(syprt(ipt)%estr%evects,Q)
        !call Canon_DM_PRT(P1,H1,nocc,1000.0_dp,Q,syprt(ipt)%estr%evals,Ef,12,int(norb))
        !call bml_import_from_dense(lt%bml_type,P1,ptrho_bml,0.0_dp,norb)
        !deallocate(P1)
        !deallocate(H1)
        !deallocate(Q)

        call prg_canon_response(ptrho_bml,ptham_bml,nocc,23.208882712264995_dp,syprt(ipt)%estr%evects,&
             &syprt(ipt)%estr%evals,ef,12,norb)

        call gpmdcov_msI("gpmdcov_get_kernel","Time for Canonincal Response construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)


        call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)

        call bml_deallocate(ptham_bml)
        call bml_deallocate(zq_bml)
        call bml_deallocate(zqt_bml)
        call bml_deallocate(ptaux_bml)

        ! call bml_print_matrix("ptrho_bml",ptrho_bml,0,10,0,10)


        mlsi = mls()
        allocate(mynumel(10))
        mynumel = 0.0_dp
        call prg_get_charges(ptrho_bml, syprt(ipt)%estr%over, syprt(ipt)%estr%hindex, ptnet_charge, mynumel,&
             syprt(ipt)%spindex, lt%mdim, lt%threshold)
        deallocate(mynumel)
        !   do l = 1,nats
        !   write(*,*)"chrges",l,ptnet_charge(l)
        !   enddo

        call gpmdcov_msIII("gpmdcov_get_kernel","Time for getting charges"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)
        call bml_deallocate(ptrho_bml)

        !Constructing 1D matrix J to prepare for all-gather
        do j=1,gpat%sgraph(ipt)%llsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          Jacob((i-1)*nats + jj) = ptnet_charge(j)
        enddo
        deallocate(ptnet_charge)
      enddo

    enddo

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
        if(l == m) then
          Ker(m,l) = Jacob(lm) - 1.0_dp
        else
          Ker(m,l) = Jacob(lm)
        endif
        !  write(234,*)m,l,Ker(m,l)
      enddo
    enddo

    ! call bml_print_matrix("JJ",Ker,0,4,0,4)
    deallocate(Jacob)

    allocate(work(nats+nats*nats))
    allocate(ipiv(nats))

    call DGETRF(nats, nats, Ker, nats, ipiv, info)
    call DGETRI(nats, Ker, nats, ipiv, work, nats+nats*nats, info)
    deallocate(work)
    deallocate(ipiv)

    do m = 1,Nats
      KSum = 0.D0
      do l = 1,Nats
        KSum = KSum + Ker(l,m)
      enddo
      write(*,*) ' KSUM ', m, ' = ',KSum
    enddo

    !   call bml_print_matrix("Ker",Ker,0,4,0,4)
  end subroutine gpmdcov_getKernel


  subroutine gpmdcov_getKernel_byBlocks(nats)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: my_coul_pot(:),ptcoul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_k(:,:),my_coul_forces_r(:,:)
    real(dp), allocatable :: ptcoul_forces_r(:,:), my_coul_pot_k(:), my_coul_pot_r(:)
    real(dp), allocatable :: Jacob(:), lwork(:), work(:), ptcoul_pot_k(:),ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:), Q(:,:),H1(:,:),P1(:,:),row2(:)
    real(dp), allocatable :: mynumel(:)
    real(dp) :: mlsi, KSUM
    integer :: nats, l, m, lm, info, atom
    integer, allocatable :: ipiv(:)
    type(system_type), allocatable    :: mysyprt(:)
    type(bml_matrix_t) :: zq_bml, myaux_bml, zqt_bml
    type(bml_matrix_t) :: ptham_bml, ptrho_bml, ptaux_bml

    call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Entering routine",myRank)

    !Each rank will do its own part
#ifdef DO_MPI
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    do ipt = 1,gpat%TotalParts
#endif
      !Charge perturbation vector. A charge +1 is added to an atomic
      !site to compute the response.
      if(.not. allocated(chargePertVect))allocate(chargePertVect(sy%nats))

      !The coulomb potentials will be computed for every perturbation.
      if(.not.allocated(my_coul_forces_k))allocate(my_coul_forces_k(3,sy%nats))
      if(.not.allocated(my_coul_forces_r))allocate(my_coul_forces_r(3,sy%nats))
      if(.not.allocated(my_coul_pot_k))allocate(my_coul_pot_k(sy%nats))
      if(.not.allocated(my_coul_pot_r))allocate(my_coul_pot_r(sy%nats))
      my_coul_pot_k = 0.0_dp
      my_coul_pot_r = 0.0_dp

      if(.not. allocated(Jacob)) allocate(Jacob(sy%nats*sy%nats))
      Jacob = 0.0_dp

      chargePertVect=0.0_dp

      do i=1,syprt(ipt)%nats
        mlsi = mls()
        chargePertVect=0.0_dp
        atom = gpat%sgraph(ipt)%core_halo_index(i)+1
        chargePertVect(atom)=1.0_dp

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Constructing response for atom ="//to_string(atom),lt%verbose,myRank)

        call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
             ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
             sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
             nl%nnIz,nl%nrnnlist,nl%nnType,my_coul_forces_r,my_coul_pot_r);

        call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
             ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
             sy%recip_vector,sy%volr,lt%coul_acc,my_coul_forces_k,my_coul_pot_k);

        if(.not.allocated(ptcoul_pot_k))allocate(ptcoul_pot_k(syprt(ipt)%nats))
        if(.not.allocated(ptcoul_pot_r))allocate(ptcoul_pot_r(syprt(ipt)%nats))
        if(.not.allocated(ptnet_charge))allocate(ptnet_charge(syprt(ipt)%nats))

        ptcoul_pot_k = 0.0_dp
        ptcoul_pot_r = 0.0_dp
        ptnet_charge = 0.0_dp

        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          ptcoul_pot_k(j) = my_coul_pot_k(jj)
          ptcoul_pot_r(j) = my_coul_pot_r(jj)
          ptnet_charge(j) = chargePertVect(jj)
        enddo

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for coulomb "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        norb = syprt(ipt)%estr%norbs

        !Get part(j)\%hamPert given mysyprt(j)\%over, and mysyprt(j)\%pot
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptham_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptrho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zq_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zqt_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptaux_bml)

        mlsi = mls()
        call gpmdcov_msI("gpmdcov_getKernel_byBlocks","Entering prg_get_hscf to construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(ptaux_bml,syprt(ipt)%estr%over,ptham_bml,syprt(ipt)%spindex,&
             syprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,lt%mdim,lt%threshold)
        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for H construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)
        deallocate(ptcoul_pot_k)
        deallocate(ptcoul_pot_r)

        call bml_multiply(syprt(ipt)%estr%zmat,syprt(ipt)%estr%evects,zq_bml, 1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)

        call prg_canon_response(ptrho_bml,ptham_bml,nocc,23.208882712264995_dp,syprt(ipt)%estr%evects,&
             &syprt(ipt)%estr%evals,ef,12,norb)

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for Canonincal Response construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)

        call bml_deallocate(ptham_bml)
        call bml_deallocate(zq_bml)
        call bml_deallocate(zqt_bml)
        call bml_deallocate(ptaux_bml)

        mlsi = mls()
        allocate(mynumel(10))
        mynumel = 0.0_dp
        call prg_get_charges(ptrho_bml, syprt(ipt)%estr%over, syprt(ipt)%estr%hindex, ptnet_charge, mynumel,&
             syprt(ipt)%spindex, lt%mdim, lt%threshold)
        deallocate(mynumel)

        call gpmdcov_msIII("gpmdcov_get_kernel_byBlocks","Time for getting charges"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)
        call bml_deallocate(ptrho_bml)

        !Constructing 1D matrix J to prepare for all-gather
        do j=1,gpat%sgraph(ipt)%llsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          Jacob((atom-1)*nats + jj) = ptnet_charge(j)
        enddo
        deallocate(ptnet_charge)
      enddo

    enddo

    if(allocated(ptcoul_pot_k))deallocate(ptcoul_pot_k)
    if(allocated(ptcoul_pot_r))deallocate(ptcoul_pot_r)
    if(allocated(ptcoul_forces_k))deallocate(ptcoul_forces_k)

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
        if(l == m) then
          Ker(m,l) = Jacob(lm) - 1.0_dp
        else
          Ker(m,l) = Jacob(lm)
        endif
      enddo
    enddo

    deallocate(Jacob)

    allocate(work(nats+nats*nats))
    allocate(ipiv(nats))

    call DGETRF(nats, nats, Ker, nats, ipiv, info)
    call DGETRI(nats, Ker, nats, ipiv, work, nats+nats*nats, info)
    deallocate(work)
    deallocate(ipiv)

    do m = 1,Nats
      KSum = 0.D0
      do l = 1,Nats
        KSum = KSum + Ker(l,m)
      enddo
      write(*,*) ' KSUM ', m, ' = ',KSum
    enddo

  end subroutine gpmdcov_getKernel_byBlocks

  subroutine gpmdcov_getKernel_byParts(nats,mysyprt)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: my_coul_pot(:),ptcoul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_k(:,:),my_coul_forces_r(:,:)
    real(dp), allocatable :: ptcoul_forces_r(:,:), my_coul_pot_k(:)
    real(dp), allocatable :: my_coul_pot_r(:)
    real(dp), allocatable :: Jacob(:), lwork(:), work(:)
    real(dp), allocatable :: ptcoul_pot_k(:),ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:), Q(:,:),H1(:,:),P1(:,:),row2(:)
    real(dp) :: mynumel(10)
    real(dp) :: mlsi, KSUM
    integer :: nats, l, m, lm, info, atom, coreSize
    integer, allocatable :: ipiv(:)
    type(system_type), allocatable, intent(inout) :: mysyprt(:)
    type(bml_matrix_t) :: zq_bml, myaux_bml, zqt_bml
    type(bml_matrix_t) :: ptham_bml, ptrho_bml, ptaux_bml

    !call gpmdcov_msI("gpmdcov_get_kernel_byParts","Entering routine",myRank)
    !Each rank will do its own part and keep it.
    if(.not.allocated(my_coul_forces_k))allocate(my_coul_forces_k(3,sy%nats))
    if(.not.allocated(my_coul_forces_r))allocate(my_coul_forces_r(3,sy%nats))
    if(.not.allocated(my_coul_pot_k))allocate(my_coul_pot_k(sy%nats))
    if(.not.allocated(my_coul_pot_r))allocate(my_coul_pot_r(sy%nats))

#ifdef DO_MPI
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    do ipt = 1,gpat%TotalParts
#endif
      !Charge perturbation vector. A charge +1 is added to an atomic
      !site to compute the response.
      if(.not.allocated(chargePertVect))allocate(chargePertVect(sy%nats))

      !The coulomb potentials will be computed for every perturbation.
      my_coul_pot_k = 0.0_dp
      my_coul_pot_r = 0.0_dp
      my_coul_forces_k = 0.0_dp
      my_coul_forces_r = 0.0_dp

      coreSize = gpat%sgraph(ipt)%llsize

      if(allocated(mysyprt(ipt)%estr%ker)) deallocate(mysyprt(ipt)%estr%ker)
      allocate(mysyprt(ipt)%estr%ker(coreSize,coreSize))

      mysyprt(ipt)%estr%ker = 0.0_dp

      do i=1,coreSize
        mlsi = mls()
        chargePertVect=0.0_dp
        atom = gpat%sgraph(ipt)%core_halo_index(i)+1
        chargePertVect(atom)=1.0_dp

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Constructing response&
             &for atom ="//to_string(atom),lt%verbose,myRank)

        call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
             ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
             sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
             nl%nnIz,nl%nrnnlist,nl%nnType,my_coul_forces_r,my_coul_pot_r);

        call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
             ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
             sy%recip_vector,sy%volr,lt%coul_acc,my_coul_forces_k,my_coul_pot_k);

        allocate(ptcoul_pot_k(mysyprt(ipt)%nats))
        allocate(ptcoul_pot_r(mysyprt(ipt)%nats))
        allocate(ptnet_charge(mysyprt(ipt)%nats))

        ptcoul_pot_k = 0.0_dp
        ptcoul_pot_r = 0.0_dp
        ptnet_charge = 0.0_dp

        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          ptcoul_pot_k(j) = my_coul_pot_k(jj)
          ptcoul_pot_r(j) = my_coul_pot_r(jj)
          ptnet_charge(j) = chargePertVect(jj)
        enddo

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for coulomb&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        norb = mysyprt(ipt)%estr%norbs

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptham_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptrho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zq_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zqt_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ptaux_bml)

        mlsi = mls()
        call gpmdcov_msI("gpmdcov_getKernel_byBlocks","Entering prg_get_hscf to&
             &construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(ptaux_bml,mysyprt(ipt)%estr%over,ptham_bml,mysyprt(ipt)%spindex,&
             mysyprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,lt%mdim,lt%threshold)
        deallocate(ptcoul_pot_r)
        deallocate(ptcoul_pot_k)

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for H construction&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_multiply(mysyprt(ipt)%estr%zmat,mysyprt(ipt)%estr%evects,zq_bml,&
             &1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)

        call prg_canon_response(ptrho_bml,ptham_bml,nocc,23.208882712264995_dp,mysyprt(ipt)%estr%evects,&
             &mysyprt(ipt)%estr%evals,ef,12,norb)

        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for Canonincal&
             &Response construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)

        call bml_deallocate(ptham_bml)
        call bml_deallocate(zq_bml)
        call bml_deallocate(zqt_bml)
        call bml_deallocate(ptaux_bml)

        mlsi = mls()
        mynumel = 0.0_dp

        call prg_get_charges(ptrho_bml, mysyprt(ipt)%estr%over,&
             &mysyprt(ipt)%estr%hindex, ptnet_charge, mynumel,&
             mysyprt(ipt)%spindex, norb, lt%threshold)
        call bml_deallocate(ptrho_bml)
        call gpmdcov_msIII("gpmdcov_get_kernel_byBlocks","Time for getting &
             &charges"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        !Constructing J to prepare for K
        do j=1,gpat%sgraph(ipt)%llsize
          mysyprt(ipt)%estr%ker(j,i) = ptnet_charge(j)
          if(i == j)mysyprt(ipt)%estr%ker(j,i) = mysyprt(ipt)%estr%ker(j,i) - 1.0_dp
        enddo
        deallocate(ptnet_charge)
        call bml_deallocate(ptrho_bml)
      enddo

      if(allocated(work))deallocate(work);allocate(work(coreSize+coreSize*coreSize))
      if(allocated(ipiv))deallocate(ipiv);allocate(ipiv(coreSize))
      call DGETRF(coreSize, coreSize, mysyprt(ipt)%estr%ker, coreSize, ipiv, info)
      call DGETRI(coreSize, mysyprt(ipt)%estr%ker, coreSize, ipiv, work, coreSize+coreSize*coreSize, info)
      deallocate(work)
      deallocate(ipiv)

      do m = 1,coreSize
        KSum = 0.0_dp
        do l = 1,coreSize
          KSum = KSum + mysyprt(ipt)%estr%ker(l,m)
        enddo
        write(*,*) ' KSUM ', m, ' = ',KSum
      enddo

    enddo
    deallocate(my_coul_forces_k)
    deallocate(my_coul_forces_r)
    deallocate(my_coul_pot_k)
    deallocate(my_coul_pot_r)

  end subroutine gpmdcov_getKernel_byParts


  subroutine gpmdcov_applyKernel(my_nguess,my_chargesOld,kernelTimesRes)
    use gpmdcov_vars
    implicit none
    real(dp), allocatable, intent(inout) :: kernelTimesRes(:)
    real(dp), allocatable, intent(in) :: my_nguess(:)
    real(dp), allocatable, intent(in) ::  my_chargesOld(:)
    real(dp), allocatable :: nguessPart(:), chargesOldPart(:)
    real(dp), allocatable :: kernelTimesResPart(:)
    integer :: nats, myj, myjj

      !> Get Coulombic potential and charges for the part.
      nats = size(my_nguess,dim=1)
      if(.not. allocated(kernelTimesRes)) allocate(kernelTimesRes(nats))

#ifdef DO_MPI
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif
        allocate(kernelTimesResPart(gpat%sgraph(ipt)%llsize))
        allocate(nguessPart(gpat%sgraph(ipt)%llsize))
        allocate(chargesOldPart(gpat%sgraph(ipt)%llsize))

        !Collapse old charges and nguess for the part
        do myj=1,gpat%sgraph(ipt)%llsize
          myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
          chargesOldPart(myj) = my_chargesOld(myjj)
          nguessPart(myj) = my_nguess(myjj)
        enddo
        kernelTimesResPart = MATMUL(syprt(ipt)%estr%ker,(nguessPart-chargesOldPart))

        !Expand old charges and nguess for the part
        do myj=1,gpat%sgraph(ipt)%llsize
          myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
          kernelTimesRes(myjj) = kernelTimesResPart(myj)
        enddo
        
        deallocate(kernelTimesResPart)
        deallocate(nguessPart)
        deallocate(chargesOldPart)

      enddo

#ifdef DO_MPI
      if (getNRanks() .gt. 1) then
        call prg_sumRealReduceN(kernelTimesRes, nats)
      endif
#endif

      call prg_wait

  end subroutine gpmdcov_applyKernel




end module gpmdcov_kernel_mod
