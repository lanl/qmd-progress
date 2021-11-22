!> Kernel Preconditioner related routined.
!! \brief This module will be used to compute quantities related to the
!! Kernel preconditioner.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_kernel_mod

  use prg_kernelparser_mod
  use gpmdcov_allocation_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_getKernel, gpmdcov_getKernel_byBlocks, gpmdcov_parseKernel
  public :: gpmdcov_applyKernel, gpmdcov_getKernel_byParts, gpmdcov_rankN_update_byParts

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

    !> If the full preconditioner needs to be always built
    logical :: buildAlways

    !> Number of rank N updates
    integer :: rankNUpdate

  end type kernel_type

  type(kernel_type), public :: kernel

contains


  !> The parser for the kernel type.
  !!
  subroutine gpmdcov_parseKernel(kernel,filename)

    implicit none
    type(kernel_type) :: kernel
    integer, parameter :: nkey_char = 1, nkey_int = 3, nkey_re = 1, nkey_log = 2
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'KernelType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'Full']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Verbose=', 'UpdateEach=', 'RankNUpdate=']
    integer :: valvector_int(nkey_int) = (/ &
         0   ,     0,    0  /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Threshold=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.00001 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'KernelMixing=','BuildAlways=']
    logical :: valvector_log(nkey_log) = (/&
         .false., .false./)

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
    kernel%rankNUpdate = valvector_int(3)

    !Reals
    kernel%threshold = valvector_re(1)

    !Logical
    kernel%kernelmixing = valvector_log(1)
    kernel%buildAlways = valvector_log(2)

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
    integer :: nats, l, m, lm, info, norbs
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
        norbs = syprt(ipt)%estr%norbs

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
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptham_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptrho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,zq_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,zqt_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptaux_bml)

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

        !allocate(P1(norbs,norbs))
        !allocate(H1(norbs,norbs))
        !allocate(Q(norbs,norbs))
        !call bml_export_to_dense(ptham_bml,H1)
        !call bml_export_to_dense(syprt(ipt)%estr%evects,Q)
        !call Canon_DM_PRT(P1,H1,nocc,1000.0_dp,Q,syprt(ipt)%estr%evals,Ef,12,int(norbs))
        !call bml_import_from_dense(lt%bml_type,P1,ptrho_bml,0.0_dp,norbs)
        !deallocate(P1)
        !deallocate(H1)
        !deallocate(Q)

        !call prg_canon_response(ptrho_bml,ptham_bml,norbs,23.208882712264995_dp,syprt(ipt)%estr%evects,&
        !    &syprt(ipt)%estr%evals,ef,12,norbs)
        call prg_canon_response(ptrho_bml,ptham_bml,norbs,beta,syprt(ipt)%estr%evects,&
             &syprt(ipt)%estr%evals,ef,12,norbs)
        

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
    integer :: nats, l, m, lm, info, atom, norbs
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

        norbs = syprt(ipt)%estr%norbs

        !Get part(j)\%hamPert given mysyprt(j)\%over, and mysyprt(j)\%pot
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptham_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptrho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,zq_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,zqt_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptaux_bml)

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

        call prg_canon_response(ptrho_bml,ptham_bml,norbs,23.208882712264995_dp,syprt(ipt)%estr%evects,&
             &syprt(ipt)%estr%evals,ef,12,norbs)

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
    real(dp), allocatable :: myOverCore(:,:),myOver(:,:),myPtRho(:,:),myPtRhoCore(:,:)
    integer, allocatable :: myHindex(:,:),myHindexCore(:,:),mySpindex(:),mySpindexCore(:)
    real(dp) :: mynumel(10)
    real(dp) :: mlsi, KSUM, trdPdMuAO, trp1, mu1
    integer :: nats, l, m, lm, info, atom, coreSize, norbsCore, norbs
    integer, allocatable :: ipiv(:)
    type(system_type), allocatable, intent(inout) :: mysyprt(:)
    type(bml_matrix_t) :: zq_bml, myaux_bml, zqt_bml
    type(bml_matrix_t) :: ptham_bml, ptrho_bml, ptaux_bml, myptrhoCore_bml
    type(bml_matrix_t) :: myOverCore_bml,dPdMuAO_bml,p1_bml,dPdMuAOS_bml,p1S_bml
    real(dp), allocatable :: dPdMuAO_dia(:),p1_dia(:),dPdMu(:)


    if(.not.allocated(my_coul_forces_k))allocate(my_coul_forces_k(3,sy%nats))
    if(.not.allocated(my_coul_forces_r))allocate(my_coul_forces_r(3,sy%nats))
    if(.not.allocated(my_coul_pot_k))allocate(my_coul_pot_k(sy%nats))
    if(.not.allocated(my_coul_pot_r))allocate(my_coul_pot_r(sy%nats))

    !Each rank will do its own part and keep it.
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
      !allocate(mysyprt(ipt)%estr%ker(gpat%sgraph(ipt)%lsize,gpat%sgraph(ipt)%lsize))

      mysyprt(ipt)%estr%ker = 0.0_dp

      do i=1, coreSize
        mlsi = mls()
        chargePertVect=0.0_dp
        atom = gpat%sgraph(ipt)%core_halo_index(i)+1
        chargePertVect(atom)=1.0_dp

        call gpmdcov_msI("gpmdcov_get_kernel_byParts","Constructing response&
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

        !Now only for core
        !allocate(ptcoul_pot_k(coreSize))
        !allocate(ptcoul_pot_r(coreSize))
        !allocate(ptnet_charge(coreSize))

        ptcoul_pot_k = 0.0_dp
        ptcoul_pot_r = 0.0_dp
        ptnet_charge = 0.0_dp

        !> Get Coulombic potential and charges for the core part.
        !do j=1,coreSize
        do j=1,mysyprt(ipt)%nats
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          ptcoul_pot_k(j) = my_coul_pot_k(jj)
          ptcoul_pot_r(j) = my_coul_pot_r(jj)
          ptnet_charge(j) = chargePertVect(jj)
        enddo

        call gpmdcov_msI("gpmdcov_get_kernel_byParts","Time for coulomb&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        norbs = mysyprt(ipt)%estr%norbs
        norbsCore = mysyprt(ipt)%estr%norbsCore

!!!   System -> Cores+halos(0)

        !H =   [H_c      ]
        !      [     H_h ]


        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptham_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptrho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,zq_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,zqt_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,ptaux_bml)

        mlsi = mls()
        call gpmdcov_msI("gpmdcov_getKernel_byParts","Entering prg_get_hscf to&
             &construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(ptaux_bml,mysyprt(ipt)%estr%over,ptham_bml,mysyprt(ipt)%spindex,&
             mysyprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,lt%mdim,lt%threshold)


        deallocate(ptcoul_pot_r)
        deallocate(ptcoul_pot_k)

        call gpmdcov_msI("gpmdcov_get_kernel_byParts","Time for H construction&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_multiply(mysyprt(ipt)%estr%zmat,mysyprt(ipt)%estr%evects,zq_bml,&
             &1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)

        allocate(dPdMu(norbs))
        !call prg_canon_response_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
        !     &norbsCore,23.208882712264995_dp,mysyprt(ipt)%estr%evects,&
        !     &mysyprt(ipt)%estr%evals,ef,12,norbs)
        call prg_canon_response_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
             &norbsCore,beta,mysyprt(ipt)%estr%evects,&
             &mysyprt(ipt)%estr%evals,ef,12,norbs)
        call bml_get_diagonal(p1_bml,p1_dia)
        trP1 = sum(p1_dia(1:norbsCore))
        write(*,*)"trP1 (After canon)",trP1
        call gpmdcov_msI("gpmdcov_get_kernel_byParts","Time for Canonincal&
             &Response construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,dPdMuAO_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,dPdMuAOS_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,p1S_bml)
        call bml_set_diagonal(dPdMuAO_bml,dPdMu)

        deallocate(dPdMu)

        call bml_multiply(zq_bml,p1_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,p1_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(p1_bml,syprt(ipt)%estr%over,p1S_bml,1.0_dp,0.0_dp,0.0_dp)

        call bml_multiply(zq_bml,dPdMuAO_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,dPdMuAO_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(dPdMuAO_bml,syprt(ipt)%estr%over,dPdMuAOS_bml,1.0_dp,0.0_dp,0.0_dp)

        call bml_get_diagonal(dPdMuAOS_bml,dPdMuAO_dia)
        call bml_get_diagonal(p1S_bml,p1_dia)
        ! write(*,*)"p1_dia",p1_dia
        ! write(*,*)"dPdMu_dia",dPdMuAO_dia
        trP1 = sum(p1_dia(1:norbsCore))
        write(*,*)"trP1",trP1
        trdPdMuAO = sum(dPdMuAO_dia(1:norbsCore))
        deallocate(dPdMuAO_dia)
        deallocate(p1_dia)
        if(abs(trdPdMuAO) < 1.0E-12)then
              mu1 = 0.0
        else
              mu1 =  -trP1/trdPdMuAO
        endif
        write(*,*)"trP1,trdPdMuAO,mu1",trP1,trdPdMuAO,mu1
        write(*,*)"sizes",size(p1_dia,dim=1),size(dPdMuAO_dia,dim=1),norbs
        !stop
        call bml_copy(p1_bml,ptrho_bml)
        !call bml_add_deprecated(2.0_dp,ptrho_bml,mu1,dPdMuAO_bml,lt%threshold)
        call bml_add(ptrho_bml,dPdMuAO_bml,2.0_dp,2.0_dp*mu1,lt%threshold)
        !call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        !call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)

        call bml_deallocate(ptham_bml)
        call bml_deallocate(zq_bml)
        call bml_deallocate(zqt_bml)
        call bml_deallocate(ptaux_bml)
        call bml_deallocate(p1_bml)
        call bml_deallocate(dPdMuAO_bml)

        mlsi = mls()
        mynumel = 0.0_dp

        !! P1_MO = P1_MO + mu_1*dPdmu -> P1_ao, Tr_core(P1_ao*S) = 0
        !! P1_MO = P1_MO + mu_1*dPdmu -> P1_ao, Tr_core(ZQ*(P1_e + mu_1 dPdMu_e)QZ'*S) = 0
        !! P1_MO = P1_MO + mu_1*dPdmu -> P1_ao, Tr_core(ZQ*P1_e*ZQ'*S) + mu_1*Tr_core(ZQ*dPdMu_e*ZQ'*S) = 0
        !! P1_ao = ZQ*P1_e*ZQ' + mu_1*ZQ*dPdMu_e*QZ'
        call prg_get_charges(ptrho_bml, mysyprt(ipt)%estr%over,&
             &mysyprt(ipt)%estr%hindex, ptnet_charge, mynumel,&
             mysyprt(ipt)%spindex, norbs, lt%threshold)

        call bml_deallocate(ptrho_bml)
        call gpmdcov_msIII("gpmdcov_get_kernel_byParts","Time for getting &
             &charges"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

      do m = 1,coreSize
        KSum = 0.0_dp
        do l = 1,coreSize
          KSum = KSum + ptnet_charge(l)
        enddo
        write(*,*) ' JSUM ', m, ' = ',KSum
      enddo           
        !write(*,*)"Sum net charge",sum(ptnet_charge(1:coreSize))
        write(*,*)"Sum net charge",sum(ptnet_charge)
        !coreSize = gpat%sgraph(ipt)%lsize
        !Constructing J to prepare for K
        !do j=1,gpat%sgraph(ipt)%lsize
        do j=1,coreSize
          mysyprt(ipt)%estr%ker(j,i) = ptnet_charge(j)
          if(i == j)mysyprt(ipt)%estr%ker(j,i) = mysyprt(ipt)%estr%ker(j,i) - 1.0_dp
        enddo
        deallocate(ptnet_charge)
        call bml_deallocate(ptrho_bml)

      enddo

      if(allocated(work))deallocate(work);allocate(work(coreSize+coreSize*coreSize))
      if(allocated(ipiv))deallocate(ipiv);allocate(ipiv(coreSize))
      call DGETRF(coreSize, coreSize, mysyprt(ipt)%estr%ker, coreSize, ipiv, info)
      call DGETRI(coreSize, mysyprt(ipt)%estr%ker, coreSize, ipiv,&
           & work, coreSize+coreSize*coreSize, info)
      deallocate(work)
      deallocate(ipiv)

      do m = 1,coreSize
        KSum = 0.0_dp
        do l = 1,coreSize
          KSum = KSum + mysyprt(ipt)%estr%ker(l,m)
        enddo
        write(*,*) ' KSUM ', m, ' = ',KSum
      enddo
! Hack to use the identity (to see the effect of low rank updates)
!      mysyprt(ipt)%estr%ker = 0.0_dp
!      do m = 1,coreSize
!        mysyprt(ipt)%estr%ker(m,m) = -0.5_dp
!      enddo
    enddo


    write(*,*)"MDBUG, alloc ker",allocated(mysyprt(ipt)%estr%ker),coreSize,size(mysyprt(ipt)%estr%ker)
    !stop
    deallocate(my_coul_forces_k)
    deallocate(my_coul_forces_r)
    deallocate(my_coul_pot_k)
    deallocate(my_coul_pot_r)
  end subroutine gpmdcov_getKernel_byParts


  subroutine gpmdcov_rankN_update_byParts(myqn,myn,mysyprt,maxRanks,KK0Res)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    use gpmdcov_allocation_mod
    implicit none
    integer :: myi,mii,myj,myjj,irank,mRanks
    integer :: norbsCore,norbs
    integer, intent(in) :: maxRanks
    integer :: natsCore,natsCoreHalo
    integer, allocatable :: maxCoresAmongParts(:)
    integer :: maxCoresAmongPartsAndRanks
    real(dp), allocatable, intent(in) :: myqn(:), myn(:)
    real(dp), allocatable :: myqnPart(:), mynPart(:)
    real(dp), allocatable :: dr(:),K0Res(:),K0ResPart(:,:),dr_save(:,:)
    real(dp), allocatable :: vi(:,:), v(:)
    real(dp), allocatable :: my_coul_pot(:),ptcoul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_k(:,:),my_coul_forces_r(:,:)
    real(dp), allocatable :: ptcoul_forces_r(:,:), my_coul_pot_k(:)
    real(dp), allocatable :: my_coul_pot_r(:),chargePertVect(:)
    real(dp), allocatable :: ptcoul_pot_k(:)
    real(dp), allocatable :: ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:)
    real(dp), allocatable :: v_core_i(:,:,:)
    real(dp), allocatable :: dPdMu(:), IdK0Res(:)
    real(dp), allocatable :: q1(:,:),dqdmu(:,:)
    real(dp), allocatable :: f(:),c_i(:),ff(:,:,:)
    real(dp), allocatable :: oij(:,:),mMat(:,:)
    real(dp), allocatable :: vRank(:)
    real(dp), intent(inout), allocatable :: KK0Res(:)
    real(dp) :: trP1(1), trdPdMu(1)
    real(dp) :: mu1_global
    type(bml_matrix_t) :: ptham_bml, ptaux_bml
    type(bml_matrix_t) :: zq_bml,zqt_bml, ptrho_bml
    type(bml_matrix_t) :: p1_bml,dPdMu_bml
    type(system_type), allocatable, intent(inout) :: mysyprt(:)
    real(dp) :: mynumel(10)
    real(dp), allocatable :: lwork(:),work(:),auxVect(:)
    integer :: info
    integer, allocatable :: ipiv(:)
   
    call gpmdcov_msI("gpmdcov_rankN_update_byParts","Updating the Kernel",lt%verbose,myRank)
 
    if(.not.allocated(K0Res))allocate(K0Res(sy%nats))

    K0Res = 0.0_dp
 
    !This will be used to capture the maximun number of core atoms 
    !among all the parts/subsystems for the present MPI rank.
    allocate(maxCoresAmongParts(nRanks))
    maxCoresAmongParts = 0
#ifdef DO_MPI
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    do iptt = 1,gpat%TotalParts
      ipt = iptt
#endif
      natsCore = gpat%sgraph(ipt)%llsize
      maxCoresAmongParts(myRank) = max(maxCoresAmongParts(myRank),natsCore)
    enddo
    write(*,*)"myn",myn(1:5)
    write(*,*)"myqn",myqn(1:5)
!stop
    call gpmdcov_reallocate_realMat(K0ResPart,maxCoresAmongParts(myRank),partsInEachRank(myRank))

    K0ResPart = 0.0_dp

    !In this loop we will compute K0Res which is the product of the
    !Preconditioner K with the residue q(n) - n
#ifdef DO_MPI
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    do iptt = 1,gpat%TotalParts
      ipt = iptt
#endif
   
      natsCore = gpat%sgraph(ipt)%llsize

     ! maxCoresAmongParts(myRank) = max(maxCoresAmongParts(myRank),natsCore)

      call gpmdcov_reallocate_realVect(myqnPart,natsCore)
      call gpmdcov_reallocate_realVect(mynPart,natsCore)

      !Get old charges and nguess for the part
      do myj=1,natsCore
        myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
        mynPart(myj) = myn(myjj)
        myqnPart(myj) = myqn(myjj)
      enddo

      write(*,*)"MDBUG",size(myqnPart),size(mynPart),size(mysyprt(ipt)%estr%ker)
      K0ResPart(1:natsCore,iptt) = MATMUL(mysyprt(ipt)%estr%ker,(myqnPart-mynPart))
      !write(*,*)mysyprt(ipt)%estr%ker(1,1:5),mysyprt(ipt)%estr%ker(norbs,norbs-5:norbs)
      !write(*,*)myqnPart(1:5),mynPart(1:5)
      !write(*,*)"K0ResPart15-1-1",K0ResPart(1:5,1),iptt
      !write(*,*)"K0ResPart15-1-2",K0ResPart(1:5,2),iptt
      !write(*,*)"IPTT",iptt

      !Expand K0resPart into K0Res
      do myj=1,natsCore
        myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
        K0Res(myjj) = K0ResPart(myj,iptt)
      enddo
      deallocate(myqnPart)
      deallocate(mynPart)

    enddo
#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
      call prg_sumRealReduceN(K0Res,sy%nats)
      !call prg_sumIntReduceN(maxCoresAmongParts,sy%nats)
      call prg_sumIntReduceN(maxCoresAmongParts,nRanks)
    endif
    call prg_wait()
#endif

    !This is the maximun core atoms among "all" of the parts/subsystems
    maxCoresAmongPartsAndRanks = maxval(maxCoresAmongParts)

    allocate(dr(sy%nats))

    if(.not.allocated(my_coul_forces_k))allocate(my_coul_forces_k(3,sy%nats))
    if(.not.allocated(my_coul_forces_r))allocate(my_coul_forces_r(3,sy%nats))
    if(.not.allocated(my_coul_pot_k))allocate(my_coul_pot_k(sy%nats))
    if(.not.allocated(my_coul_pot_r))allocate(my_coul_pot_r(sy%nats))
    my_coul_pot_k = 0.0_dp
    my_coul_pot_r = 0.0_dp
    my_coul_forces_k = 0.0_dp
    my_coul_forces_r = 0.0_dp

    dr = K0Res

    write(*,*)"K0Res"
    write(*,*)K0Res(1:5)
    write(*,*)K0Res(sy%nats - 5:sy%nats)

    mRanks = maxRanks !Number of total rank updates 
    allocate(vi(sy%nats,mRanks))
    allocate(dr_save(sy%nats,mranks))
    vi = 0.0_dp
    allocate(v_core_i(maxCoresAmongPartsAndRanks,partsInEachRank(myRank),mRanks))
    v_core_i = 0.0_dp
    allocate(c_i(mRanks))
    c_i = 0.0_dp
    allocate(ff(maxCoresAmongPartsAndRanks,partsInEachRank(myRank),mRanks))
    ff = 0.0_dp
    
    !Here we enter the loop for the rank updates (do not confuse with MPI rank)
    do irank = 1, mRanks
      vi(:,irank) = dr/norm2(dr)
      !Gram-Schmidt OrthoNormalization
      if(irank > 1)then
        do kk = 1,irank-1
          vi(:,irank) = vi(:,irank) - DOT_PRODUCT(vi(:,kk),vi(:,irank))*vi(:,kk) 
        enddo
        vi(:,irank) = vi(:,irank)/norm2(vi(:,irank))
      endif
      chargePertVect = vi(:,irank)
      write(*,*)"chargePertVect",chargepertVect(1:5)
  
      call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,my_coul_forces_r,my_coul_pot_r);

      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,my_coul_forces_k,my_coul_pot_k);

write(*,*)"my_coul_pot_r",my_coul_pot_r(1:5)
write(*,*)"my_coul_pot_k",my_coul_pot_k(1:5)

      call gpmdcov_reallocate_realMat(q1,maxCoresAmongPartsAndRanks,partsInEachRank(myRank))
      call gpmdcov_reallocate_realMat(dqdmu,maxCoresAmongPartsAndRanks,partsInEachRank(myRank))
      
      trP1 = 0.0_dp; trdPdMu = 0.0_dp

      write(*,*)"K0ResPart15-3",K0ResPart(1:5,1)


#ifdef DO_MPI
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do iptt = 1,gpat%TotalParts
         ipt = iptt
#endif
        norbs = mysyprt(ipt)%estr%norbs 
        
        !We obtain the atoms in the Core+Halo
        !Can also be obtained from gpat%sgraph(ipt)%lsize
        natsCoreHalo = mysyprt(ipt)%nats
        natsCore = gpat%sgraph(ipt)%llsize

        call gpmdcov_reallocate_realVect(ptcoul_pot_k,natsCoreHalo)
        call gpmdcov_reallocate_realVect(ptcoul_pot_r,natsCoreHalo)
        call gpmdcov_reallocate_realVect(ptnet_charge,natsCoreHalo)

        ptcoul_pot_k = 0.0_dp
        ptcoul_pot_r = 0.0_dp
        ptnet_charge = 0.0_dp

        !Get Coulombic potential and charges for the part (C+H)
        do myj=1,natsCoreHalo
          myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
          ptcoul_pot_k(myj) = my_coul_pot_k(myjj)
          ptcoul_pot_r(myj) = my_coul_pot_r(myjj)
          ptnet_charge(myj) = chargePertVect(myjj)
        enddo

        !Extract the perturbation ovre the core part only
        do myj=1,natsCore
          myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
          v_core_i(myj,iptt,irank) = vi(myjj,irank) 
        enddo

        !COMMENT: MAYBE v_core_i should not be normalized!!!!!!

        !Computing the perturbative Hamiltonian (H1) from the coulombic perturbations
        !generated by dr=K0*res. At the C+H level
        write(*,*)"NUMORBS",norbs
        call gpmdcov_reallocate_denseBmlRealMat(ptham_bml,norbs)
        call gpmdcov_reallocate_denseBmlRealMat(ptaux_bml,norbs)
        !ptaux_bml corresponds to H0 which is 0 in this case.
        call prg_get_hscf_v2(ptaux_bml,mysyprt(ipt)%estr%over,ptham_bml,mysyprt(ipt)%spindex,&
             mysyprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,norbs,lt%threshold)

        !Compute transformations ZQ and (ZQ)^t transformation that takes from the canonical nonorthogonal
        !to the orthogonal eigenbasis. 
        call gpmdcov_reallocate_denseBmlRealMat(zq_bml,norbs)
        call bml_multiply(mysyprt(ipt)%estr%zmat,mysyprt(ipt)%estr%evects,zq_bml,1.0_dp,0.0_dp,lt%threshold)
        call gpmdcov_reallocate_denseBmlRealMat(zqt_bml,norbs)
        call bml_transpose(zq_bml,zqt_bml)

        !Take H1 to the ortho-eigen basis set.
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,0.0_dp)

        !Construct the "bare" response P1 and the derivative with respect to the
        !chemical potential (dPdMu). Everything in the ortho-eigen basis set
        call gpmdcov_reallocate_realVect(dPdMu,norbs)
        norbsCore = mysyprt(ipt)%estr%norbsCore
        call gpmdcov_reallocate_denseBmlRealMat(p1_bml,norbs)
        write(*,*)"DIMHQ1",bml_get_N(ptham_bml),size(mysyprt(ipt)%estr%evals),bml_get_N(p1_bml),norbsCore,bml_get_N(mysyprt(ipt)%estr%evects)&
&, bml_get_N(p1_bml),norbs
         
        call prg_canon_response_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
             &norbsCore,beta,mysyprt(ipt)%estr%evects,&
             &mysyprt(ipt)%estr%evals,ef,12,norbs)

        !At this point ptham is not needed anymore
        if(bml_allocated(ptham_bml)) call bml_deallocate(ptham_bml)

        !Transform P1 back to the nonortho-canonical basis set. 
        call bml_multiply(zq_bml,p1_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,p1_bml,2.0_dp,0.0_dp,0.0_dp)

        !Since dPdMu is represented by a vector, we will convert it into a BML matrix
        call gpmdcov_reallocate_denseBmlRealMat(dPdMu_bml,norbs)
        call bml_set_diagonal(dPdMu_bml,dPdMu)

        !Transform dPdMu back to the nonortho-canonical basis set
        call bml_multiply(zq_bml,dPdMu_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,dPdMu_bml,2.0_dp,0.0_dp,0.0_dp)

        !At this point ZQ and ZQt are not needed anymore
        call bml_deallocate(zq_bml)
        call bml_deallocate(zqt_bml)

        !Here we compute the charges response (q1) from P1 and we store it on 
        !a vector q1 that stores all the previous q1s from past iranks iterations
        !We also compute the partial trace contribution (trP1) from this mpi
        !execution and the current part (ipt).
        ptnet_charge = 0.0_dp
        mynumel = 0.0_dp
        call prg_get_charges(p1_bml, mysyprt(ipt)%estr%over,&
             &mysyprt(ipt)%estr%hindex, ptnet_charge, mynumel,&
             mysyprt(ipt)%spindex, norbs, lt%threshold)
        q1(:,iptt) = ptnet_charge(1:natsCore)
        trP1 = trP1 +  sum(ptnet_charge(1:natsCore))
        call bml_deallocate(p1_bml)

        !Here we compute the charges response (dqdmu) from dPdMu and we store
        !them on a matrix dqdmu that stores all the previous dqdmus from past
        !irank iterations.
        !We also compute the partial trace contribution (trdPdMu) from this node
        !and the current part (ipt).
        ptnet_charge = 0.0_dp
        call prg_get_charges(dpdmu_bml, mysyprt(ipt)%estr%over,&
             &mysyprt(ipt)%estr%hindex, ptnet_charge, mynumel,&
             mysyprt(ipt)%spindex, norbs, lt%threshold)
        dqdmu(:,iptt) = ptnet_charge(1:natsCore)
        trdPdMu = trdPdMu + sum(ptnet_charge(1:natsCore))
        call bml_deallocate(dPdMu_bml)
        deallocate(ptnet_charge)

     enddo !End of the loop over the parts

     !We will gather all the partial traces from all the MPI Ranks.  
#ifdef DO_MPI
     if (getNRanks() .gt. 1) then
        call prg_sumRealReduceN(trP1, 1)
        call prg_sumRealReduceN(trdPdMu, 1)
      endif
#endif


      write(*,*)"trP1,trdPdMu",trP1,trdPdMu
      mu1_Global = -trP1(1)/trdPdMu(1)
      write(*,*)"mu1_global",mu1_global
      q1 = q1 + mu1_Global*dqdmu
      write(*,*)"q1",q1(1:5,1)

      call gpmdcov_reallocate_realVect(f,maxCoresAmongPartsAndRanks)
      !deallocate(maxCoresAmongParts)
      f = 0.0_dp
      dr = 0.0_dp
#ifdef DO_MPI
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do iptt = 1,gpat%TotalParts
         ipt = iptt
#endif
        natsCore = gpat%sgraph(ipt)%llsize
        
        f = q1(:,iptt) - v_core_i(:,iptt,irank)
        ff(:,iptt,irank) = MATMUL(syprt(ipt)%estr%ker,f(1:natsCore))
        c_i(irank) = c_i(irank) + DOT_PRODUCT(ff(1:natsCore,iptt,irank),K0ResPart(1:natsCore,iptt))
        write(*,*)ff(1:5,iptt,irank)
        write(*,*)"K0ResPart15-4",K0ResPart(1:5,1)
        write(*,*)"natsCore",natsCore
        do myj=1,natsCore
          myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
          dr(myjj) = ff(myj,iptt,irank)
        enddo
      enddo
 
#ifdef DO_MPI
      if (getNRanks() .gt. 1) then
        call prg_sumRealReduceN(dr, sy%nats)
      endif
#endif
    dr_save(:,iRank) = dr

   enddo !Rank enddo

#ifdef DO_MPI
   if (getNRanks() .gt. 1) then
     call prg_sumRealReduceN(c_i, mRanks)
   endif
#endif

   allocate(auxVect(mRanks*mRanks))
   auxVect = 0.0_dp
   do myi = 1,mRanks
      do myj = 1,mRanks
        do iptt=1,partsInEachRank(myRank)
          auxVect((myi-1)*mRanks + myj) = & 
          &auxVect((myi-1)*mRanks + myj) + DOT_PRODUCT(ff(:,iptt,myi),ff(:,iptt,myj))
        enddo
      enddo
    enddo


#ifdef DO_MPI
   if (getNRanks() .gt. 1) then
     call prg_sumRealReduceN(auxVect, mRanks*mRanks)
   endif
#endif


   allocate(oij(mRanks,mRanks))
   allocate(mMat(mRanks,mRanks))
   do myi = 1,mRanks
      do myj = 1,mRanks
        oij(myi,myj) = auxVect((myi-1)*mRanks + myj)
      enddo
   enddo

   write(*,*)"OIJ",oij
   mMat = oij

   
   if(allocated(work))deallocate(work);allocate(work(mRanks+mRanks*mRanks))
   if(allocated(ipiv))deallocate(ipiv);allocate(ipiv(mRanks))
   call DGETRF(mRanks,mRanks,mMat,mRanks,ipiv,info)
   call DGETRI(mRanks,mMat,mRanks,ipiv,work,mRanks+mRanks*mRanks,info)
   deallocate(work)
   deallocate(ipiv)

   !Reconstruct the full v from the v_core_is
   allocate(vRank(sy%nats))
   do iptt=1,partsInEachRank(myRank)
     ipt= reshuffle(iptt,myRank)
     do myj=1,gpat%sgraph(ipt)%llsize
       myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
       vRank(myjj) = v_core_i(myj,iptt,irank)
     enddo
   enddo

#ifdef DO_MPI
   if (getNRanks() .gt. 1) then
     call prg_sumRealReduceN(vRank, sy%nats)
   endif
#endif
   
   call gpmdcov_reallocate_realVect(KK0Res,sy%nats)
   call gpmdcov_reallocate_realVect(IdK0Res,sy%nats)
   KK0Res = 0.0_dp
   IdK0Res = 0.0_dp
   do myi = 1,mRanks
     do myj = 1,mRanks
       KK0Res = KK0Res + vi(:,myi)*mMat(myi,myj)*c_i(myj)
       IdK0Res = IdK0Res + dr_save(:,myi)*mMat(myi,myj)*c_i(myj)
     enddo
   enddo
!   dn2dt2 = -(w**2)*KK0Res
   write(*,*)"RELERR",norm2(IdK0Res - K0Res)/norm2(IdK0Res)
   write(*,*)mMat
   write(*,*)c_i(1:mRanks)
    write(*,*)"KKDIFF"
    write(*,*)KK0Res(1:5) - K0Res(1:5)
    write(*,*)KK0Res(1:5)
    write(*,*)K0Res(1:5)
    write(*,*)"d1/ff1",dr_save(1:5,1)
    !write(*,*)"d1/ff1",dr_save(1:5,2)

   do myi = 1,mRanks
      do myj = 1,mRanks
          write(*,*) DOT_PRODUCT(dr_save(:,myi),dr_save(:,myj))
      enddo
    enddo

   deallocate(my_coul_pot_k)
   deallocate(my_coul_pot_r)
   deallocate(my_coul_forces_k)
   deallocate(my_coul_forces_r)
   deallocate(K0Res)
   deallocate(K0ResPart)

!stop

  end subroutine gpmdcov_rankN_update_byParts


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CACULATES THE PRECONDITIONED MULTI-RANK APPROXIMATION OF KERNEL ACTING ON THE
!RESIDUAL !!
!!                  THEORY GIVEN IN Niklasson, JCP 152, 104103 (2020) [*]
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine KK0Res(KRes,KK0,Res,FelTol,L,LL,H0,mu0,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,MaxIt,eps,nr_mom,&
       &HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START,H_INDEX_END,&
       &H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

!! Res = q[n] - n
!! KK0 is preconditioner
!! KRes rank-L approximation of (K0*J)^(-1)*K0*(q[n]-n) with (K0*J)^(-1) as in Eq. (41) in Ref. [*]
!! QQ, ee Eigenvectors and eigen values of H0
!! Fe_vec Fermi occupation factors

integer, parameter          :: PREC = 8
integer,    intent(in)      :: Nr_atoms, HDIM, Nocc,Max_Nr_Neigh
real(PREC), parameter       :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
real(PREC), parameter       :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
real(PREC), intent(in)      :: Coulomb_acc, TIMERATIO, FelTol
real(PREC), intent(in)      :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
real(PREC), intent(in)      :: KK0(Nr_atoms,Nr_atoms), Res(Nr_atoms)
integer,    intent(in)      :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)      :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
real(PREC), intent(in)      :: H(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM),H0(HDIM,HDIM) 
real(PREC)                  :: H1(HDIM,HDIM), HX(HDIM,HDIM), K0Res(Nr_atoms)
real(PREC)                  :: D0(HDIM,HDIM)
real(PREC)                  :: X(HDIM,HDIM)
real(PREC), intent(in)      :: T, eps
real(PREC), intent(inout)   :: mu0
integer,    intent(in)      :: nr_mom, MaxIt
character(10), intent(in)   :: Element_Type(Nr_atoms)
integer,    intent(in)      :: nrnnlist(Nr_atoms),nnType(Nr_atoms,Max_Nr_Neigh), LL
real(PREC), intent(in)      :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)      :: nnRy(Nr_atoms,Max_Nr_Neigh),nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC)                  :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Pot(Nr_atoms)
real(PREC)                  :: Coulomb_Pot_Real_I,Coulomb_Pot_k(Nr_atoms),dq_dv(Nr_atoms)
real(PREC)                  :: Coulomb_Pot_Real_dq_v(Nr_atoms),Coulomb_Pot_dq_v(Nr_atoms)
real(PREC)                  :: Coulomb_Force_Real_I(3),Coulomb_Force_k(3,Nr_atoms)
real(PREC)                  :: D_dq_v(HDIM,HDIM), H_dq_v(HDIM,HDIM),dq_v(Nr_atoms), ZQ(HDIM,HDIM)
real(PREC), intent(out)     :: KRes(Nr_atoms)
real(PREC), intent(in)      :: QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM)
integer                     :: I,J,K
integer, intent(out)        :: L  !! Number of rank updates used to reach below threshold
real(PREC)                  :: Fel,mu1
real(PREC)                  :: vi(Nr_atoms,Nr_atoms), v(Nr_atoms),q_tmp(Nr_atoms)
real(PREC)                  :: fi(Nr_atoms,Nr_atoms)
real(PREC)                  :: dr(Nr_Atoms), proj_tmp, IdentRes(Nr_atoms)
real(PREC), allocatable     :: O(:,:),M(:,:)

     K0Res = MATMUL(KK0,Res)  !! EXTRA STEP FOR PRECONDITIONING
     dr = K0Res               !! Preconditioned residual
     I = 0                    !! Count number of rank updates
     Fel = 1.D0
     !do while (Fel > FelTol)   !! Fel = "Error" in Swedish, Could potentially
     !also use a highest number of allowed rank updates
     do while ((Fel > FelTol).AND.(I < LL))   !! Fel = "Error" in Swedish, Could potentially also use a highest number of allowed rank updates
        I = I + 1

        vi(:,I) = dr/norm2(dr)
        do J = 1,I-1
           vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)  !!Orthogonalized v_i as in Eq. (42) Ref. [*]
        enddo
        vi(:,I) = vi(:,I)/norm2(vi(:,I))
        v(:) = vi(:,I)  ! v_i

        !! Calculated dq_dv, which is the response in q(n) from change in
        !directional input charge n = v
       ! call V_Kernel_Fermi_No_Diag(D0,dq_dv,v,H0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,&
       !& MaxIt,eps,nr_mom,HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,&
       !& H_INDEX_START,H_INDEX_END,H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

        dr = dq_dv - v       !! dr = df/dlambda, last row in Eq. (42) Ref[*]
        dr = MATMUL(KK0,dr)  !! dr = K0*(df/dlambda), last row in Eq. (42) Ref[*]
        fi(:,I) = dr  ! fv_i

        L = I
        allocate(O(L,L), M(L,L))
        do K = 1,L
        do J = 1,L
           O(K,J) = dot_product(fi(:,K),fi(:,J))  ! O_KJ = < fv_i(K) | fv_i(J) > see below Eq. (31)
        enddo
        enddo
 !       call Invert(O,M,L)                        ! M = O^(-1)
        IdentRes = 0.D0*K0Res
        KRes = 0.D0
        do K = 1,L
        do J = 1,L
           proj_tmp = M(K,J)*dot_product(fi(:,J),K0Res)
           IdentRes = IdentRes + proj_tmp*fi(:,K)
           KRes = KRes + proj_tmp*vi(:,K)            !! KRes becomes the rank-L approximate of (K0*J)^(-1)*K0*(q[n]-n) 
        enddo
        enddo
        Fel = norm2(IdentRes-K0Res)/norm2(IdentRes)  !! RELATIVE RESIDUAL ERROR ESTIMATE Eq. (48) Ref. [*]
        write(*,*) '## I, L, Fel = ',I,L,Fel         !! Fel goes down with L
        deallocate(O, M)
     enddo

end subroutine KK0Res

  subroutine gpmdcov_applyKernel(my_nguess,my_chargesOld,kernelTimesRes,control)
    use gpmdcov_vars
    implicit none
    real(dp), allocatable, intent(inout) :: kernelTimesRes(:)
    real(dp), allocatable, intent(in) :: my_nguess(:)
    real(dp), allocatable, intent(in) ::  my_chargesOld(:)
    real(dp), allocatable :: nguessPart(:), chargesOldPart(:)
    real(dp), allocatable :: kernelTimesResPart(:)
    integer :: nats, myj, myjj
    logical, optional  :: control

    if(present(control))then
      write(*,*)"minikernel"
      write(*,*) syprt(1)%estr%ker
      write(*,*)"my_chargesOld", my_chargesOld
      write(*,*) "my_nguess",my_nguess
    endif
    !> Get Coulombic potential and charges for the part.
    nats = size(my_nguess,dim=1)
    if(.not. allocated(kernelTimesRes)) allocate(kernelTimesRes(nats))
    kernelTimesRes = 0.0_dp
    write(*,*)"nats",nats
    write(*,*)"size my_chargesOld",size(my_chargesOld,dim=1)
    write(*,*)"size my_nguess",size(my_nguess,dim=1)

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
        !         write(*,*)"myj myjj llsize nats",myj,myjj,gpat%sgraph(ipt)%llsize,size(my_chargesOld,dim=1),size(my_nguess,dim=1)
        chargesOldPart(myj) = my_chargesOld(myjj)
        nguessPart(myj) = my_nguess(myjj)
      enddo
      !        if(present(control))stop
      write(*,*)"MDBUG1",size(nguessPart),size(chargesOldPart),size(syprt(ipt)%estr%ker)
      kernelTimesResPart = MATMUL(syprt(ipt)%estr%ker,(nguessPart-chargesOldPart))

     !  write(*,*)syprt(ipt)%estr%ker(1,1:5)
     ! write(*,*)nguessPart(1:5),chargesOldPart(1:5)
     ! stop


      !Expand  
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
    call prg_wait
#endif


  end subroutine gpmdcov_applyKernel

  subroutine MMult(alpha,A,B,beta,C,TTA,TTB,HDIM)

    implicit none
    integer,      parameter     :: PREC = 8
    integer,      intent(in)    :: HDIM
    real(PREC),   intent(inout)    :: A(HDIM,HDIM), B(HDIM,HDIM), alpha, beta
    real(PREC),   intent(inout) :: C(HDIM,HDIM)
    character(1), intent(in)    :: TTA, TTB

    !write(*,*) ' FORE A(1,1) = ', A(1,1), A(2,2)
    !DO I = 1,HDIM
    !  A(I,I) = floor(A(I,I)*100000000.D0)/100000000.D0
    !ENDDO
    !write(*,*) ' EFTER A(1,1) = ', A(1,1), A(2,2)

    if (PREC.eq.4) then
      call SGEMM(TTA, TTB, HDIM, HDIM, HDIM, alpha, &
           A, HDIM, B, HDIM, beta, C, HDIM)
    else
      call DGEMM(TTA, TTB, HDIM, HDIM, HDIM, alpha, &
           A, HDIM, B, HDIM, beta, C, HDIM)

    endif

  end subroutine MMult


end module gpmdcov_kernel_mod
