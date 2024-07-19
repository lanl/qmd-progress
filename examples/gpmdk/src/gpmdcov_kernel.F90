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
  use gpmdcov_writeout_mod

#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif
    
  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_getKernel, gpmdcov_getKernel_byBlocks, gpmdcov_parseKernel
  public :: gpmdcov_applyKernel, gpmdcov_getKernel_byParts, gpmdcov_rankN_update_byParts

  !> General latte input variables type.
  !!
  type, public :: kernel_type

    !> Kernel type (if full or other)
    character(100) :: kernelType
    
    !> Kernel initial mixing method 
    character(100) :: initialMixingWith

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
    
    !> If Update after build
    logical :: updateAfterBuild

    !> XLBO level. This adds an extra scf which increases stability
    logical :: xlbolevel1

  end type kernel_type

  type(kernel_type), public :: kernel

contains


  !> The parser for the kernel type.
  !!
  subroutine gpmdcov_parseKernel(kernel,filename)

    implicit none
    type(kernel_type) :: kernel
    integer, parameter :: nkey_char = 2, nkey_int = 3, nkey_re = 1, nkey_log = 4
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         'KernelType=','InitialMixingWith=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'Full','Lin']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Verbose=', 'UpdateEach=', 'RankNUpdate=']
    integer :: valvector_int(nkey_int) = (/ &
         0   ,     0,    0  /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Threshold=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.00001 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         'KernelMixing=','BuildAlways=','UpdateAfterBuild=',"XLBOLevel1="]
    logical :: valvector_log(nkey_log) = (/&
         .false., .false.,.false.,.false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'KERNEL{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    kernel%kernelType = valvector_char(1)
    kernel%initialMixingWith = valvector_char(2)

    !Integers
    kernel%verbose = valvector_int(1)
    kernel%updateEach = valvector_int(2)
    kernel%rankNUpdate = valvector_int(3)

    !Reals
    kernel%threshold = valvector_re(1)

    !Logical
    kernel%kernelmixing = valvector_log(1)
    kernel%buildAlways = valvector_log(2)
    kernel%updateAfterBuild = valvector_log(3)
    kernel%xlbolevel1 = valvector_log(4)

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
    real(dp), allocatable :: my_coul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_r(:,:), my_coul_pot_k(:), my_coul_pot_r(:)
    real(dp), allocatable :: Jacob(:), work(:), ptcoul_pot_k(:),ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:)
    real(dp), allocatable :: mynumel(:)
    real(dp) :: mlsi, KSUM
    integer :: nats, l, m, lm, info, norbs
    integer, allocatable :: ipiv(:)
    type(bml_matrix_t) :: zq_bml, zqt_bml
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
      call get_ewald_list_real_dcalc_vect(sy%spindex,sy%splist,sy%coordinate&
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

        call bml_multiply(syprt(ipt)%estr%zmat,syprt(ipt)%estr%evects,zq_bml, 1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)

     !   call prg_canon_response(ptrho_bml,ptham_bml,norbs,beta,syprt(ipt)%estr%evects,&
     !        &syprt(ipt)%estr%evals,ef,12,norbs)

            call prg_canon_response(ptrho_bml,ptham_bml,1.0_dp,beta,&
             &syprt(ipt)%estr%evals,ef,12,norbs)

        call gpmdcov_msI("gpmdcov_get_kernel","Time for Canonincal Response construction "&
                &//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

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
      !write(*,*) ' KSUM ', m, ' = ',KSum
    enddo

  end subroutine gpmdcov_getKernel


  subroutine gpmdcov_getKernel_byBlocks(nats)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: ptcoul_forces_k(:,:)
    real(dp), allocatable :: my_coul_forces_k(:,:),my_coul_forces_r(:,:)
    real(dp), allocatable ::  my_coul_pot_k(:), my_coul_pot_r(:)
    real(dp), allocatable :: Jacob(:), work(:), ptcoul_pot_k(:),ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:)
    real(dp), allocatable :: mynumel(:)
    real(dp) :: mlsi
    integer :: nats, l, m, lm, info, atom, norbs
    integer, allocatable :: ipiv(:)
    type(bml_matrix_t) :: zq_bml, zqt_bml
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

        call get_ewald_list_real_dcalc_vect(sy%spindex,sy%splist,sy%coordinate&
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

        mlsi = mls()
        call bml_multiply(syprt(ipt)%estr%zmat,syprt(ipt)%estr%evects,zq_bml, 1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)
        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for tranf to eigen basis " &
                &//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        mlsi = mls()
        !call prg_canon_response(ptrho_bml,ptham_bml,norbs,23.208882712264995_dp,syprt(ipt)%estr%evects,&
        !     &syprt(ipt)%estr%evals,ef,12,norbs)
        call prg_canon_response(ptrho_bml,ptham_bml,1.0_dp,beta,&
            &syprt(ipt)%estr%evals,ef,12,norbs)
        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for Canonincal Response construction " &
                &//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        mlsi = mls()
        call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)
        call gpmdcov_msI("gpmdcov_get_kernel_byBlocks","Time for back to canonical basis " &
                &//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

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

        call gpmdcov_msIII("gpmdcov_get_kernel_byBlocks","Time for getting charges"//to_string(mls() &
                &- mlsi)//" ms",lt%verbose,myRank)
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

   ! do m = 1,Nats
   !   KSum = 0.D0
   !   do l = 1,Nats
   !     KSum = KSum + Ker(l,m)
   !   enddo
   !   write(*,*) ' KSUM ', m, ' = ',KSum
   ! enddo

  end subroutine gpmdcov_getKernel_byBlocks

  subroutine gpmdcov_getKernel_byParts(mysyprt,mysyprtk)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: my_coul_forces_k(:,:),my_coul_forces_r(:,:)
    real(dp), allocatable :: my_coul_pot_k(:)
    real(dp), allocatable :: my_coul_pot_r(:)
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: ptcoul_pot_k(:),ptcoul_pot_r(:)
    real(dp), allocatable :: ptnet_charge(:)
    real(dp) :: mynumel(10)
    real(dp) :: mlsi, trdPdMuAO, trp1, mu1, mls_v
    integer :: info, atom, coreSize, norbsCore, norbs
    integer, allocatable :: ipiv(:)
    type(system_type), allocatable, intent(inout) :: mysyprt(:)
    type(system_type), allocatable, intent(inout) :: mysyprtk(:)
    type(bml_matrix_t) :: zq_bml, zqt_bml
    type(bml_matrix_t) :: ptham_bml, ptrho_bml, ptaux_bml
    type(bml_matrix_t) :: dPdMuAO_bml,p1_bml,dPdMuAOS_bml,p1S_bml
    real(dp), allocatable :: dPdMuAO_dia(:),p1_dia(:),dPdMu(:)

    mls_v = mls() 

    if(.not.allocated(my_coul_forces_k))allocate(my_coul_forces_k(3,sy%nats))
    if(.not.allocated(my_coul_forces_r))allocate(my_coul_forces_r(3,sy%nats))
    if(.not.allocated(my_coul_pot_k))allocate(my_coul_pot_k(sy%nats))
    if(.not.allocated(my_coul_pot_r))allocate(my_coul_pot_r(sy%nats))
    if(allocated(mysyprtk))deallocate(mysyprtk)
    allocate(mysyprtk(gpat%TotalParts))
    getKernel_byParts_cont = getKernel_byParts_cont + 1 

    call gpmdcov_msI("gpmdcov_get_kernel_byParts", to_string(myRank)//"has "//&
             &to_string(partsInEachRank(myRank))//" Parts",lt%verbose,1)
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

      if(allocated(mysyprtk(ipt)%estr%ker)) deallocate(mysyprtk(ipt)%estr%ker)
      allocate(mysyprtk(ipt)%estr%ker(coreSize,coreSize))

      mysyprtk(ipt)%estr%ker = 0.0_dp

      do i=1, coreSize
        mlsi = mls()
        chargePertVect=0.0_dp
        atom = gpat%sgraph(ipt)%core_halo_index(i)+1
        chargePertVect(atom)=1.0_dp

        call gpmdcov_msII("gpmdcov_get_kernel_byParts","Constructing response&
             & for atom ="//to_string(atom),lt%verbose,myRank)

        call get_ewald_list_real_dcalc_vect(sy%spindex,sy%splist,sy%coordinate&
             ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
             sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
             nl%nnIz,nl%nrnnlist,nl%nnType,my_coul_forces_r,my_coul_pot_r);
        
        call gpmdcov_msII("gpmdcov_get_kernel_byParts","Time for coulomb real &
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)
        
        mlsi = mls()

        call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
             ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
             sy%recip_vector,sy%volr,lt%coul_acc,my_coul_forces_k,my_coul_pot_k);

        call gpmdcov_msII("gpmdcov_get_kernel_byParts","Time for coulomb recip &
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        allocate(ptcoul_pot_k(mysyprt(ipt)%nats))
        allocate(ptcoul_pot_r(mysyprt(ipt)%nats))
        allocate(ptnet_charge(mysyprt(ipt)%nats))

        ptcoul_pot_k = 0.0_dp
        ptcoul_pot_r = 0.0_dp
        ptnet_charge = 0.0_dp

        !> Get Coulombic potential and charges for the core part.
        do j=1,mysyprt(ipt)%nats
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          ptcoul_pot_k(j) = my_coul_pot_k(jj)
          ptcoul_pot_r(j) = my_coul_pot_r(jj)
          ptnet_charge(j) = chargePertVect(jj)
        enddo


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
        call gpmdcov_msIII("gpmdcov_getKernel_byParts","Entering prg_get_hscf to&
             &construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(ptaux_bml,mysyprt(ipt)%estr%over,ptham_bml,mysyprt(ipt)%spindex,&
             mysyprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,lt%mdim,lt%threshold)


        deallocate(ptcoul_pot_r)
        deallocate(ptcoul_pot_k)

        call gpmdcov_msIII("gpmdcov_get_kernel_byParts","Time for H construction&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)
        
        mlsi = mls()
        call bml_multiply(mysyprt(ipt)%estr%zmat,mysyprt(ipt)%estr%evects,zq_bml,&
             &1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,lt%threshold)
        call gpmdcov_msIII("gpmdcov_get_kernel_byParts","Time for trasnf to eigenbasis&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        mlsi = mls()
        allocate(dPdMu(norbs))
        !call prg_canon_response2_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
        !     &norbsCore,beta,mysyprt(ipt)%estr%evects,&
        !     &mysyprt(ipt)%estr%evals,ef,12,norbs,lt%threshold)
        call prg_canon_response_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
             &norbsCore,beta,mysyprt(ipt)%estr%evects,&
             &mysyprt(ipt)%estr%evals,ef,12,norbs)
        call bml_get_diagonal(p1_bml,p1_dia)
        trP1 = sum(p1_dia(1:norbsCore))
        call gpmdcov_msII("gpmdcov_get_kernel_byParts","Time for Canonincal&
             &Response construction "//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,dPdMuAO_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,dPdMuAOS_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norbs,norbs,p1S_bml)
        call bml_set_diagonal(dPdMuAO_bml,dPdMu)

        deallocate(dPdMu)

        mlsi = mls()
        call bml_multiply(zq_bml,p1_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,p1_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(p1_bml,syprt(ipt)%estr%over,p1S_bml,1.0_dp,0.0_dp,0.0_dp)

        call bml_multiply(zq_bml,dPdMuAO_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,dPdMuAO_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(dPdMuAO_bml,syprt(ipt)%estr%over,dPdMuAOS_bml,1.0_dp,0.0_dp,0.0_dp)
        call gpmdcov_msIII("gpmdcov_get_kernel_byParts","Time for trasnf to canonical&
             &"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

        call bml_get_diagonal(dPdMuAOS_bml,dPdMuAO_dia)
        call bml_get_diagonal(p1S_bml,p1_dia)
        trP1 = sum(p1_dia(1:norbsCore))
        !write(*,*)"trP1",trP1
        trdPdMuAO = sum(dPdMuAO_dia(1:norbsCore))
        deallocate(dPdMuAO_dia)
        deallocate(p1_dia)
        if(abs(trdPdMuAO) < 1.0E-12)then
              mu1 = 0.0
        else
              mu1 =  -trP1/trdPdMuAO
        endif
        !write(*,*)"trP1,trdPdMuAO,mu1",trP1,trdPdMuAO,mu1
        call bml_copy(p1_bml,ptrho_bml)
        call bml_add(ptrho_bml,dPdMuAO_bml,2.0_dp,2.0_dp*mu1,lt%threshold)

        call bml_deallocate(ptham_bml)
        call bml_deallocate(zq_bml)
        call bml_deallocate(zqt_bml)
        call bml_deallocate(ptaux_bml)
        call bml_deallocate(p1_bml)
        call bml_deallocate(dPdMuAO_bml)
        call bml_deallocate(dPdMuAOS_bml)
        call bml_deallocate(p1S_bml)

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
        call gpmdcov_msII("gpmdcov_get_kernel_byParts","Time for getting &
             &charges"//to_string(mls() - mlsi)//" ms",lt%verbose,myRank)

    !  do m = 1,coreSize
    !    KSum = 0.0_dp
    !    do l = 1,coreSize
    !      KSum = KSum + ptnet_charge(l)
    !    enddo
    !    write(*,*) ' JSUM ', m, ' = ',KSum
    !  enddo           
        
        !Constructing J to prepare for K
        do j=1,coreSize
          mysyprtk(ipt)%estr%ker(j,i) = ptnet_charge(j)
          if(i == j)mysyprtk(ipt)%estr%ker(j,i) = mysyprtk(ipt)%estr%ker(j,i) - 1.0_dp
          !Regularized Kernel
          !if(i == j)mysyprtk(ipt)%estr%ker(j,i) = mysyprtk(ipt)%estr%ker(j,i) - 1.2_dp
        enddo
        deallocate(ptnet_charge)
        call bml_deallocate(ptrho_bml)

      enddo

      if(allocated(work))deallocate(work);allocate(work(coreSize+coreSize*coreSize))
      if(allocated(ipiv))deallocate(ipiv);allocate(ipiv(coreSize))
      call DGETRF(coreSize, coreSize, mysyprtk(ipt)%estr%ker, coreSize, ipiv, info)
      call DGETRI(coreSize, mysyprtk(ipt)%estr%ker, coreSize, ipiv,&
           & work, coreSize+coreSize*coreSize, info)
      deallocate(work)
      deallocate(ipiv)

    !  do m = 1,coreSize
    !    KSum = 0.0_dp
    !    do l = 1,coreSize
    !      KSum = KSum + mysyprtk(ipt)%estr%ker(l,m)
    !    enddo
    !    write(*,*) ' KSUM ', m, ' = ',KSum
    !  enddo

! Hack to use the identity (to see the effect of low rank updates)
      !mysyprtk(ipt)%estr%ker = 0.0_dp
      !write(*,*)"WARNING!!!!!!! using Identity as precond"
      !do m = 1,coreSize
      !  mysyprtk(ipt)%estr%ker(m,m) = -0.2_dp
      !enddo
    enddo

    deallocate(my_coul_forces_k)
    deallocate(my_coul_forces_r)
    deallocate(my_coul_pot_k)
    deallocate(my_coul_pot_r)

    call gpmdcov_msI("gpmdcov_get_kernel_byParts - Rank "//to_string(myRank) ,"Time for get_kernel_byParts &
             &"//to_string(mls() - mls_v)//" ms",lt%verbose,myRank)


    
  end subroutine gpmdcov_getKernel_byParts


  subroutine gpmdcov_rankN_update_byParts(myqn,myn,mysyprt,mysyprtk,maxRanks,KK0Res)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    use gpmdcov_allocation_mod
    implicit none
    integer :: myi,myj,myjj,irank,mRanks
    integer :: norbsCore,norbs
    integer, save :: maxorbs = -1
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
    real(dp) :: mu1_global, error
    type(bml_matrix_t),save :: ptham_bml, ptaux_bml
    type(bml_matrix_t),save :: zq_bml,zqt_bml
    type(bml_matrix_t),save :: p1_bml,dPdMu_bml
    type(system_type), allocatable, intent(inout) :: mysyprt(:)
    type(system_type), allocatable, intent(inout) :: mysyprtk(:)
    real(dp) :: mynumel(10), mls_v
    real(dp), allocatable :: work(:),auxVect(:)
    integer :: info
    integer, allocatable :: ipiv(:)
   
    call gpmdcov_msI("gpmdcov_rankN_update_byParts","Updating the Kernel",lt%verbose,myRank)
    rankN_update_byParts_cont = rankN_update_byParts_cont + 1 
 
    mls_v = mls()

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

      call gpmdcov_reallocate_realVect(myqnPart,natsCore)
      call gpmdcov_reallocate_realVect(mynPart,natsCore)

      !Get old charges and nguess for the part
      do myj=1,natsCore
        myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
        mynPart(myj) = myn(myjj)
        myqnPart(myj) = myqn(myjj)
      enddo

      K0ResPart(1:natsCore,iptt) = MATMUL(mysyprtk(ipt)%estr%ker,(myqnPart-mynPart))

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

    mRanks = maxRanks !Number of total rank updates 
    allocate(vi(sy%nats,mRanks))
    allocate(dr_save(sy%nats,mRanks))
    vi = 0.0_dp
    allocate(v_core_i(maxCoresAmongPartsAndRanks,partsInEachRank(myRank),mRanks))
    v_core_i = 0.0_dp
    allocate(c_i(mRanks))
    c_i = 0.0_dp
    allocate(ff(maxCoresAmongPartsAndRanks,partsInEachRank(myRank),mRanks))
    ff = 0.0_dp

    if (maxorbs == -1) then
       maxorbs = min(maxCoresAmongPartsAndRanks*20,sy%nats)*4
    endif
      
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
#ifdef USE_NVTX
      call nvtxStartRange("Ewald",1)
#endif
      call get_ewald_list_real_dcalc_vect(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,my_coul_forces_r,my_coul_pot_r);

      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,my_coul_forces_k,my_coul_pot_k);
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      call gpmdcov_reallocate_realMat(q1,maxCoresAmongPartsAndRanks,partsInEachRank(myRank))
      call gpmdcov_reallocate_realMat(dqdmu,maxCoresAmongPartsAndRanks,partsInEachRank(myRank))
      
      trP1 = 0.0_dp; trdPdMu = 0.0_dp

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

        !Extract the perturbation over the core part only
        do myj=1,natsCore
          myjj = gpat%sgraph(ipt)%core_halo_index(myj)+1
          v_core_i(myj,iptt,irank) = vi(myjj,irank) 
        enddo
#ifdef USE_NVTX
        call nvtxStartRange("BML_math",2)
#endif
        !COMMENT: MAYBE v_core_i should not be normalized!!!!!!

        !Computing the perturbative Hamiltonian (H1) from the coulombic perturbations
        !generated by dr=K0*res. At the C+H level
#ifdef USE_NVTX
        call nvtxStartRange("Reallocation event",7)
#endif
        if (.not.bml_allocated(ptham_bml))then
           call gpmdcov_reallocate_denseBmlRealMat(ptham_bml,maxorbs)
           call gpmdcov_reallocate_denseBmlRealMat(zq_bml,maxorbs)
           call gpmdcov_reallocate_denseBmlRealMat(zqt_bml,maxorbs)
        endif        
        call gpmdcov_bml_set_N(ptham_bml,norbs)
        call gpmdcov_bml_set_N(zq_bml,norbs)
        call gpmdcov_bml_set_N(zqt_bml,norbs)
!        call gpmdcov_reallocate_denseBmlRealMat(ptham_bml,norbs)
#ifdef USE_NVTX
        call nvtxEndRange
#endif
        call gpmdcov_reallocate_denseBmlRealMat(ptaux_bml,norbs)
        !ptaux_bml corresponds to H0 which is 0 in this case.
        call prg_get_hscf_v2(ptaux_bml,mysyprt(ipt)%estr%over,ptham_bml,mysyprt(ipt)%spindex,&
             mysyprt(ipt)%estr%hindex,tb%hubbardu,ptnet_charge,&
             ptcoul_pot_r,ptcoul_pot_k,norbs,lt%threshold)

        !Compute transformations ZQ and (ZQ)^t transformation that takes from the canonical nonorthogonal
        !to the orthogonal eigenbasis. 
        call bml_multiply(mysyprt(ipt)%estr%zmat,mysyprt(ipt)%estr%evects,zq_bml,1.0_dp,0.0_dp,lt%threshold)
        call bml_transpose(zq_bml,zqt_bml)

        !Take H1 to the ortho-eigen basis set.
        call bml_multiply(zqt_bml,ptham_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,0.0_dp)

        !Construct the "bare" response P1 and the derivative with respect to the
        !chemical potential (dPdMu). Everything in the ortho-eigen basis set
        call gpmdcov_reallocate_realVect(dPdMu,norbs)
        norbsCore = mysyprt(ipt)%estr%norbsCore
        call gpmdcov_reallocate_denseBmlRealMat(p1_bml,norbs)
#ifdef USE_NVTX
        call nvtxEndRange
#endif
        !call prg_canon_response2_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
        !     &norbsCore,beta,mysyprt(ipt)%estr%evects,&
        !     &mysyprt(ipt)%estr%evals,ef,12,norbs,lt%threshold)
#ifdef USE_NVTX
        call nvtxStartRange("prg_canon_response_p1_dpdmu",3)
#endif
        call prg_canon_response_p1_dpdmu(p1_bml,dPdMu,ptham_bml,&
            &norbsCore,beta,mysyprt(ipt)%estr%evects,&
            &mysyprt(ipt)%estr%evals,ef,12,norbs)
#ifdef USE_NVTX
        call nvtxEndRange
#endif
#ifdef USE_NVTX
        call nvtxStartRange("BML_math_2",4)
#endif
        !At this point ptham is not needed anymore
        !if(bml_allocated(ptham_bml)) call bml_deallocate(ptham_bml)

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
        !call bml_deallocate(zq_bml)
        !call bml_deallocate(zqt_bml)
#ifdef USE_NVTX
        call nvtxEndRange
#endif
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
        call bml_deallocate(ptaux_bml)
        deallocate(ptnet_charge)

     enddo !End of the loop over the parts

     !We will gather all the partial traces from all the MPI Ranks.  
#ifdef DO_MPI
     if (getNRanks() .gt. 1) then
        call prg_sumRealReduceN(trP1, 1)
        call prg_sumRealReduceN(trdPdMu, 1)
      endif
#endif

      mu1_Global = -trP1(1)/trdPdMu(1)
      q1 = q1 + mu1_Global*dqdmu

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
        ff(:,iptt,irank) = MATMUL(syprtk(ipt)%estr%ker,f(1:natsCore))
        c_i(irank) = c_i(irank) + DOT_PRODUCT(ff(1:natsCore,iptt,irank),K0ResPart(1:natsCore,iptt))
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
       vRank(myjj) = v_core_i(myj,iptt,mRanks)
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

   error = norm2(IdK0Res - K0Res)/norm2(K0Res)
   if(myRank == 1)write(*,*)"Error Rank-Update",error,mRanks

   !do myi = 1,mRanks
   !   do myj = 1,mRanks
   !       write(*,*) DOT_PRODUCT(dr_save(:,myi),dr_save(:,myj))
   !   enddo
   ! enddo
#ifdef USE_NVTX
        call nvtxStartRange("Deallocations",5)
#endif


  !Deallocations 
    if(allocated(maxCoresAmongParts))deallocate(maxCoresAmongParts)
    if(allocated(myqnPart))deallocate(myqnPart)
    if(allocated(mynPart))deallocate(mynPart)
    if(allocated(dr))deallocate(dr)
    if(allocated(v))deallocate(v)
    if(allocated(my_coul_pot))deallocate(my_coul_pot)
    if(allocated(ptcoul_forces_k))deallocate(ptcoul_forces_k)
    if(allocated(ptcoul_forces_r))deallocate(ptcoul_forces_r)
    if(allocated(my_coul_pot_k))deallocate(my_coul_pot_k)
    if(allocated(chargePertVect))deallocate(chargePertVect)
    if(allocated(ptcoul_pot_k))deallocate(ptcoul_pot_k)
    if(allocated(ptcoul_pot_r))deallocate(ptcoul_pot_r)
    if(allocated(ptnet_charge))deallocate(ptnet_charge)
    if(allocated(v_core_i))deallocate(v_core_i)
    if(allocated(dPdMu))deallocate(dPdMu)
    if(allocated(IdK0Res))deallocate(IdK0Res)
    if(allocated(q1))deallocate(q1)
    if(allocated(dqdmu))deallocate(dqdmu)
    if(allocated(f))deallocate(f)
    if(allocated(c_i))deallocate(c_i)
    if(allocated(ff))deallocate(ff)
    if(allocated(oij))deallocate(oij)
    if(allocated(mMat))deallocate(mMat)
    if(allocated(vRank))deallocate(vRank)
    if(allocated(dr_save))deallocate(dr_save)
    if(allocated(vi))deallocate(vi)
    if(allocated(maxCoresAmongParts))deallocate(maxCoresAmongParts)
    if(allocated(my_coul_pot_k))deallocate(my_coul_pot_k)
    if(allocated(my_coul_pot_r))deallocate(my_coul_pot_r)
    if(allocated(my_coul_forces_k))deallocate(my_coul_forces_k)
    if(allocated(my_coul_forces_r))deallocate(my_coul_forces_r)
    if(allocated(K0Res))deallocate(K0Res)
    if(allocated(K0ResPart))deallocate(K0ResPart)
#ifdef USE_NVTX
    call nvtxEndRange
#endif
   call gpmdcov_msI("gpmdcov_get_kernel_byParts","Time for get_kernel_byParts &
             &"//to_string(mls() - mls_v)//" ms",lt%verbose,myRank)

  end subroutine gpmdcov_rankN_update_byParts


  subroutine gpmdcov_applyKernel(my_nguess,my_chargesOld,mysyprtk,kernelTimesRes)
    use gpmdcov_vars
    implicit none
    real(dp), allocatable, intent(inout) :: kernelTimesRes(:)
    real(dp), allocatable, intent(in) :: my_nguess(:)
    real(dp), allocatable, intent(in) ::  my_chargesOld(:)
    real(dp), allocatable :: nguessPart(:), chargesOldPart(:)
    real(dp), allocatable :: kernelTimesResPart(:)
    integer :: nats, myj, myjj
    type(system_type), allocatable, intent(in) :: mysyprtk(:)

    call gpmdcov_msI("gpmdcov_applyKernel","Applying the kernel...",lt%verbose,myRank)   
 
    !> Get Coulombic potential and charges for the part.
    nats = size(my_nguess,dim=1)
    if(.not. allocated(kernelTimesRes)) allocate(kernelTimesRes(nats))
    kernelTimesRes = 0.0_dp

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
      kernelTimesResPart = MATMUL(mysyprtk(ipt)%estr%ker,(nguessPart-chargesOldPart))

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

call gpmdcov_msI("gpmdcov_applyKernel","End of Applying the kernel...",lt%verbose,myRank)   
  end subroutine gpmdcov_applyKernel

end module gpmdcov_kernel_mod
