!> Kernel related routined.
!! \brief This module will be used to compute quantities related to the Kernel.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_kernel_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_getKernel

contains

  subroutine gpmdcov_getKernel(nats)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use prg_response_mod
    real(dp), allocatable :: chargePertVect(:)
    real(dp), allocatable :: coul_pot(:)
    real(dp), allocatable :: Jacob(:)
    real(dp) :: mlsi
    integer :: nats

    if(.not. allocated(chargePertVect))allocate(chargePertVect(nats))
    if(.not. allocated(coul_pot))allocate(coul_pot(nats))

    allocate(Jacob(nats*nats))
    do i=1,nats
      write(*,*)"Constructing response for atom =",i
      chargePertVect=0.0_dp
      chargePertVect(i)=1.0_dp
      call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);

      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,chargePertVect,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);

      coul_pot = coul_pot_k + coul_pot_r

      !#ifdef DO_MPI
      !    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
        !#else
        !     do ipt = 1,gpat%TotalParts
        !#endif
        norb = syprt(ipt)%estr%norbs

        if(.not.allocated(syprt(ipt)%estr%coul_pot_k))then
          allocate(syprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
          allocate(syprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))
        endif

        syprt(ipt)%estr%coul_pot_k = 0.0_dp
        syprt(ipt)%estr%coul_pot_r = 0.0_dp
        syprt(ipt)%net_charge = 0.0_dp

        mlsi = mls()
        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          syprt(ipt)%estr%coul_pot_k(j) = coul_pot_k(jj)
          syprt(ipt)%estr%coul_pot_r(j) = coul_pot_r(jj)
          syprt(ipt)%net_charge(j) = chargePertVect(jj)
        enddo
        write(*,*)"Time for coulomb =",mls() - mlsi

        !Get part(j)\%hamPert given syprt(j)\%over, and syprt(j)\%pot
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%ham)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%mat)
        call bml_copy_new(syprt(ipt)%estr%mat,syprt(ipt)%estr%ham)

        mlsi = mls()
        call gpmdcov_msI("gpmdcov_getKernel","Entering prg_get_hscf to construct perturbative ham ...",lt%verbose,myRank)
        call prg_get_hscf(syprt(ipt)%estr%mat,syprt(ipt)%estr%over,syprt(ipt)%estr%ham,syprt(ipt)%spindex,&
             syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
             syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)
        write(*,*)"Time for hscf construction =",mls() - mlsi

        call prg_orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
                & lt%threshold,lt%bml_type,lt%verbose)

        call prg_toEigenspace(syprt(ipt)%estr%oham,syprt(ipt)%estr%ham,syprt(ipt)%estr%evects,lt%threshold,lt%verbose)
 !       call bml_print_matrix("hamPert",syprt(ipt)%estr%ham,0,10,0,10)

        mlsi = mls()
        call prg_canon_response(syprt(ipt)%estr%mat,syprt(ipt)%estr%ham,nocc,beta,syprt(ipt)%estr%evects,&
                &syprt(ipt)%estr%aux(1,:),ef,12,norb)
        write(*,*)"Time for hscf canon =",mls() - mlsi

        call prg_get_charges(syprt(ipt)%estr%mat, syprt(ipt)%estr%over, syprt(ipt)%estr%hindex, syprt(ipt)%net_charge, tb%numel,&
         syprt(ipt)%spindex, lt%mdim, lt%threshold)

        !Constructing 1D matrix J to prepare for all-gather
        do j=1,gpat%sgraph(ipt)%llsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          Jacob(i*(nats - 1) + jj) = syprt(ipt)%net_charge(j)
        enddo

      enddo

    enddo

#ifdef DO_MPI
  if (getNRanks() .gt. 1) then
    call prg_sumRealReduceN(Jacob, nats*nats)
  endif
#endif

   write(*,*)Jacob
      call prg_wait
      stop

  end subroutine gpmdcov_getKernel

end module gpmdcov_kernel_mod
