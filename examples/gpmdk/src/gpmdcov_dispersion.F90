!> Apply dispersion corrections 
!!  
!! \author C. Negre 

module gpmdcov_dispersion_mod

  implicit none 
  integer, parameter  ::  dp1 = kind(1.0d0)

  !> Type to store the dispersion potential paramenters for every pair in the system.
  !! This type will store all the possible dispersion potential parameters
  !! Allocation: If disppot is of type disppot_type, then
  !! \verbatim disppot(nsp,nsp) \endverbatim
  !! and for any particular pair of species i,j then:
  !! \verbatim disppot(i,j)%potparams(maxparams) \endverbatim
  !! Where maxparams is the maximum number of different parameters.
  !! - Example: If we want to know the parameter number 2 for the pair i,j we can do:
  !! \verbatim  disppot(i,j)%potparams(2)
  type, public :: disppot_type
    real(dp1),allocatable :: potparams(:)
  end type disppot_type

  public :: gpmdcov_read_disp_params, gpmdcov_get_disp_contrib

contains

  !> Read parameters for dispersion corrections 
  !! 
  !! \param fname_in File name 
  !! \param parampath Parameters file path
  !! \param splist Species list. Runs through all the species index.
  !! \return disppot_inout Structuee with the dispersion pair potential.
  !! Allocation: If disppot_inout is of type disppot_type, then
  !! \verbatim disppot_inout(nsp,nsp) \endverbatim
  !! and for any particular pair of species i,j then:
  !! \verbatim disppot_inout(i,j)%params(maxparams) \endverbatim
  !! \param verbose Verbose level
  !! \param rank MPI rank id. If application if not MPI enabled use 1 instead.
  !! 
  subroutine gpmdcov_read_disp_params(fname_in,parampath,splist,disppot_inout,verbose,rank)
    use gpmdcov_writeout_mod
    use prg_openfiles_mod

    implicit none 
    character(*), intent(in) :: fname_in, parampath
    character(200)  :: fname 
    character(2), intent(in) ::  splist(:)
    type(disppot_type), allocatable, intent(inout) :: disppot_inout(:,:)
    real(dp1), allocatable :: dispCoeff(:)
    integer, intent(in) :: verbose, rank
    integer :: funit, npairs
    integer :: ii,jj,kk,nsp
    character(1) :: dummyc
    character(2) :: ele1,ele2
    character(200) :: msg


    nsp = size(splist,dim=1)

    call gpmdcov_status_message("gpmdcov_read_disp_params","Reading dispersion parameters ...",verbose,1,rank)

    fname = adjustl(trim(parampath))//"/"//adjustl(trim(fname_in))

    call prg_open_file(funit, adjustl(trim(fname)))

    msg = "Using paramters at: "//fname
    call gpmdcov_info_message("gpmdcov_read_disp_params",adjustl(trim(msg)),verbose,2,rank)

    read(funit,*)dummyc,npairs
    read(funit,*)dummyc

    if(allocated(disppot_inout)) deallocate(disppot_inout)

    allocate(disppot_inout(nsp,nsp))

    allocate(DispCoeff(3))


    do jj=1,nsp
      do kk=1,nsp
        allocate(disppot_inout(jj,kk)%potparams(3))
        disppot_inout(jj,kk)%potparams = 0.0_dp1
      enddo
    enddo

    do ii=1,npairs
      read(funit,*)ele1,ele2,(dispCoeff(jj),jj=1,3)
      write(*,*)ele1,ele2,(dispCoeff(jj),jj=1,3)
      do jj=1,nsp
        do kk=1,nsp
          if(trim(adjustl(ele1)).eq.trim(adjustl(spList(jj))))then
            if(trim(adjustl(ele2)).eq.trim(adjustl(spList(kk))))then
              disppot_inout(jj,kk)%potparams(:) = DispCoeff(:)
              disppot_inout(kk,jj)%potparams(:) = DispCoeff(:)
            endif
          endif
        enddo
      enddo
    enddo

    close(funit)


    if(rank == 1)then 
      write(*,*)"Using the following dispersion potential parameters:"
      do jj=1,nsp
        do kk=1,nsp
          write(*,*)splist(jj),splist(kk),disppot_inout(jj,kk)%potparams(:)
        enddo
      enddo
    endif

  end subroutine gpmdcov_read_disp_params 

  !> Compute dispersion contribution 
  !! 
  subroutine gpmdcov_get_disp_contrib(coords,lattice_vectors,nnIx,nnIy,&
         nnIz,nrnnlist,nnType,spindex,disppot,DispForces,EDisp,verbose,myRank)
    use gpmdcov_writeout_mod
    implicit none
    integer                              ::  i, ii, j, jj
    integer                              ::  nats, nni
    integer                              ::  nr_shift_X, nr_shift_Y, nr_shift_Z
    integer, intent(in)               ::  nnIx(:,:),nnIy(:,:),nnIz(:,:)
    integer, intent(in)                  ::  spindex(:), verbose, myrank
    integer, intent(in)                  ::  nrnnlist(:), nnType(:,:)
    real(dp1)                             ::  DC(3)
    real(dp1)                             ::  Lx, Ly, Lz
    real(dp1)                             ::  DispCoef(3)
    real(dp1)                             ::  RXb, RYb, RZb
    real(dp1)                             ::  Ra(3), Rb(3), univpot
    real(dp1)                             ::  jcontrib(3),forcei(3),dR2, dr, rab(3), X
    real(dp1), allocatable, intent(inout)  ::  DispForces(:,:)
    real(dp1), intent(in)                 ::  coords(:,:), lattice_vectors(:,:)
    real(dp1), intent(inout)              ::  EDisp
    real(dp1) :: myexp, dmyexp, pot, dpot
    type(disppot_type), intent(inout)       ::  disppot(:,:)
          
    
    
    call gpmdcov_status_message("gpmdcov_get_disp_contrib","Getiing the disperion contribution ...",verbose,1,myrank)


    nats = size(coords,dim=2)
    if(.not.allocated(DispForces))then
      allocate(DispForces(3,nats))
    endif


    DispForces = 0.0_dp1

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)


    !$omp parallel do default(none) private(i) &
    !$omp private(Ra,Rb,RXb,RYb,RZb,Rab,dR,X,dR2,DC) &
    !$omp private(myexp,dmyexp,pot,dpot,jcontrib) &
    !$omp private(j,jj,ii) &
    !$omp private(DispCoef,forcei) &
    !$omp shared(nats,coords,spindex,disppot,Lx,Ly,Lz) &
    !$omp shared(DispForces,nrnnlist,nnType) &
    !$omp reduction (+:univpot)
    do i = 1, nats
      Ra(1) = coords(1,i); Ra(2) = coords(2,i); Ra(3) = coords(3,i)
      ii=spindex(i)
      forcei = 0.0_dp1
      do nni = 1,nrnnlist(i)
        j = nnType(nni,i)

        if(i.ne.j)then

          jj=spindex(j)

          DispCoef = disppot(ii,jj)%potparams;

          Rb(1) = coords(1,j)
          Rb(2) = coords(2,j)
          Rb(3) = coords(3,j)

          ! ***NOTE: The following is only valid for orthogonal unit cells

          rab(1) = modulo((Rb(1) - Ra(1) + Lx/2.0_dp1),Lx) - Lx/2.0_dp1
          rab(2) = modulo((Rb(2) - Ra(2) + Ly/2.0_dp1),Ly) - Ly/2.0_dp1
          rab(3) = modulo((Rb(3) - Ra(3) + Lz/2.0_dp1),Lz) - Lz/2.0_dp1


          dR2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
          dR = sqrt(dR2)

          DC = Rab/dR;

          x = dR - DispCoef(3)

          myexp = exp(-DispCoef(2)*x)
          dmyexp = -DispCoef(2)*myexp

          pot = DispCoef(1)*(1.0_dp1 - myexp)**2 
          dpot = -2.0_dp1*DispCoef(1)*(1.0_dp1 - myexp)*dmyexp
          jcontrib = -DC*dpot
        
          forcei = forcei - jcontrib
          univpot = univpot + pot;

        endif

      enddo
      DispForces(:,i) = forcei
    enddo
    ! $omp end parallel do

    EDisp = 0.5*univpot



  end subroutine gpmdcov_get_disp_contrib

end module gpmdcov_dispersion_mod
