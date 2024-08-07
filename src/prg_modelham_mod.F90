!> The prg_hamiltonian module.
!! \brief This module will create a model Hamiltonian for benchmarking purposes.
!!
module prg_modelham_mod

  use bml
  use prg_kernelparser_mod
  use prg_system_mod


  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  !> General ModelHam type
  type, public :: mham_type
    integer :: norbs, seed
    character(100) :: jobname
    character(100) :: bml_type
    real(dp) :: ea
    real(dp) :: eb
    real(dp) :: dab
    real(dp) :: daiaj
    real(dp) :: dbibj
    real(dp) :: dec, rcoeff
    logical :: reshuffle
  end type mham_type

  public :: prg_parse_mham, prg_twolevel_model, prg_twolevel_model3d

contains

  !> Model Ham parse.
  subroutine prg_parse_mham(mham,filename)

    implicit none
    type(mham_type), intent(inout) :: mham
    integer, parameter :: nkey_char = 2, nkey_int = 2, nkey_re = 7, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         'JobName=', 'BMLType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'GetModelHam', 'Dense' ]

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'NOrbs=', 'Seed=']
    integer :: valvector_int(nkey_int) = (/ &
         10, 100  /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'EpsilonA=', 'EpsilonB=', 'DeltaAB=','DeltaAiAj=','DeltaBiBj=','Decay=','RCoeff=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0, 0.0, -1.0, 0.0, -1.0, -100.0, 0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         'Reshuffle=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'MHAM{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    mham%JobName = valvector_char(1)

    if(valvector_char(2) == "Dense")then
      mham%bml_type = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
      mham%bml_type = BML_MATRIX_ELLPACK
    elseif(valvector_char(2) == "CSR")then
      mham%bml_type = BML_MATRIX_CSR
    elseif(valvector_char(2) == "Ellblock")then
      mham%bml_type = BML_MATRIX_ELLBLOCK
    endif

    !Integers
    mham%norbs = valvector_int(1)
    mham%seed = valvector_int(2)

    !Reals
    mham%ea = valvector_re(1)
    mham%eb = valvector_re(2)
    mham%dab = valvector_re(3)
    mham%daiaj = valvector_re(4)
    mham%dbibj = valvector_re(5)
    mham%dec = valvector_re(6)
    mham%rcoeff = valvector_re(7)

    !Logicals
    mham%reshuffle = valvector_log(1)

  end subroutine prg_parse_mham

  !> Construct a two-level model Hamiltonian
  !!
  !! \param ea First onsite energy
  !! \param eb Second onsite energy
  !! \param dab Onsite Hamiltonian element
  !! \param daiaj Intersite first level Hamiltonian elements
  !! \param dbibj Intersite second level Hamiltonian elements
  !! \param dec Decay constant
  !! \param rcoeff Random coefficient
  !! \param reshuffle If rows needs to be reshuffled
  !! \param seed Random seed
  !! \param h_bml Output hamiltonian matrix
  !! \param verbose Verbosity level
  subroutine prg_twolevel_model(ea, eb, dab, daiaj, dbibj, dec, rcoeff, reshuffle, &
       & seed, h_bml, verbose)
    real(dp), intent(in) :: ea, eb, dab, daiaj, dbibj, rcoeff
    integer, intent(in) :: verbose, seed
    integer, allocatable :: seedin(:)
    logical, intent(in) :: reshuffle
    type(bml_matrix_t),intent(inout) :: h_bml
    real(dp), allocatable :: diagonal(:), row(:), rowi(:), rowj(:)
    type(bml_matrix_t) :: ht_bml
    integer :: norbs, i, j, ssize
    real(dp) :: dec, dist, ran

    norbs = bml_get_N(h_bml)
    allocate(diagonal(norbs))
    allocate(row(norbs))

    call random_seed()
    call random_seed(size=ssize)
    allocate(seedin(ssize))
    seedin = seed
    call random_seed(PUT=seedin)

    do i=1,norbs
      if(mod(i,2) == 0)then
        call random_number(ran)
        diagonal(i) = ea + rcoeff*(2.0_dp*ran - 1.0_dp)
      else
        call random_number(ran)
        diagonal(i) = eb + rcoeff*(2.0_dp*ran - 1.0_dp)
      endif
    enddo

    do i=1,norbs
      do j=1,norbs
        if(abs(real(i-j,dp)) <= norbs/2.0d0) then
          dist = max(abs(real(i-j,dp))- 2.0_dp,0.0_dp)
        else
          dist = max((-abs(real(i-j,dp))+norbs) - 2.0_dp,0.0_dp)
        endif
        !A-A type
        if((mod(i,2) .ne. 0) .and. (mod(j,2) .ne. 0))then
          call random_number(ran)
          row(j) = (daiaj + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
          !A-B type
        elseif((mod(i,2) == 0) .and. (mod(j,2) == 0))then
          call random_number(ran)
          row(j) = (dbibj + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
          !B-B type
        else
          call random_number(ran)
          row(j) = (dab + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
        endif
        ! write(*,*)i,j,row(j),mod(i,2),mod(j,2),abs(real(i-j,dp)+norbs),abs(real(i-j,dp)+norbs)

      enddo
      call bml_set_row(h_bml,i,row)
    enddo

    call bml_set_diagonal(h_bml,diagonal)

    !Symmetrization
    call bml_copy_new(h_bml,ht_bml)
    call bml_transpose(h_bml,ht_bml)
    if(verbose.gt.0)then
      call bml_print_matrix("h_bml",h_bml,0,10,0,10)
      call bml_print_matrix("ht_bml",ht_bml,0,10,0,10)
    endif
    call bml_add(h_bml,ht_bml,0.5d0,0.5d0,0.0d0)

    call bml_deallocate(ht_bml)

    if(reshuffle)then
      allocate(rowj(norbs))
      allocate(rowi(norbs))
      do i=1,norbs
        call random_number(ran)
        j = int(floor(ran*norbs+1))
        call bml_get_row(h_bml,i,rowi)
        call bml_get_row(h_bml,j,rowj)
        call bml_set_row(h_bml,i,rowj)
        call bml_set_row(h_bml,j,rowi)
      enddo
      deallocate(rowi)
      deallocate(rowj)
    endif

  end subroutine prg_twolevel_model

  !> Construct a two-level model Hamiltonian
  !!
  !! \param ea First onsite energy
  !! \param eb Second onsite energy
  !! \param dab Onsite Hamiltonian element
  !! \param daiaj Intersite first level Hamiltonian elements
  !! \param dbibj Intersite second level Hamiltonian elements
  !! \param dec Decay constant
  !! \param rcoeff Random coefficient
  !! \param reshuffle If rows needs to be reshuffled
  !! \param seed Random seed
  !! \param h_bml Output hamiltonian matrix
  !! \param verbose Verbosity level
  subroutine prg_twolevel_model3d(ea, eb, dab, daiaj, dbibj, dec, rcoeff, reshuffle, &
       & seed, h_bml, verbose)
    real(dp), intent(in) :: ea, eb, dab, daiaj, dbibj, rcoeff
    integer, intent(in) :: verbose, seed
    integer, allocatable :: seedin(:)
    logical, intent(in) :: reshuffle
    type(bml_matrix_t),intent(inout) :: h_bml
    real(dp), allocatable :: diagonal(:), row1(:), row2(:), rowi(:), rowj(:)
    type(bml_matrix_t) :: ht_bml
    integer :: norbs, i, j, ssize, i1, j1, k1, i2, j2, k2, n1d
    real(dp) :: dec, dist, di, dj, dk, ran, tmp

    norbs = bml_get_N(h_bml)

    allocate(diagonal(norbs))
    allocate(row1(norbs))
    allocate(row2(norbs))

    tmp = norbs
    tmp = tmp*0.5
    tmp = tmp**(1./3.)
    n1d = int(tmp)
    print*,'n1d = ',n1d

    call random_seed()
    call random_seed(size=ssize)
    allocate(seedin(ssize))
    seedin = seed
    call random_seed(PUT=seedin)

    do i=1,norbs
      if(mod(i,2) == 0)then
        call random_number(ran)
        diagonal(i) = ea + rcoeff*(2.0_dp*ran - 1.0_dp)
      else
        call random_number(ran)
        diagonal(i) = eb + rcoeff*(2.0_dp*ran - 1.0_dp)
      endif
    enddo

    do i1=1,n1d
      do j1=1,n1d
        do k1=1,n1d
          i = (i1-1) + n1d*(j1-1) + n1d*n1d*(k1-1)
          do i2=1,n1d
            do j2=1,n1d
              do k2=1,n1d
                j = (i2-1) + n1d*(j2-1) + n1d*n1d*(k2-1)
                if(abs(real(i1-i2,dp)) <= n1d/2.0d0) then
                  di = max(abs(real(i1-i2,dp))- 2.0_dp,0.0_dp)
                else
                  di = max((-abs(real(i1-i2,dp))+n1d)- 2.0_dp,0.0_dp)
                endif
                if(abs(real(j1-j2,dp)) <= n1d/2.0d0) then
                  dj = max(abs(real(j1-j2,dp))- 2.0_dp,0.0_dp)
                else
                  dj = max((-abs(real(j1-j2,dp))+n1d)- 2.0_dp,0.0_dp)
                endif
                if(abs(real(k1-k2,dp)) <= n1d/2.0d0) then
                  dk = max(abs(real(k1-k2,dp))- 2.0_dp,0.0_dp)
                else
                  dk = max((-abs(real(k1-k2,dp))+n1d)- 2.0_dp,0.0_dp)
                endif
                dist = dsqrt(di+di+dj*dj+dk*dk)

                !assign matrix elements for all interactions between two atoms
                !A-A type
                call random_number(ran)
                row1(2*j+1) = (daiaj + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
                !A-B type
                call random_number(ran)
                row1(2*j+2) = (dab + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
                row2(2*j+1) = (dab + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
                !B-B type
                call random_number(ran)
                row2(2*j+2) = (dbibj + rcoeff*(2.0_dp*ran - 1.0_dp))*exp(dec*dist)
                ! write(*,*)i,j,row(j),mod(i,2),mod(j,2),abs(real(i-j,dp)+norbs),abs(real(i-j,dp)+norbs)
              enddo
            enddo
          enddo
          call bml_set_row(h_bml,2*i+1,row1)
          call bml_set_row(h_bml,2*i+2,row2)
        enddo
      enddo
    enddo

    call bml_set_diagonal(h_bml,diagonal)

    !Symmetrization (necessary when random factor used)
    call bml_copy_new(h_bml,ht_bml)
    call bml_transpose(h_bml,ht_bml)
    if(verbose.gt.0)then
      call bml_print_matrix("h_bml",h_bml,0,10,0,10)
      call bml_print_matrix("ht_bml",ht_bml,0,10,0,10)
    endif
    call bml_add(h_bml,ht_bml,0.5d0,0.5d0,0.0d0)

    call bml_deallocate(ht_bml)

    if(reshuffle)then
      allocate(rowj(norbs))
      allocate(rowi(norbs))
      do i=1,norbs
        call random_number(ran)
        j = int(floor(ran*norbs+1))
        call bml_get_row(h_bml,i,rowi)
        call bml_get_row(h_bml,j,rowj)
        call bml_set_row(h_bml,i,rowj)
        call bml_set_row(h_bml,j,rowi)
      enddo
      deallocate(rowi)
      deallocate(rowj)
    endif

  end subroutine prg_twolevel_model3d

end module prg_modelham_mod
