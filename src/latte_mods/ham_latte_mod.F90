!> Hamiltonian module. 
!! \ingroup LATTE 
!! \brief Routines in this module are used to build the Hamiltonian and Svelap matrix from the system type.
!! 
module ham_latte_mod

  use bml 
  use tbparams_latte_mod
  use prg_extras_mod

  implicit none

  private 

  integer, parameter :: dp = kind(1.0d0)

  public :: get_hindex, get_hsmat, get_SKBlock

contains 

  !> Gets the Hamiltonian indices for every atom in the system.
  !! \param spindex Species index for every atom in the system.
  !! \param norbi Number of orbital for every species in the species list.
  !! \param hindex Start and end index for every atom in the system.
  !! \param norb Output parameter corresponding to the number of orbitals of the system.
  !! \note norb is always greater or equal nats.
  !! \note spindex can be gather using the system type: 
  !! \verbatim spindex(:) = system%spindex(:) \endverbatim
  !! \note norbi can be gather from the tbparams type as it is a property 
  !! that depends strictly on the parametrization;
  !! \verbatim norbi(:) = tbparams%norbi(:) \endverbatim
  !! \todo add verbosity
  !! 
  subroutine get_hindex(spindex,norbi,hindex,norb,verbose)
    implicit none
    integer, intent(in) :: spindex(:)
    integer, intent(in) :: norbi(:)
    integer, optional, intent(in) :: verbose
    integer, intent(inout) :: norb
    integer :: nats, cnt, i
    integer, allocatable, intent(inout) :: hindex(:,:)


    nats = size(spindex,dim=1)
    cnt = 1
    
    if(present(verbose).and.verbose >= 1)then 
      write(*,*)""; write(*,*)"In get hindex ..."
    endif

    if(.not.allocated(hindex))then 
      allocate(hindex(2,nats))
    endif

    cnt = 1
    do i = 1,nats
      hindex(1,i) = cnt
      cnt = cnt + norbi(spindex(i)) 
      hindex(2,i) = cnt - 1
    enddo

    norb = cnt-1;

    if(present(verbose).and.verbose >= 1)then 
      write(*,*)"Number of orbitals =",norb
    endif
    
  end subroutine get_hindex  

  !> Constructs Hamiltonian and Overlap Matrix
  !! \brief Construction of the Hamiltonian and Overlap matrix. 
  !! \param ham_bml Hamiltonian in bml format. 
  !! \param over_bml Overlap in bml format. 
  !! \param coordinate Coordinates of the system. This can be getter from the system type: 
  !! \verbatim system%coordinate \endverbatim
  !! \param lattice_vector Lattice vectors for the system.
  !! \param spindex Species indices (see system_type).
  !! \param norbi Number of orbital for every species.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex)
  !! \param onsitesH Onsites energies for every pair of equal type. Allocation:
  !! \verbatim onsitesH(maxints,nsp) \endverbatim
  !! \param onsitesS Same as the Hamiltonian but for the overlap matrix. 
  !! In this case the "onsite" elements are equal to 1.0. This was done to maintain
  !! a consistency and have the possibility of generalize the SK block construction.
  !! elements are or the same type.   
  !! \param intPairsH,intPairsS See in intPairs_type
  !! \param threshold Threshold value for matrix elements. 
  subroutine get_hsmat(ham_bml,over_bml,coordinate,lattice_vector,spindex,&
      norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,threshold)
    implicit none
    character(20)                       ::  bml_type, bml_dmode
    integer                             ::  dimi, dimj, i, ii
    integer                             ::  j, jj, nats, norb
    integer, intent(in)                 ::  hindex(:,:), norbi(:), spindex(:)
    integer                             ::  maxnorbi
    real(dp)                            ::  ra(3), rb(3)
    real(dp), allocatable               ::  block(:,:,:), ham(:,:), over(:,:)
    real(dp), intent(in)                ::  coordinate(:,:), lattice_vector(:,:), onsitesH(:,:), onsitesS(:,:)
    real(dp), intent(in)                ::  threshold
    type(bml_matrix_t), intent(inout)   ::  ham_bml, over_bml
    type(intpairs_type), intent(in)     ::  intPairsH(:,:), intPairsS(:,:)

    nats = size(spindex,dim=1)
    norb = bml_get_N(ham_bml)
    bml_type = bml_get_type(ham_bml)
    bml_dmode = bml_get_distribution_mode(ham_bml)

!     if(bml_get_N(ham_bml).LE.0) then  !Carefull we need to clean S and H before rebuilding them!!!
      call bml_deallocate(ham_bml)
      call bml_deallocate(over_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,ham_bml, &
        bml_dmode)    
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,over_bml, &
        bml_dmode)    
!     endif

    maxnorbi = maxval(norbi)

    if (.not.allocated(block)) then
       allocate(block(maxnorbi,maxnorbi,nats))
    endif

   !$omp parallel do default(none) firstprivate(j) &
   !$omp private(ra,rb,dimi,dimj,ii,jj) &
   !$omp shared(nats,coordinate,hindex,spindex, intPairsS,intPairsH,threshold,lattice_vector,norbi,onsitesH,onsitesS,ham_bml,over_bml) &
   !$omp shared(block)
    do i = 1, nats
      ra(:) = coordinate(:,i)
      dimi = hindex(2,i)-hindex(1,i)+1
      do j = 1, nats
        rb(:) = coordinate(:,j) 
        dimj = hindex(2,j)-hindex(1,j)+1
        !Hamiltonian block for a-b atom pair
        call get_SKBlock(spindex(i),spindex(j),coordinate(:,i),&
        coordinate(:,j),lattice_vector,norbi,&
        onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,block,i)

        do jj=1,dimj
          do ii=1,dimi
!               ham(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = block(ii,jj)
             if(abs(block(ii,jj,i)).gt.threshold)then 
               call bml_set_element_new(ham_bml,hindex(1,i)-1+ii,&
                hindex(1,j)-1+jj,block(ii,jj,i))
             endif
          enddo
        enddo  

        call get_SKBlock(spindex(i),spindex(j),coordinate(:,i),&
          coordinate(:,j),lattice_vector,norbi,&
          onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,block,i)

        do jj=1,dimj
          do ii=1,dimi
!             over(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = block(ii,jj)
             if(abs(block(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(over_bml,hindex(1,i)-1+ii,&
                hindex(1,j)-1+jj,block(ii,jj,i))
             endif
          enddo
        enddo     

        !  write(*,*)spindex(i),spindex(j)   
        !  write(*,'(100F10.5)')intPairsS(spindex(i),spindex(j))%intParams
        ! call print_matrix("block",block(:,:,i),1,4,1,4)

      enddo
    enddo
   !$omp end parallel do

!     bml_type=bml_get_type(ham_bml) !Get the bml type
!     call bml_convert_from_dense(bml_type,over,over_bml,threshold,norb) !Dense to dense_bml
!     call bml_convert_from_dense(bml_type,ham,ham_bml,threshold,norb) !Dense to dense_bml

!     call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
!     call bml_print_matrix("over_bml",over_bml,0,6,0,6)

  end subroutine get_hsmat

  !> Function to calculate the bond integral for a given distance 
  !! and coefficients set. 
  !! \param dr distance between atoms.
  !! \param f parameters (coefficients) for the bond integral.  
  real(dp) function bondIntegral(dr,f)
    implicit none 
    real(dp) :: rmod
    real(dp) :: polynom
    real(dp) :: rminusr1
    real(dp) :: x
    real(dp), intent(in) :: dr
    real(dp), intent(in) :: f(16)

    if(dr <= f(7))then 
      rmod = dr - f(6);
      polynom = rmod*(f(2) + rmod*(f(3) + rmod*(f(4) + f(5)*rmod)));
      x = exp(polynom);
    elseif(dr > f(7).and.dr < f(8))then 
      rminusr1 = dr - f(7)
      x = f(9) + rminusr1*(f(10) + rminusr1*(f(11) + rminusr1*(f(12) + rminusr1*(f(13) + rminusr1*f(14)))))
    else
      x = 0
    end if
    bondintegral = f(1)*x

  end function bondintegral

  !> Standard Slater-Koster sp-parameterization for an atomic block between a pair of atoms
  !! \param sp1 Species index for atom 1. This can be obtained from the 
  !! system type as following:
  !! \verbatim sp1 = system%spindex(atom1) \endverbatim
  !! \param sp2 Species index for atom 2. 
  !! \param coorda Coordinates for atom 1.
  !! \param coordb Coordinates for atom 2.
  !! \param lattice_vectors Lattice vectors for the system. This can be obtained from 
  !! the system type as following: 
  !! \verbatim lattice_vectors = system%lattice_vectors \endverbatim
  !! \param norbi Number of orbitals for every species in the system. This can be obtained from 
  !! the tbparams type as following: 
  !! \verbatim norbi = tbparams%norbi \endverbatim
  !! \param onsites Onsites energies for every pair of equal type. Two different variants
  !! onsitesH and onsitesS will be used as inputs (see get_hsmat routine) Allocation:
  !! \verbatim onsites(maxints,nsp) \endverbatim
  !! \param intParams See intpairs_type.
  !! \param block Output parameter SK block. 
  !! \param atnum Input atom number
  subroutine get_SKBlock(sp1,sp2,coorda,coordb,lattice_vectors&
      ,norbi,onsites,intParams,intParamsr,block,atnum)
    implicit none
    integer                              ::  dimi, dimj, i, nr_shift_X
    integer                              ::  nr_shift_Y, nr_shift_Z
    integer, intent(in)                  ::  norbi(:), sp1, sp2, atnum
    real(dp)                             ::  HPPP, HPPS, HSPS, HSPSR, HSSS
    real(dp)                             ::  L, LBox(3), M, N
    real(dp)                             ::  PPSMPP, PXPX, PXPY, PXPZ
    real(dp)                             ::  PYPX, PYPY, PYPZ, PZPX
    real(dp)                             ::  PZPY, PZPZ, dr, ra(3)
    real(dp)                             ::  rab(3), rb(3), rxb, ryb
    real(dp)                             ::  rzb
    real(dp), allocatable, intent(inout)  ::  block(:,:,:)
    real(dp), intent(in)                 ::  coorda(:), coordb(:), intParams(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  onsites(:,:), intParamsr(:,:)

    ra = coorda
    rb = coordb

    dimi= norbi(sp1)
    dimj= norbi(sp2)

!     write(*,*)atom_type_a, atom_type_b,dimi,dimj

!!    if(allocated(block))then 
!!      deallocate(block)
!!    endif    

!!    allocate(block(dimi,dimj))
    block(:,:,atnum)=0.0_dp

!     call write_matrix_to_screen("block",block,size(block,dim=1),size(block,dim=2))

    RXb = Rb(1); RYb = Rb(2); RZb = Rb(3)

! For cubic lattice
!     LBox(1) = lattice_vectors(1,1)
!     LBox(2) = lattice_vectors(2,2)
!     LBox(3) = lattice_vectors(3,3)   

    !Periodic BC shifts in X, Y and Z. Costs a lot extra!
    do nr_shift_x = -1,1  
      do nr_shift_y = -1,1
        do nr_shift_z = -1,1

              rb(1) = RXb + nr_shift_x*lattice_vectors(1,1) ! shifts for pbc
              ! rb(1) = rb(1) + nr_shift_y*lattice_vectors(2,1) ! shifts for pbc            
              ! rb(1) = rb(1) + nr_shift_z*lattice_vectors(3,1) ! shifts for pbc            

              rb(2) = RYb + nr_shift_y*lattice_vectors(2,2) ! shifts for pbc
              ! rb(2) = rb(2) + nr_shift_x*lattice_vectors(1,2) ! shifts for pbc            
              ! rb(2) = rb(2) + nr_shift_z*lattice_vectors(3,2) ! shifts for pbc            

              rb(3) = RZb + nr_shift_z*lattice_vectors(3,3) ! shifts for pbc
              ! rb(3) = rb(3) + nr_shift_y*lattice_vectors(2,3) ! shifts for pbc            
              ! rb(3) = rb(3) + nr_shift_x*lattice_vectors(1,3) ! shifts for pbc            

          Rab = Rb-Ra;  ! OBS b - a !!!
          dR = sqrt(Rab(1)**2+ Rab(2)**2+ Rab(3)**2)

          if(dR.lt.6.5_dp)then
          if(dR .LT.1e-12)then !same position and thus the same type sp1 = sp2
            do i=1,dimi
              block(i,i,atnum) = onsites(i,sp1)
            enddo
          else

            L = Rab(1)/dR;  !Direction cosines
            M = Rab(2)/dR;
            N = Rab(3)/dR;
            if(dimi == dimj.and.dimi == 1)then        !s-s  overlap 1 x 1 block
              HSSS = BondIntegral(dR,intParams(:,1))  !Calculate the s-s bond integral
              block(1,1,atnum) = block(1,1,atnum) + HSSS
            elseif(dimi < dimj.and.dimi == 1)then    !s-sp overlap 1 x 4 block
              HSSS = BondIntegral(dR,intParams(:,1))
              block(1,1,atnum) = block(1,1,atnum) + HSSS        
              HSPS = BondIntegral(dR,intParams(:,2))
              block(1,2,atnum) = block(1,2,atnum) + L*HSPS
              block(1,3,atnum) = block(1,3,atnum) + M*HSPS
              block(1,4,atnum) = block(1,4,atnum) + N*HSPS
            elseif(dimi > dimj.and.dimj == 1)then ! sp-s overlap 4 x 1 block
              HSSS = BondIntegral(dR,intParams(:,1))
              block(1,1,atnum) = block(1,1,atnum) + HSSS
              HSPS = BondIntegral(dR,intParams(:,2))
              block(2,1,atnum) = block(2,1,atnum) - L*HSPS
              block(3,1,atnum) = block(3,1,atnum) - M*HSPS
              block(4,1,atnum) = block(4,1,atnum) - N*HSPS
            elseif(dimi == dimj.and.dimj == 4)then !sp-sp overlap
              HSSS = BondIntegral(dR,intParams(:,1))
              HSPS = BondIntegral(dR,intParams(:,2))
              HSPSR = BondIntegral(dR,intParamsr(:,2))                            
              HPPS = BondIntegral(dR,intParams(:,3))
              HPPP = BondIntegral(dR,intParams(:,4))
              PPSMPP = HPPS - HPPP
              PXPX = HPPP + L*L*PPSMPP
              PXPY = L*M*PPSMPP
              PXPZ = L*N*PPSMPP
              PYPX = M*L*PPSMPP
              PYPY = HPPP + M*M*PPSMPP
              PYPZ = M*N*PPSMPP
              PZPX = N*L*PPSMPP
              PZPY = N*M*PPSMPP
              PZPZ = HPPP + N*N*PPSMPP
              block(1,1,atnum) = block(1,1,atnum) + HSSS
              block(1,2,atnum) = block(1,2,atnum) + L*HSPS
              block(1,3,atnum) = block(1,3,atnum) + M*HSPS
              block(1,4,atnum) = block(1,4,atnum) + N*HSPS
              block(2,1,atnum) = block(2,1,atnum) - L*HSPSR  !Change spindex
              block(2,2,atnum) = block(2,2,atnum) + PXPX
              block(2,3,atnum) = block(2,3,atnum) + PXPY
              block(2,4,atnum) = block(2,4,atnum) + PXPZ
              block(3,1,atnum) = block(3,1,atnum) - M*HSPSR  !Change spindex
              block(3,2,atnum) = block(3,2,atnum) + PYPX
              block(3,3,atnum) = block(3,3,atnum) + PYPY
              block(3,4,atnum) = block(3,4,atnum) + PYPZ
              block(4,1,atnum) = block(4,1,atnum) - N*HSPSR  !Change spindex
              block(4,2,atnum) = block(4,2,atnum) + PZPX
              block(4,3,atnum) = block(4,3,atnum) + PZPY
              block(4,4,atnum) = block(4,4,atnum) + PZPZ
            endif
          endif
        endif
        enddo
      enddo
    enddo
!            write(*,*)"block",dr,block
!            stop
  end subroutine get_SKBlock

end module ham_latte_mod
