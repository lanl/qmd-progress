!> A module to compute Extended huckel Hamiltonian and Overlap matrices.
!! \brief This module will compute H ans S from the EH parameter.
!! This code was rewritten from: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwj6ponmzejLAhUjsoMKHT3xADAQFggcMAA&url=https%3A%2F%2Fwww.researchgate.net%2Ffile.PostFileLoader.html%3Fid%3D550c2cf8cf57d7c7218b45dc%26assetKey%3DAS%253A273739621568514%25401442276019154&usg=AFQjCNGITp8x-EHwSX4O5BcpcKOp8-4FNA&sig2=1lNCl1qo0P6QjcTlIlqiuA
!! This code should not be taken as part of PROGRESS or LATTE library.
!! @ingroup EXTERNAL
!!
module huckel_latte_mod

  use prg_openfiles_mod
  use prg_ptable_mod
  use bml

  implicit none     

  private 

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: SQRT3    =    1.73205080756888
  real(dp), parameter :: SQRT6    =    2.44948974278318
  real(dp), parameter :: SQRT10   =    3.16227766016838
  real(dp), parameter :: SQRT15   =    3.87298334620742
  real(dp), parameter :: AUI      =    1.889644746

  public :: get_hshuckel

  !> The single orb huckel type
  type :: single_orb  
    integer :: orb_of_e
    real(dp) :: VSIP
    real(dp) :: expo
    real(dp) :: exp2(2)
    real(dp) :: coef2(2)
  end type 

  !> Huckel atom type
  type :: atom_parameter 
    character(3) :: symbol
    integer :: valence_electron
    type(single_orb) :: orb(4)
  end type 

  !> atom type  
  type :: atom
    integer atomtype;
    real(dp) x;
    real(dp) y;
    real(dp) z;
  end type 

contains

  !> Get the the Huckel H and S. 
  !! \param extension Extension of the file
  !! 

  subroutine get_hshuckel(ham_bml,over_bml,coordinates,spindex,spatnum&
      ,parampath,bml_type,mdim,threshold& 
      ,nsp,splist,basis,numel,onsite_energ,&
      norbi,hubbardu)
    implicit none
    character(len=*), intent(in) :: parampath
    type(atom_parameter), allocatable :: period(:)
    type(bml_matrix_t), intent(inout) :: ham_bml
    type(bml_matrix_t), intent(inout) :: over_bml
    integer, intent(inout) :: nsp
    integer, intent(in) :: spindex(:)
    integer, intent(in) :: spatnum(:)
    character(2), allocatable, intent(inout) :: splist(:)
    character(4),allocatable, intent(inout) :: basis(:)
    real(dp),allocatable, intent(inout) :: numel(:)
    real(dp),allocatable, intent(inout) :: onsite_energ(:,:)
    integer,allocatable, intent(inout) :: norbi(:)
    real(dp),allocatable, intent(inout) :: hubbardu(:)
    real(dp), intent(in) :: coordinates(:,:)
    type(atom), allocatable :: molecule(:)
    integer :: no_atoms, i,j, indexi, indexj
    integer :: noorb_i, atom_row, no_orbitals
    integer :: atom_col, noorb_j, ii, jj
    real(dp) :: delx, dely, delz, S(16,16),H(16,16)
    character(len=*), intent(in) :: bml_type
    integer, intent(inout) :: mdim
    real(dp), intent(in) :: threshold
    real(dp), allocatable :: row(:)

    write(*,*)"In get_hshuckel ..."

    allocate(period(103))

    call read_atomic_parameters(period,parampath)

    no_atoms = size(spindex,dim=1)

    nsp = size(spatnum,dim=1)

    if(.not.allocated(norbi)) allocate(norbi(nsp))
    if(.not.allocated(splist)) allocate(splist(nsp))
    if(.not.allocated(basis)) allocate(basis(nsp))
    if(.not.allocated(numel)) allocate(numel(nsp))
    if(.not.allocated(onsite_energ))allocate(onsite_energ(4,nsp))
    if(.not.allocated(hubbardu))allocate(hubbardu(nsp))
    !     if(.not.allocated(mass))allocate(mass(nsp))

    allocate(molecule(no_atoms))

    do i=1,no_atoms
      molecule(i)%atomtype=spatnum(spindex(i))
      molecule(i)%x=coordinates(1,i)
      molecule(i)%y=coordinates(2,i)
      molecule(i)%z=coordinates(3,i)
    enddo     

    !get norb and no_orbitals
    indexi=0
    do i=1, no_atoms

      indexj=1
      atom_row=molecule(i)%atomtype
      noorb_i=0
      if(period(atom_row)%orb(1)%orb_of_e.ne.0) noorb_i=noorb_i+1
      if(period(atom_row)%orb(2)%orb_of_e.ne.0) noorb_i=noorb_i+3
      if(period(atom_row)%orb(3)%orb_of_e.ne.0) noorb_i=noorb_i+5
      if(period(atom_row)%orb(4)%orb_of_e.ne.0) noorb_i=noorb_i+7
      indexi=indexi+noorb_i

    enddo
    no_orbitals = indexi

    ! /*****************CALCULATIONS*******************/

    !Allocate bml's  
    if(mdim.lt.0)mdim=no_orbitals
    if(bml_get_N(ham_bml).le.0)then 
      call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,no_orbitals,ham_bml)
      call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,no_orbitals,over_bml)
    endif

    indexi=0;

    do i=1, no_atoms

      indexj=1
      atom_row=molecule(i)%atomtype
      noorb_i=0
      if(period(atom_row)%orb(1)%orb_of_e.ne.0) noorb_i=noorb_i+1
      if(period(atom_row)%orb(2)%orb_of_e.ne.0) noorb_i=noorb_i+3
      if(period(atom_row)%orb(3)%orb_of_e.ne.0) noorb_i=noorb_i+5
      if(period(atom_row)%orb(4)%orb_of_e.ne.0) noorb_i=noorb_i+7

      do j=1,no_atoms

        delx=molecule(j)%x-molecule(i)%x
        dely=molecule(j)%y-molecule(i)%y
        delz=molecule(j)%z-molecule(i)%z
        atom_col=molecule(j)%atomtype
        noorb_j=0
        if(period(atom_col)%orb(1)%orb_of_e.ne.0) noorb_j=noorb_j+1
        if(period(atom_col)%orb(2)%orb_of_e.ne.0) noorb_j=noorb_j+3
        if(period(atom_col)%orb(3)%orb_of_e.ne.0) noorb_j=noorb_j+5
        if(period(atom_col)%orb(4)%orb_of_e.ne.0) noorb_j=noorb_j+7

        S=0
        H=0

        if(i==j)then 

          do ii=1,16
            S(ii,ii)=1;
            if(ii==1) H(ii,ii)=period(atom_col)%orb(1)%VSIP
            if(4>=ii.and.ii>1) H(ii,ii)=period(atom_col)%orb(2)%VSIP;
            if(9>=ii.and.ii>4) H(ii,ii)=period(atom_col)%orb(3)%VSIP;
            if(16>=ii.and.ii>9) H(ii,ii)=period(atom_col)%orb(4)%VSIP;
          enddo
        else

          call overlap(period,atom_row,atom_col,delx,dely,delz,S,H);

          !           if(atom_row==6.and.atom_col==16)then 
          !         do ii=1, noorb_i
          !           do jj=1,noorb_j
          !             write(*,'(2I,F7.4)')ii,jj, S(ii,jj)
          !           enddo   
          !         enddo
          !           stop
          !           
          !           endif 
        endif 

        do ii=1, noorb_i
          do jj=1,noorb_j
            call bml_set_element_new(ham_bml,indexi+ii,indexj+jj-1,H(ii,jj))
            call bml_set_element_new(over_bml,indexi+ii,indexj+jj-1,S(ii,jj))
          enddo   
        enddo

        indexj=indexj+noorb_j

      enddo

      indexi=indexi+noorb_i

    enddo

    no_orbitals=indexi


    !Get the number of orbital for every species.
    do i=1,nsp
      noorb_i=0
      atom_row = spatnum(i)
      if(period(atom_row)%orb(1)%orb_of_e.ne.0) noorb_i=noorb_i+1
      if(period(atom_row)%orb(2)%orb_of_e.ne.0) noorb_i=noorb_i+3
      if(period(atom_row)%orb(3)%orb_of_e.ne.0) noorb_i=noorb_i+5
      if(period(atom_row)%orb(4)%orb_of_e.ne.0) noorb_i=noorb_i+7
      norbi(i)=noorb_i
      if(norbi(i).eq.1)basis(i)="s"
      if(norbi(i).eq.4)basis(i)="sp"
      if(norbi(i).eq.9)basis(i)="spd"
      if(norbi(i).eq.16)basis(i)="spdf"
      splist(i) = period(atom_row)%symbol
      numel(i) = period(atom_row)%valence_electron
      onsite_energ(1,i) = period(atom_row)%orb(1)%VSIP
      onsite_energ(2,i) = period(atom_row)%orb(2)%VSIP
      onsite_energ(3,i) = period(atom_row)%orb(3)%VSIP
      onsite_energ(4,i) = period(atom_row)%orb(4)%VSIP
      hubbardu(i) = element_ip(atom_row)-element_ea(atom_row)
    enddo

    call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)

  end subroutine get_hshuckel

  subroutine read_atomic_parameters(period,parampath)
    character(len=*), intent(in) :: parampath
    character(100) :: filename
    integer :: io_unit, tsize
    integer :: i,j
    type(atom_parameter), intent(inout) :: period(:)
    real(dp) :: table(24)

    filename = adjustl(trim(parampath))//"/parameters.dat"

    call open_file_to_read(io_unit,filename)

    tsize = 103 

    do i=1,tsize
      period(i)%symbol="";
      read(io_unit,*)period(i)%symbol,period(i)%valence_electron&
        ,(period(i)%orb(j)%orb_of_e,j=1,4),(table(j),j=1,24)

      do j=1,4      
        period(i)%orb(j)%VSIP = table((j-1)*6+1)
        period(i)%orb(j)%expo = table((j-1)*6+2)
        period(i)%orb(j)%exp2(1) = table((j-1)*6+3)
        period(i)%orb(j)%exp2(2) = table((j-1)*6+4)
        period(i)%orb(j)%coef2(1) = table((j-1)*6+5)
        period(i)%orb(j)%coef2(2) = table((j-1)*6+6)

        period(6)%orb(j)%coef2(2)=0

        if(period(i)%orb(j)%coef2(2)==0 .or. j<3)then
          period(i)%orb(j)%exp2(1)=period(i)%orb(j)%expo
          period(i)%orb(j)%coef2(1)=1
          period(i)%orb(j)%coef2(2)=0       
        endif          

      enddo           

    enddo
    close(io_unit)

  end subroutine read_atomic_parameters

  subroutine overlap(period,atom_row,atom_col,delx,dely,delz,S,H)

    integer :: atom_row,atom_col
    real(dp) :: delx,dely,delz
    real(dp), intent(inout) :: S(:,:),H(:,:)
    type(atom_parameter), intent(inout) :: period(:)      
    real(dp) :: rt2,r,t,ca,cb,sa,sb,ca2,sa2,cb2,sb2,cbsb
    real(dp) :: casa,cb2sb2,s2b,sa3,c2b,s3b,c3b,s2a,c2a,VSIP_row,VSIP_col
    real(dp) ::  pt(9),dt(25),ft(49),ptr(9),dtr(25),ftr(49)
    real(dp) :: sigma,pi,delta,phi
    integer :: j,k,n1,n2,nsrow,nprow,ndrow,nfrow,nscol,npcol,ndcol,nfcol

    ptr=pt-1
    dtr=dt-1
    ftr=ft-1

    rt2=delx*delx+dely*dely
    r=sqrt(rt2+delz*delz)

    if(rt2 <= 1e-10)then
      cb=1.0
      sb=0.0
      sa=0.0
    else
      t=sqrt(rt2)
      cb=delx/t
      sb=dely/t
      sa=t/r
    endif

    ca=delz/r;

    ptr(1)=sa*cb;
    ptr(2)=sa*sb;
    ptr(3)=ca;
    ptr(4)=ca*cb;
    ptr(5)=ca*sb;
    ptr(6)=-sa;
    ptr(7)=-sb;
    ptr(8)=cb;
    ptr(9)=0.0;

    if(period(atom_row)%orb(3)%orb_of_e+period(atom_col)%orb(3)%orb_of_e>0) then

      ca2=ca*ca
      sa2=sa*sa
      cb2=cb*cb
      sb2=sb*sb
      cbsb=cb*sb
      casa=ca*sa
      cb2sb2=cb2-sb2
      dtr(1)=SQRT3*0.5*sa2*cb2sb2
      dtr(2)=1.0-1.5*sa2
      dtr(3)=SQRT3*cbsb*sa2
      dtr(4)=SQRT3*casa*cb
      dtr(5)=SQRT3*casa*sb
      dtr(6)=casa*cb2sb2
      dtr(7)=-SQRT3*casa
      dtr(8)=2.0*casa*cbsb
      dtr(9)=cb*(ca2-sa2)
      dtr(10)=sb*(ca2-sa2)
      dtr(11)=-2.0*sa*cbsb
      dtr(12)=0.0
      dtr(13)=sa*cb2sb2
      dtr(14)=-ptr(5)
      dtr(15)=ptr(4)

      if(period(atom_row)%orb(3)%orb_of_e*period(atom_col)%orb(3)%orb_of_e>0)then 

        dtr(16)=0.5*(1.0+ca2)*cb2sb2
        dtr(17)=0.5*SQRT3*sa2
        dtr(18)=cbsb*(1.0+ca2)
        dtr(19)=-casa*cb
        dtr(20)=-casa*sb
        dtr(21)=-2.0*ca*cbsb
        dtr(22)=0.0
        dtr(23)=ca*cb2sb2
        dtr(24)=ptr(2)
        dtr(25)=-ptr(1)

      endif

    endif

    if(period(atom_row)%orb(4)%orb_of_e+period(atom_col)%orb(4)%orb_of_e>0)then 

      s2b=2.0*sb*cb
      sa3=sa2*sa
      c2b=(cb2-sb2)
      s3b=(c2b*sb + s2b*cb)
      c3b=(c2b*cb - s2b*sb)
      s2a=2.0*sa*ca
      ftr(1)=0.5*ca*(5.0*ca2 -3.0)
      ftr(2)=SQRT6*0.25*cb*sa*(5.0*ca2 -1.0)
      ftr(3)=SQRT6*0.25*sb*sa*(5.0*ca2 -1.0)
      ftr(4)=SQRT15*0.5*s2b*ca*sa2
      ftr(5)=SQRT15*0.5*c2b*ca*sa2
      ftr(6)=SQRT10*0.25*c3b*sa3
      ftr(7)=SQRT10*0.25*s3b*sa3
      ftr(8)=-SQRT6*0.25*sa*(5.0*ca2 -1.0)
      ftr(9)=0.25*cb*ca*(15.0*ca2 -11.0)
      ftr(10)=0.25*sb*ca*(15.0*ca2 -11.0)
      ftr(11)=SQRT10*0.25*s2b*sa*(3.0*ca2 - 1.0)
      ftr(12)=SQRT10*0.25*c2b*sa*(3.0*ca2 - 1.0)
      ftr(13)=SQRT15*0.25*c3b*ca*sa2
      ftr(14)=SQRT15*0.25*s3b*ca*sa2
      ftr(15)=0.0
      ftr(16)=-0.25*sb*(5.0*ca2 -1.0)
      ftr(17)=0.25*cb*(5.0*ca2 -1.0)
      ftr(18)=SQRT10*0.25*c2b*s2a
      ftr(19)=-SQRT10*s2b*s2a*0.25
      ftr(20)=-SQRT15*0.25*s3b*sa2
      ftr(21)=SQRT15*0.25*c3b*sa2

      if(period(atom_row)%orb(3)%orb_of_e*period(atom_col)%orb(3)%orb_of_e>0)then

        c2a=ca2-sa2
        ftr(22)=0.0
        ftr(23)=SQRT10*0.5*sb*ca*sa
        ftr(24)=-SQRT10*0.5*cb*ca*sa
        ftr(25)=c2b*c2a
        ftr(26)=-s2b*c2a
        ftr(27)=-SQRT6*0.25*s3b*s2a
        ftr(28)=SQRT6*0.25*c3b*s2a
        ftr(29)=SQRT15*0.5*ca*sa2
        ftr(30)=SQRT10*0.25*cb*sa*(1.0 -3.0*ca2)
        ftr(31)=SQRT10*0.25*sb*sa*(1.0 -3.0*ca2)
        ftr(32)=0.5*s2b*ca*(3.0*ca2 -1.0)
        ftr(33)=0.5*c2b*ca*(3.0*ca2 -1.0)
        ftr(34)=SQRT6*0.25*c3b*sa*(1.0 + ca2)
        ftr(35)=SQRT6*0.25*s3b*sa*(1.0 + ca2)

      endif

      if(period(atom_row)%orb(4)%orb_of_e+period(atom_col)%orb(4)%orb_of_e>0)then 

        ftr(36)=-SQRT10*0.25*sa3
        ftr(37)=SQRT15*0.25*cb*ca*sa2
        ftr(38)=SQRT15*0.25*sb*ca*sa2
        ftr(39)=-SQRT6*0.25*s2b*sa*(1.0 + ca2)
        ftr(40)=-SQRT6*0.25*c2b*sa*(1.0 + ca2)
        ftr(41)=0.25*c3b*ca*(3.0 + ca2)
        ftr(42)=0.25*s3b*ca*(3.0 + ca2)
        ftr(43)=0.0
        ftr(44)=-SQRT15*0.25*sb*sa2
        ftr(45)=SQRT15*0.25*cb*sa2
        ftr(46)=-SQRT6*0.25*c2b*s2a
        ftr(47)=SQRT6*0.25*s2b*s2a
        ftr(48)=-0.25*s3b*(1.0 +3.0*ca2)
        ftr(49)=0.25*c3b*(1.0 + 3.0*ca2)

      endif

    endif

    r=AUI*r

    nsrow=period(atom_row)%orb(1)%orb_of_e
    nprow=period(atom_row)%orb(2)%orb_of_e;
    ndrow=period(atom_row)%orb(3)%orb_of_e;
    nfrow=period(atom_row)%orb(4)%orb_of_e;
    nscol=period(atom_col)%orb(1)%orb_of_e;
    npcol=period(atom_col)%orb(2)%orb_of_e;
    ndcol=period(atom_col)%orb(3)%orb_of_e;
    nfcol=period(atom_col)%orb(4)%orb_of_e;

    ! /*      (s_row:s_col)   */
    if(nsrow*nscol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nscol,nsrow,1,1);
      S(1,1)=sigma
    endif

    ! /*      (s_row:p_col)   */
    if(nsrow*npcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,npcol,nsrow,2,1);
      sigma=-sigma;
      do k=2,4
        S(1,k)=ptr(k-1)*sigma
      enddo        
      !         write(*,*)'sp',S(1,2),S(1,3),S(1,4)
    endif

    ! /*      (s_row:d_col)   */
    if(nsrow*ndcol>0)then

      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,ndcol,nsrow,3,1);
      do k=1,6 
        S(1,4+k)=dtr(k)*sigma
      enddo
      !       write(*,*)'sd',S(1,5),S(1,6),S(1,7),S(1,8),S(1,9),S(1,10)
    endif

    ! /*      (s_row:f_col)   */
    if(nsrow*nfcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nfcol,nsrow,4,1);
      sigma=-sigma;
      do k=1,8
        S(1,9+k)=ftr(k)*sigma
      enddo
    endif

    ! /*      (p_row:s_col)   */
    if(nprow*nscol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nscol,nprow,1,2);
      do j=2,4 
        S(j,1)=ptr(j-1)*sigma
      enddo
      !       write(*,*)'ps',S(2,1),S(3,1),S(4,1)
    endif

    ! /*      (p_row:p_col)   */
    if(nprow*npcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,npcol,nprow,2,2);
      sigma=-sigma
      do j=2,4
        do k=2,4
          S(j,k)=ptr(j-1)*ptr(k-1)*sigma+(ptr(j+3-1)*ptr(k+3-1)+ptr(j+6-1)*ptr(k+6-1))*pi
          S(k,j)=S(j,k)
          !           write(*,*)'pp',S(j,k)
        enddo
      enddo

    endif

    ! /*      (p_row:d_col)   */
    if(nprow*ndcol>0)then

      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,ndcol,nprow,3,2);
      pi=-pi
      do j=2,4
        do k=1,6
          S(j,4+k)=ptr(j-1)*dtr(k)*sigma+(ptr(j+3-1)*dtr(k+5)+ptr(j+6-1)*dtr(k+10))*pi;
          !           write(*,*)'pd',j,4+k,S(j,4+k)
        enddo
      enddo

    endif

    ! /*      (p_row:f_col)   */
    if(nprow*nfcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nfcol,nprow,4,2);
      sigma=-sigma
      do j=2,4
        do k=1,8
          S(j,9+k)=ptr(j-1)*ftr(k)*sigma+(ptr(j+3-1)*ftr(k+7)+ptr(j+6-1)*ftr(k+14))*pi;
        enddo
      enddo
    endif

    ! /*      (d_row:s_col)   */
    if(ndrow*nscol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nscol,ndrow,1,3);
      do j=1,6
        S(4+j,1)=dtr(j)*sigma;
      enddo
      !       write(*,*)'ds',S(4+j,1)
    endif        

    ! /*      (d_row:p_col)   */
    if(ndrow*npcol>0)then

      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,npcol,ndrow,2,3);
      sigma=-sigma
      do j=1,6
        do k=2,4
          S(4+j,k)=dtr(j)*ptr(k-1)*sigma+(dtr(j+5)*ptr(k+3-1)+dtr(j+10)*ptr(k+6-1))*pi
        enddo
      enddo
      !               write(*,*)'dp',S(4+j,k)
    endif

    ! /*      (d_row:d_col)   */
    if(ndrow*ndcol>0)then

      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,ndcol,ndrow,3,3);
      pi=-pi
      do j=1,5
        do k=1,5
          S(4+j,4+k)=dtr(j)*dtr(k)*sigma+(dtr(j+5)*dtr(k+5)+dtr(j+10)*dtr(k+10))*pi+(dtr(j+15)*dtr(k+15)+dtr(j+20)*dtr(k+20))*delta
          S(4+k,4+j)=S(4+j,4+k)
        enddo
      enddo
      !       write(*,*)"dd",S(4+j,4+k)
    endif

    ! /*      (d_row:f_col)   */
    if(ndrow*nfcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nfcol,ndrow,4,3);
      sigma=-sigma
      delta=-delta
      do j=1,6
        do k=1,8
          S(4+j,9+k)=dtr(j)*ftr(k)*sigma+(dtr(j+5)*ftr(k+7)+dtr(j+10)*ftr(k+14))*pi+(dtr(j+15)*ftr(k+21)+dtr(j+20)*ftr(k+28))*delta;
        enddo
      enddo
    endif

    ! /*      (f_row:s_col)   */
    if(nfrow*nscol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nscol,nfrow,1,4);
      do j=1,8
        S(9+j,1)=ftr(j)*sigma
      enddo
    endif

    ! /*      (f_row:p_col)   */
    if(nfrow*npcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,npcol,nfrow,2,4);
      sigma=-sigma;
      do j=1,8
        do k=2,4
          S(9+j,k)=ftr(j)*ptr(k-1)*sigma+(ftr(j+7)*ptr(k+3-1)+ftr(j+14)*ptr(k+6-1))*pi;
        enddo
      enddo
    endif

    ! /*      (f_row:d_col)   */
    if(nfrow*ndcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,ndcol,nfrow,3,4);
      pi=-pi
      do j=1,8
        do k=1,6
          S(9+j,4+k)=ftr(j)*dtr(k)*sigma+(dtr(k+5)*ftr(j+7)+dtr(k+10)*ftr(j+14))*pi+(ftr(j+21)*dtr(k+15)+ftr(j+28)*dtr(k+20))*delta
        enddo
      enddo
    endif

    ! /*      (f_row:f_col)   */
    if(nfrow*nfcol>0)then
      call mov(period,sigma,pi,delta,phi,atom_col,atom_row,r,nfcol,nfrow,4,4);
      sigma=-sigma
      delta=-delta
      do j=1,8
        do k=1,8
          S(9+j,9+k)=ftr(j)*ftr(k)*sigma+(ftr(j+7)*ftr(k+7)+ftr(j+14)*ftr(k+14))*pi+(ftr(j+21)*ftr(k+21)+ftr(j+28)*ftr(k+28))*delta+(ftr(j+35)*ftr(k+35)+ftr(j+42)*ftr(k+42))*phi
          S(9+k,9+j)=S(9+j,9+k)
        enddo
      enddo
    endif

    do j=1,16

      if(j==1)                VSIP_row=period(atom_row)%orb(1)%VSIP;
      if(4>=j.and.j>1)       VSIP_row=period(atom_row)%orb(2)%VSIP;
      if(9>=j.and.j>4)       VSIP_row=period(atom_row)%orb(3)%VSIP;
      if(16>=j.and.j>9)      VSIP_row=period(atom_row)%orb(4)%VSIP;

      do k=1,16


        if(k==1)                VSIP_col=period(atom_col)%orb(1)%VSIP;
        if(4>=k.and.k>1)       VSIP_col=period(atom_col)%orb(2)%VSIP;
        if(9>=k.and.k>4)       VSIP_col=period(atom_col)%orb(3)%VSIP;
        if(16>=k.and.k>9)      VSIP_col=period(atom_col)%orb(4)%VSIP;

        H(j,k)=(VSIP_row+VSIP_col)*S(j,k)*0.875_dp
      enddo
    enddo

  end subroutine overlap

  subroutine mov(period,sigma,pi,delta,phi,atom_col,atom_row,rr,n1,n2,l1,l2)

    implicit none 

    real(dp) :: sigma,pi,delta,phi, rr
    integer :: atom_col, atom_row, n1,n2,l1,l2
    type(atom_parameter), intent(in) :: period(:)  
    integer i,nn,ia,ib,lc,ld,ik,il,ij,maxcal;
    real(dp) :: sk1,sk2,rll(4),xx,yy
    real(dp) :: aa(30),bb(30),a(30),b(30)

    !nn=(l1<l2)?l1+1:l2+1;
    if(l1 < l2)then 
      nn = l1+1
    else
      nn = l2+1
    endif

    sigma=0.0
    pi=0.0
    delta=0.0
    phi=0.0

    ia=1
    ib=1
    maxcal=n1+n2


    do i=1,30

      a(i)=0.0
      b(i)=0.0

    enddo
    do i=1,3

      rll(i)=0.0
    enddo

    if (period(atom_col)%orb(l1)%coef2(2).ne.0) ia=2
    if (period(atom_row)%orb(l2)%coef2(2).ne.0) ib=2

    do ik=1,ia

      do il=1,ib

        sk1=period(atom_col)%orb(l1)%exp2(ik)

        sk2=period(atom_row)%orb(l2)%exp2(il)

        call abfns(a,b,sk1,sk2,rr,maxcal)

        do ij=1,nn-1
          !           write(*,*)"nn",nn
          rll(ij)=lovlap(a,b,sk1,sk2,rr,l1,l2,ij,n1,n2)

        enddo

        xx=period(atom_col)%orb(l1)%coef2(ik)
        yy=period(atom_row)%orb(l2)%coef2(il)

        sigma=sigma+xx*yy*rll(1)
        pi=pi+xx*yy*rll(2)
        delta=delta+xx*yy*rll(3)
        phi=phi+xx*yy*rll(4)

      enddo

    enddo

  end subroutine mov 

  subroutine abfns(a,b,sk1,sk2,rr,maxcal)

    real(dp) :: a(:),b(:),sk1,sk2,rr
    integer :: maxcal
    integer :: i,j,il,ix,ir,is,k,in
    real(dp) ::  rho1,rho2,c,d,h,r,ra,rho22,t,tr,temp

    j=maxcal+1

    rho1=.5*(sk1+sk2)*rr
    rho2=.5*(sk1-sk2)*rr
    if(rho1>165.or.rho2>165)then

      do i=1,20

        a(i)=0
        b(i)=0
      enddo
      return
    endif

    c=exp(-rho1)
    a(1)=c/rho1

    do i=2,j
      a(i)=(real(i-1)*a(i-1)+c)/rho1

    enddo

    ix=j
    !ir=(rho2>0)? (2*rho2) : -(2*rho2);    
    if(rho2>0)then
      ir = 2*rho2
    else
      ir = -(2*rho2)
    endif

    !is=(ir+1<19)?(ir+1):19;
    if(ir+1<19)then
      is = ir+1
    else
      is = 19
    endif

    if(rho2==0.0_dp)then
      do i=1,ix,2         
        b(i)=2.0/real(i)
        b(i+1)=0
      enddo
      return
    endif

    d=exp(rho2)
    h=1/d

    r=d-h

    if(r>0)then
      temp = r
    else
      temp = -r
    endif

    !         temp=(r>0)?r:-r;
    if(temp<0.1)then
      ra=rho2
      rho22=rho2*rho2
      t=rho2
      do i=2,50,2
        t=t*rho22/real(i*i+i)
        ra=ra+t
        if(t<1e-30) exit
      enddo
      r=ra+ra
    endif

    b(1)=r/rho2

    do i=2,ix,is

      if(ir.ne.0)then
        il=is-1
        if(1.le.il)then
          do j=1,il
            k=i+j-1
            !(float)k/2 .ne. int(k/2))
            if(mod(k,2).ne.0)then 
              b(k)=(r+real(k-1)*b(k-1))/rho2
            else
              b(k)=-(d+h-real(k-1)*b(k-1))/rho2
            endif
          enddo
        endif
      endif

      in=i+is-1
      if(in-ix > 0)then                    
        return
      endif

      !(float)in/2 ==(int)(in/2)) 
      if(mod(in,2)==0)then

        tr=rho2
        b(in)=-2.0*tr/real(in+1)
        do j=1,500
          tr=tr*rho2*rho2/real((2*j)*(2*j+1))
          temp=tr/b(in)
          !temp=(temp>0)?temp:-temp;          
          if(temp>0)then
            temp = temp
          else
            temp = -temp
          endif

          if(temp-1.0e-7 <= 0)exit

          b(in)=b(in)-2.0*tr/real(in+1+2*j)
        enddo

      else

        tr=1.0

        b(in)=2.0*tr/real(in)

        do j=1,500

          tr=tr*rho2*rho2/real((2*j)*(2*j-1))
          temp=tr/b(in)
          if(temp>0)then
            temp = temp
          else
            temp = -temp
          endif

          !temp=(temp>0)?temp:-temp;
          if(temp-1.0e-7 <= 0)exit

          b(in)=b(in)+2.0*tr/real(in+2*j)
        enddo
      endif
    enddo

  end subroutine abfns


  real(dp) function lovlap(a,b,sk1,sk2,r,ll1,ll2,mm1,n1,n2)

    real(dp) :: a(30),b(30),sk1,sk2,r 
    integer :: l1,l2,ll1,ll2,m1,mm1,n1,n2
    real(dp) fact(25)
    integer m2,jend,kend,ieb,i,j,k,i1,i2,i3,i4,i5,i6,iev,iab,ibb,icb,idb,ju,ku,ip,ir,ex
    real(dp) strad,value,rhoa,rhob,rhoap,rhoab,rhopo,terma,term,value1,value2,value3,value4,con1,con12,te1
    integer bincoe(8,8)

    !         [8][8]={1,1,1,1,1,1,1,1, 0,1,2,3,4,5,6,7, 0,0,1,3,6,10,15,21, 0,0,0,1,4,10,20,35, 0,0,0,0,1,5,15,35, 0,0,0,0,0,1,6,21, 0,0,0,0,0,0,1,7, 0,0,0,0,0,0,0,1};

    bincoe(1,:) = (/ 1,1,1,1,1,1,1,1/)
    bincoe(2,:) = (/0,1,2,3,4,5,6,7/)
    bincoe(3,:) = (/0,0,1,3,6,10,15,21/)
    bincoe(4,:) = (/0,0,0,1,4,10,20,35/)
    bincoe(5,:) = (/0,0,0,0,1,5,15,35/)
    bincoe(6,:) = (/0,0,0,0,0,1,6,21/)
    bincoe(7,:) = (/0,0,0,0,0,0,1,7/)
    bincoe(8,:) = (/0,0,0,0,0,0,0,1/)

    strad=0.0

    fact(1)=1  !fact(0)=1

    do i=2,25        
      fact(i)=real(i-1)*fact(i-1) !fact is shifted
    enddo

    m1=mm1-1
    m2=m1
    rhoa=r*sk1
    rhob=r*sk2
    rhoap=rhoa**(2*n1+1)
    rhoab=rhob**(2*n2+1)
    rhopo=rhoap*rhoab

    l1=ll1-1
    l2=ll2-1

    terma=0.5**(l1+l2+1)*sqrt(real((l1+l1+1)*(l2+l2+1))*fact(l1-m1+1)&
      *fact(l2-m1+1)/(fact(n1+n1+1)*fact(n2+n2+1)*fact(l1+m1+1)*fact(l2+m1+1))*rhopo);  

    jend=1+int((l1-m1)/2)
    kend=1+int((l2-m2)/2)
    ieb=m1+1

    do j=1,jend

      ju=j-1
      iab=n1-l1+ju+ju+1
      icb=l1-m1-ju-ju+1
      con1=fact(l1+l1-ju-ju+1)/Real(fact(l1-m1-ju-ju+1)*fact(ju+1)*fact(l1-ju+1))
      do k=1,kend

        ku=k-1
        con12=con1*fact(l2+l2-ku-ku+1)/Real(fact(l2-m2-ku-ku+1)*fact(ku+1)*fact(l2-ku+1))
        iev=ju+ku+l2
        !if(int(iev/2)!=(float)iev/2) con12=-con12;
        if(mod(iev,2).ne.0) con12=-con12
        ibb=n2-l2+ku+ku+1
        idb=l2-m2-ku-ku+1
        value=0.0
        do i6=1,ieb

          do i5=1,ieb

            value1=bincoe(i6,ieb)*bincoe(i5,ieb)
            iev=i5+i6
            !if((int)(iev/2)!=(float)iev/2) value1=-value1;
            if(mod(iev,2).ne.0) value1=-value1
            do i4=1,idb

              value1=-value1
              value2=bincoe(i4,idb)*value1

              do i3=1,icb

                value3=bincoe(i3,icb)*value2

                do i2=1,ibb

                  value3=-value3
                  value4=bincoe(i2,ibb)*value3

                  do i1=1,iab

                    term=value4*bincoe(i1,iab)
                    ir=i1+i2+ieb+ieb-i6-i6-i3+idb-i4+icb-1
                    ip=iab-i1+ibb-i2+2*ieb-2*i5+icb-i3+idb-i4+1
                    value=value+a(ip)*b(ir)*term
                  enddo
                enddo
              enddo
            enddo
          enddo

        enddo
        strad=strad+value*con12
      enddo
    enddo               

    strad=strad*terma
    lovlap = strad

  end function lovlap

end module huckel_latte_mod
