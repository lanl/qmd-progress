!> A module to rotate the coordinates of a sybsystem in chemical systems.
!! \brief It works by specifying two orientations and a rotation point.
!! @ingroup PROGRESS
!!
module prg_syrotation_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> Rotation type
  type, public :: rotation_type
     character(20) :: jobname
     character(50) :: typeofrot
     !> Atomic point to determine the initial orientation
     integer :: patom1
     !> Atomic point to determine initial orientation
     integer :: patom2
     !> Atomic point to determine the rotation center
     integer :: catom
     !> Atomic point to determine a second rotation center
     integer :: catom2
     !> Point to determine initial orientation
     real(dp) :: pq1(3)
     !> Point to determine final orientation
     real(dp) :: pq2(3)
     !> Initial orientation
     real(dp) :: v1(3)
     !>Final orientation
     real(dp) :: v2(3)
     !> Center of rotation
     real(dp) :: vQ(3)
     !> First and last rotated atom in the list
     integer :: rotate_atoms(2)
  end type rotation_type

  public :: prg_parse_rotation, prg_rotate

contains

  !> The parser for rotation
  !!
  subroutine prg_parse_rotation(rot, filename)

    use prg_kernelparser_mod
    implicit none
    integer, parameter                ::  nkey_char = 2, nkey_int = 6, nkey_re = 15, nkey_log = 1
    character(20)                     ::  jobname
    character(len=*)                  ::  filename
    type(rotation_type), intent(inout)   ::  rot

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'Jobname=','TypeOfRotation=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob','DirToDir']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Atom1=','Atom2=','CAtom=','CAtom2=','ListToRot1=','ListToRot2=']
    integer :: valvector_int(nkey_int) = (/ &
         0,0,0,0,1,10000000/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'PQ1X=','PQ1Y=','PQ1Z=','PQ2X=','PQ2Y=','PQ2Z=','V1X=','V1Y=', &
         &'V1Z=','V2X=','V2Y=','V2Z=','VQX=','VQY=','VQZ=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'DUMMY=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'ROTATION{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    rot%jobname= valvector_char(1)
    rot%typeofrot= valvector_char(2)

    rot%patom1= valvector_int(1)
    rot%patom2= valvector_int(2)
    rot%catom= valvector_int(3)
    rot%catom2= valvector_int(4)
    rot%rotate_atoms(1)= valvector_int(5)
    rot%rotate_atoms(2)= valvector_int(6)

    rot%pq1(1)= valvector_re(1)
    rot%pq1(2)= valvector_re(2)
    rot%pq1(3)= valvector_re(3)

    rot%pq2(1)= valvector_re(4)
    rot%pq2(2)= valvector_re(5)
    rot%pq2(3)= valvector_re(6)

    rot%v1(1)= valvector_re(7)
    rot%v1(2)= valvector_re(8)
    rot%v1(3)= valvector_re(9)

    rot%v2(1)= valvector_re(10)
    rot%v2(2)= valvector_re(11)
    rot%v2(3)= valvector_re(12)

    rot%vQ(1)= valvector_re(13)
    rot%vQ(2)= valvector_re(14)
    rot%vQ(3)= valvector_re(15)

  end subroutine prg_parse_rotation


  !> Rotation routine.
  !! \brief It works by indicating the orientations (v1 and v1) and a rotation center.
  !! The orientation can be passed either directly by setting v1 and v2 or by
  !! indicating two points pQ1 and pQ2. Orientation can also be specified with an
  !! atom position if patom1 and patom2 indices are not zero this atoms are used
  !! to determine the initial and final orientation.
  !! \param rot Rotation type
  !! \param r Coordinates to be rotated
  !! \param verbose Verbosity level
  !!
  !! Example:
  !! \verbatim
  !!        rot%patom1 = 4
  !!        rot%patom2 = 0
  !!        rot%catom2 = 6
  !!        rot%v2 = 0.0 ; rot%v2(1) = 1
  !!        call prg_rotate(rot,r)
  !! \endverbatim
  !! The latter will orient the system such that atom 4 points to the (1,0,0)
  !! direction.
  !!
  subroutine prg_rotate(rot,r,verbose)
    integer                          ::  Natoms, catom, catom2, i
    integer                          ::  k, l, patom1, patom2
    integer                          ::  rotate_atoms(2)
    integer, intent(in)              ::  verbose
    real(dp)                         ::  Mat(3,3), P1MinusC1(3), Vtr(3)
    real(dp)                         ::  alpha, angle, d, deformation
    real(dp)                         ::  dotprod, dv1, dv2, hinge(3)
    real(dp)                         ::  hinge2, hingeDotP1MinusC1, pQ1(3), pQ2(3)
    real(dp)                         ::  pi, v1(3), v11(3), v1mod
    real(dp)                         ::  v1x, v1y, v1z, v1zversor
    real(dp)                         ::  v2(3), v22(3), vN(3), vQ(3)
    real(dp)                         ::  vad(3), vdirn(3), vp(3), vr11(3)
    real(dp)                         ::  vr2(3), vr22(3)
    real(dp), allocatable            ::  rr(:,:), v(:,:), vr(:,:), vr1(:,:)
    real(dp), intent(inout)          ::  r(:,:)
    type(rotation_type), intent(in)  ::  rot

    pi=3.1415926535897993d0
    alpha= 0.0_dp
    alpha= (alpha/180.0_dp)*Pi

    natoms = size(r, dim=2)

    allocate(v(3,natoms))
    allocate(vr(3,natoms))
    allocate(vr1(3,natoms))
    allocate(rr(3,natoms))

    vr=0.0_dp
    vr1=0.0_dp
    rr=0.0_dp
    v=0.0_dp

    patom1=rot%patom1 		!Rotation atom indices
    patom2=rot%patom2
    catom=rot%catom       !Rotation center
    catom2=rot%catom2			!Second rotation center

    pq1(1)=rot%pq1(1)     !Initial director point
    pq1(2)=rot%pq1(2)
    pq1(3)=rot%pq1(3)

    pq2(1)=rot%pq2(1)     !Final director point
    pq2(2)=rot%pq2(2)
    pq2(3)=rot%pq2(3)

    v1(1)=rot%v1(1) 	    !Initial anchored vector
    v1(2)=rot%v1(2)
    v1(3)=rot%v1(3)

    v2(1)=rot%v2(1) 	    !Final anchored vector
    v2(2)=rot%v2(2)
    v2(3)=rot%v2(3)

    vQ(1)=rot%vQ(1)	      !Rotation center
    vQ(2)=rot%vQ(2)
    vQ(3)=rot%vQ(3)

    rotate_atoms(1)=rot%rotate_atoms(1)   !index for first rotated atom
    rotate_atoms(2)=rot%rotate_atoms(2)

    if(rotate_atoms(2).gt.natoms) rotate_atoms(2)=natoms

    if(patom1.ne.0)then
       pq1(1)=r(1,patom1)
       pq1(2)=r(2,patom1)
       pq1(3)=r(3,patom1)

       pq2(1)=r(1,patom2)
       pq2(2)=r(2,patom2)
       pq2(3)=r(3,patom2)

       vQ(1)=r(1,catom)
       vQ(2)=r(2,catom)
       vQ(3)=r(3,catom)
    endif

    if(catom2.ne.0)then
       hinge=r(:,catom2)-r(:,catom)
       hinge2=hinge(1)**2+hinge(2)**2+hinge(3)**2
       p1MinusC1=pq1(:)-r(:,catom)
       hingeDotp1MinusC1=hinge(1)*p1MinusC1(1) + hinge(2)*p1MinusC1(2) + hinge(3)*p1MinusC1(3)
       VQ(:)=r(:,catom) + hinge(:)*hingeDotp1MinusC1/hinge2
       write(*,*)VQ(:)
       write(*,*)hinge(:)
       write(*,*)hingeDotp1MinusC1/hinge2
    endif

    if(V1(1).eq.0.0_dp.and.V1(2).eq.0.0_dp.and.V1(3).eq.0.0_dp)then
       v1=pq1-vQ
    endif

    if(V2(1).eq.0.0_dp.and.V2(2).eq.0.0_dp.and.V2(3).eq.0.0_dp)then
       v2=pq2-vQ
    endif

    vtr(1)=0.0_dp	    !Translation
    vtr(2)=0.0_dp
    vtr(3)=0.0_dp

    dv1=dsqrt(v1(1)**2+v1(2)**2+v1(3)**2)  !Angle between v1 y v2
    dv2=dsqrt(v2(1)**2+v2(2)**2+v2(3)**2)
    dotprod=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    alpha=-acos(dotprod/(dv1*dv2))

    if(alpha.eq.0.0_dp)write(*,*)'V1 and V2 are in the same direction'
    if(abs(alpha-Pi).lt.1.0E-10)write(*,*)'V1 and V2 are in the same direction'

    if(verbose >= 1) write(*,*)'alpha=',alpha

    do i=rotate_atoms(1),rotate_atoms(2)
       r(1,i) = r(1,i)-VQ(1)
       r(2,i) = r(2,i)-VQ(2)
       r(3,i) = r(3,i)-VQ(3)
    enddo

    vad=v1+v2

    vN(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vN(2)=-v1(1)*v2(3)+v2(1)*v1(3)
    vN(3)=v1(1)*v2(2)-v2(1)*v1(2)

    vp(1)=vn(2)*vad(3)-vn(3)*vad(2)
    vp(2)=-vn(1)*vad(3)+vad(1)*vn(3)
    vp(3)=vn(1)*vad(2)-vad(1)*vn(2)

    d=dsqrt(vad(1)**2+vad(2)**2+vad(3)**2)
    vad=vad/d

    d=dsqrt(vp(1)**2+vp(2)**2+vp(3)**2)
    vp=vp/d

    d=dsqrt(vn(1)**2+vn(2)**2+vn(3)**2)
    vn=vn/d

    !Projection onto vad vp and vn

    do i=rotate_atoms(1),rotate_atoms(2)
       v(1,i)=r(1,i)*vad(1)+r(2,i)*vad(2)+r(3,i)*vad(3)
       v(2,i)=r(1,i)*vp(1)+r(2,i)*vp(2)+r(3,i)*vp(3)
       v(3,i)=r(1,i)*vn(1)+r(2,i)*vn(2)+r(3,i)*vn(3)
    enddo

    !Rotation matrix

    Mat(1,1)=cos(alpha)
    Mat(2,1)=sin(alpha)
    Mat(3,1)=0.0_dp
    Mat(1,2)=-sin(alpha)
    Mat(2,2)=cos(alpha)
    Mat(3,2)=0.0_dp
    Mat(1,3)=0.0_dp
    Mat(2,3)=0.0_dp
    Mat(3,3)=1.0_dp

    !Rotation

    do l=rotate_atoms(1),rotate_atoms(2)
       do i=1,3
          do k=1,3
             vr(i,l)=vr(i,l)+Mat(k,i)*v(k,l)
          enddo
       enddo
    enddo

    !Basis set change

    do i=rotate_atoms(1),rotate_atoms(2)
       vr1(1,i)=vr(1,i)*vad(1) + vr(2,i)*vp(1) + vr(3,i)*vN(1)
       vr1(2,i)=vr(1,i)*vad(2) + vr(2,i)*vp(2) + vr(3,i)*vN(2)
       vr1(3,i)=vr(1,i)*vad(3) + vr(2,i)*vp(3) + vr(3,i)*vN(3)
    enddo

    if(alpha.eq.0.0_dp)Vr1=r
    if(abs(alpha-Pi).lt.1.0E-10)Vr1=-r

    do i=rotate_atoms(1),rotate_atoms(2)
       rr(1,i)=vr1(1,i)+vQ(1)+vtr(1)
       rr(2,i)=vr1(2,i)+vQ(2)+vtr(2)
       rr(3,i)=vr1(3,i)+vQ(3)+vtr(3)
    enddo

    do i=rotate_atoms(1),rotate_atoms(2)
       r(1,i)=rr(1,i)
       r(2,i)=rr(2,i)
       r(3,i)=rr(3,i)
    enddo

    deformation=0.0d0

    do i=rotate_atoms(1),rotate_atoms(2)
       d=dsqrt((r(1,1)-r(1,i))**2+(r(2,1)-r(2,i))**2+(r(3,1)-r(3,i))**2)
       deformation=deformation + d
       d=dsqrt((rr(1,1)-rr(1,i))**2+(rr(2,1)-rr(2,i))**2+(rr(3,1)-rr(3,i))**2)
       deformation=deformation - d
    enddo

    if(deformation >= 0.001_dp)then
       write(*,*)'Rotation failed ...'
       write(*,*)'Deformation=',deformation
    endif

  end subroutine prg_rotate

end module prg_syrotation_mod
