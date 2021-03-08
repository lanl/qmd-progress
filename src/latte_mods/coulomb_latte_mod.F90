!> Module to compute the coulombic contribution of the Forces.
!! \ingroup LATTE
!!
module coulomb_latte_mod

  use bml
  use prg_extras_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: get_ewald_real, get_ewald_recip, get_coulcut
  public :: get_ewald_list_real, get_ewald_list_real_dcalc

contains

  !> This routine prg_initializes the Coulombic rcut parameter.
  !! \param coul_acc Coulomb accuracy.
  !! \param timeratio Estimated ration between real and k-space time efficiency.
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param coulcut Coulombic cutoff.
  subroutine get_coulcut(coulomb_acc,timeratio,nats,lattice_vectors,coulcut)
    integer, intent(in)    ::  nats
    real(dp)               ::  Lx, Ly, Lz, calpha
    real(dp)               ::  coulvol, pi, sqrtx
    real(dp), intent(in)   ::  coulomb_acc, lattice_vectors(:,:), timeratio
    real(dp), intent(out)  ::  coulcut

    pi = 3.14159265358979323846264338327950_dp

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)

    coulvol = lx*ly*lz;
    sqrtx = sqrt(-log(coulomb_acc));

    calpha = sqrt(pi)*((timeratio*nats/(coulvol**2.0_dp))**(1.0_dp/6.0_dp))
    coulcut = sqrtx/calpha;
    if (coulcut > 50) then
       coulcut = 50
    endif

  end subroutine get_coulcut

  !> This routine computes the real space contribution of the Ewald summation.
  !! \param spindex Species index list.
  !! \param splist Element symbol for every species.
  !! \param coordinates Coordinates for every atom in the system.
  !! \param charges Charges for every atom in the system.
  !! \param hubbardu Hubbard parameter U for every species.
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param volr Volume of the system cell.
  !! \param coul_acc Coulomb accuracy.
  !! \param coul_forces_r Coulombic forces (real space contribution)
  !! \param coul_pot_r Coulombic potential (real space contribution)
  subroutine get_ewald_real(spindex,splist,coordinates,charges,hubbardu&
       ,lattice_vectors,volr,coul_acc,coul_forces_r,coul_pot_r)

    character(2), intent(in)             ::  splist(:)
    integer                              ::  atomi, i, j, nats
    integer                              ::  nr_shift_x, nr_shift_y, nr_shift_z
    integer, intent(in)                  ::  spindex(:)
    real(dp)                             ::  a2xa3(3), ca, calpha, calpha2
    real(dp)                             ::  coul_acc, coulcut, coulcut2, coulombv
    real(dp)                             ::  coulvol, dc(3), dr, expti
    real(dp)                             ::  exptj, fcoul(3), force, keconst
    real(dp)                             ::  magr, magr2, numrep_erfc, pi
    real(dp)                             ::  r0b(3), ra(3), rab(3), rb(3)
    real(dp)                             ::  relperm, sa, sb, sc
    real(dp)                             ::  sd, se, sf, sqrtp
    real(dp)                             ::  sqrtpi, sqrtx, ssa, ssb
    real(dp)                             ::  ssc, ssd, sse, tfact
    real(dp)                             ::  ti, ti2, ti2mtj2, ti3
    real(dp)                             ::  ti4, ti6, timeratio, tj
    real(dp)                             ::  tj2, tj2mti2, tj3, tj4
    real(dp)                             ::  tj6, z
    real(dp), allocatable, intent(inout)  ::  coul_forces_r(:,:), coul_pot_r(:)
    real(dp), intent(in)                 ::  charges(:), coordinates(:,:), hubbardu(:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  volr

    !Estimated ratio between real & k space
    timeratio = 10.0_dp
    pi = 3.14159265358979323846264338327950_dp

    nats = size(charges,dim=1)

    if(.not.allocated(coul_forces_r))allocate(coul_forces_r(3,nats))
    if(.not.allocated(coul_pot_r))allocate(coul_pot_r(nats))

    coul_pot_r = 0.0_dp
    coul_forces_r = 0.0_dp

    sqrtx = sqrt(-log(coul_acc));
    calpha = sqrt(pi)*((timeratio*nats/(volr**2))**(1.0_dp/6.0_dp));
    coulcut = sqrtx/calpha;
    calpha2 = calpha*calpha;
    if (coulcut > 50.0_dp) then
       coulcut = 50.0_dp;
       calpha = sqrtx/coulcut;
    endif
    coulcut2 = coulcut*coulcut
    calpha2 = calpha*calpha

    relperm = 1.0_dp;

    !! The 14.399 factor corresponds to 1/(4*pi*epsilon0) in eV*Ang
    keconst = 14.3996437701414_dp*relperm;
    tfact  = 16.0_dp/(5.0_dp*keconst);

    sqrtpi = sqrt(pi);

    nr_shift_x = 0
    nr_shift_y = 0
    nr_shift_z = 0

    do i =1,nats

       fcoul = 0.0_dp
       coulombv = 0.0_dp

       ti = tfact*hubbardu(spindex(i));

       ti2 = ti*ti;
       ti3 = ti2*ti;
       ti4 = ti2*ti2;
       ti6 = ti4*ti2;

       ssa = ti;
       ssb = ti3/48.0_dp;
       ssc = 3.0_dp*ti2/16.0_dp;
       ssd = 11.0_dp*ti/16.0_dp;
       sse = 1.0_dp;

       ra = coordinates(:,i);

       do j = 1,nats

          r0b = coordinates(:,j);

          do nr_shift_x = -1,1  ! periodic bc shifts in x, y and z. costs a lot extra!
             do nr_shift_y = -1,1
                do nr_shift_z = -1,1
                   !rb = r0b + nr_x*v1 + nr_y*v2 + nr_z*v3
                   rb(1) = r0b(1) + nr_shift_x*lattice_vectors(1,1) ! shifts for pbc
                   rb(1) = rb(1) + nr_shift_y*lattice_vectors(2,1) ! shifts for pbc
                   rb(1) = rb(1) + nr_shift_z*lattice_vectors(3,1) ! shifts for pbc

                   rb(2) = r0b(2) + nr_shift_y*lattice_vectors(2,2) ! shifts for pbc
                   rb(2) = rb(2) + nr_shift_x*lattice_vectors(1,2) ! shifts for pbc
                   rb(2) = rb(2) + nr_shift_z*lattice_vectors(3,2) ! shifts for pbc

                   rb(3) = r0b(3) + nr_shift_z*lattice_vectors(3,3) ! shifts for pbc
                   rb(3) = rb(3) + nr_shift_y*lattice_vectors(2,3) ! shifts for pbc
                   rb(3) = rb(3) + nr_shift_x*lattice_vectors(1,3) ! shifts for pbc

                   rab = rb-ra  ! obs b - a !!!
                   dr = sqrt(rab(1)**2+rab(2)**2+rab(3)**2)
                   magr = dr
                   magr2 = dr*dr

                   if (dr <= coulcut .and. dr > 1e-12) then
                      tj = tfact*hubbardu(spindex(j))
                      dc = rab/dr
                      z = abs(calpha*magr)
                      numrep_erfc = erfc(z)
                      ca = numrep_erfc/magr
                      coulombv = coulombv + charges(j)*ca
                      ca = ca + 2.0_dp*calpha*exp( -calpha2*magr2 )/sqrtpi
                      force = -keconst*charges(i)*charges(j)*ca/magr
                      expti = exp(-ti*magr)

                      if (splist(spindex(i)) == splist(spindex(j)))then
                         coulombv = coulombv - charges(j)*expti*(ssb*magr2 + ssc*magr + ssd + sse/magr)
                         force = force + (keconst*charges(i)*charges(j)*expti)*((sse/magr2 - 2*ssb*magr - ssc) +&
                              ssa*(ssb*magr2 + ssc*magr + ssd + sse/magr))
                      else
                         tj2 = tj*tj
                         tj3 = tj2*tj
                         tj4 = tj2*tj2
                         tj6 = tj4*tj2
                         exptj = exp( -tj*magr )
                         ti2mtj2 = ti2 - tj2
                         tj2mti2 = -ti2mtj2
                         sa = ti
                         sb = tj4*ti/(2 * ti2mtj2 * ti2mtj2)
                         sc = (tj6 - 3*tj4*ti2)/(ti2mtj2 * ti2mtj2 * ti2mtj2)
                         sd = tj
                         se = ti4*tj/(2 * tj2mti2 * tj2mti2)
                         sf = (ti6 - 3*ti4*tj2)/(tj2mti2 * tj2mti2 * tj2mti2)

                         coulombv = coulombv - (charges(j)*(expti*(sb - (sc/magr)) + exptj*(se - (sf/magr))))
                         force = force + keconst*charges(i)*charges(j)*((expti*(sa*(sb - (sc/magr)) - (sc/magr2))) +&
                              (exptj*(sd*(se - (sf/magr)) - (sf/magr2))))

                      endif

                      fcoul = fcoul + dc*force

                   endif

                enddo
             enddo
          enddo

       enddo
       coul_forces_r(:,i) = fcoul
       coul_pot_r(i) = coulombv
    enddo

    coul_pot_r = keconst*coul_pot_r

  end subroutine get_ewald_real

  !> This routine computes the real space contribution of the Ewald summation using a neighbor list.
  !! \param spindex Species index list.
  !! \param splist Element symbol for every species.
  !! \param coordinates Coordinates for every atom in the system.
  !! \param charges Charges for every atom in the system.
  !! \param hubbardu Hubbard parameter U for every species.
  !! \param volr Volume of the system cell.
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param coul_acc Coulomb accuracy.
  !! \param nnRx See neighlist_type structure.
  !! \param nnRy See neighlist_type structure.
  !! \param nnRz See neighlist_type structure.
  !! \param nrnnlist See neighlist_type structure.
  !! \param nnType See neighlist_type structure.
  !! \param coul_forces_r Coulombic forces (real space contribution)
  !! \param coul_pot_r Coulombic potential (real space contribution)
  subroutine get_ewald_list_real(spindex,splist,coordinates,charges,hubbardu&
       ,lattice_vectors,volr,coul_acc,timeratio,nnRx,nnRy,nnRz,nrnnlist,nnType&
       ,coul_forces_r,coul_pot_r)

    character(2), intent(in)             ::  splist(:)
    integer                              ::  atomi, i, j, nats
    integer                              ::  nr_shift_x, nr_shift_y, nr_shift_z
    integer                              ::  nnI
    integer, intent(in)                  ::  spindex(:)
    real(dp)                             ::  a2xa3(3), ca, calpha, calpha2
    real(dp)                             ::  coul_acc, coulcut, coulcut2, coulombv
    real(dp)                             ::  coulvol, dc(3), dr, expti
    real(dp)                             ::  exptj, fcoul(3), force, keconst
    real(dp)                             ::  magr, magr2, numrep_erfc, pi
    real(dp)                             ::  r0b(3), ra(3), rab(3), rb(3)
    real(dp)                             ::  relperm, sa, sb, sc
    real(dp)                             ::  sd, se, sf, sqrtp
    real(dp)                             ::  sqrtpi, sqrtx, ssa, ssb
    real(dp)                             ::  ssc, ssd, sse, tfact
    real(dp)                             ::  ti, ti2, ti2mtj2, ti3
    real(dp)                             ::  ti4, ti6, tj
    real(dp)                             ::  tj2, tj2mti2, tj3, tj4
    real(dp)                             ::  tj6, z
    real(dp), allocatable, intent(inout)  ::  coul_forces_r(:,:), coul_pot_r(:)
    real(dp), intent(in)                 ::  charges(:), coordinates(:,:), hubbardu(:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  timeratio
    real(dp), intent(in)                 ::  nnRx(:,:), nnRy(:,:), nnRz(:,:)
    integer, allocatable, intent(in)     ::  nrnnlist(:), nnType(:,:)
    real(dp), intent(in)                 ::  volr

    !Estimated ration between real & k space
    pi = 3.14159265358979323846264338327950_dp

    nats = size(charges,dim=1)

    if(.not.allocated(coul_forces_r))allocate(coul_forces_r(3,nats))
    if(.not.allocated(coul_pot_r))allocate(coul_pot_r(nats))

    coul_pot_r = 0.0_dp
    coul_forces_r = 0.0_dp

    sqrtx = sqrt(-log(coul_acc));
    calpha = sqrt(pi)*((timeratio*nats/(volr**2))**(1.0_dp/6.0_dp));
    coulcut = sqrtx/calpha;
    calpha2 = calpha*calpha;
    if (coulcut > 50.0_dp) then
       coulcut = 50.0_dp;
       calpha = sqrtx/coulcut;
    endif
    coulcut2 = coulcut*coulcut
    calpha2 = calpha*calpha

    relperm = 1.0_dp;

    !! The 14.399 factor corresponds to 1/(4*pi*epsilon0) in eV*Ang
    keconst = 14.3996437701414_dp*relperm;
    tfact  = 16.0_dp/(5.0_dp*keconst);

    sqrtpi = sqrt(pi);

    nr_shift_x = 0
    nr_shift_y = 0
    nr_shift_z = 0

    !$omp parallel do default(none) private(i) &
    !$omp private(fcoul,coulombv) &
    !$omp private(ti,ti2,ti3,ti4,ti6,ssa,ssb,ssc,ssd,sse) &
    !$omp private(tj,tj2,tj3,tj4,tj6,ti2mtj2,sa,sb,sc,sd,se,sf) &
    !$omp private(ra,rb,nnI,dr,rab,magr,magr2,j) &
    !$omp private(dc,z,numrep_erfc,ca,force,expti,exptj,tj2mti2) &
    !$omp shared(nats,hubbardu,spindex,coordinates,sqrtpi,keconst ) &
    !$omp shared(nrnnlist,coulcut,nnType,tfact,nnrx,nnry,nnrz,splist) &
    !$omp shared(coul_forces_r, coul_pot_r, calpha, charges, calpha2)
    do i =1,nats

       fcoul = 0.0_dp
       coulombv = 0.0_dp

       ti = tfact*hubbardu(spindex(i));

       ti2 = ti*ti;
       ti3 = ti2*ti;
       ti4 = ti2*ti2;
       ti6 = ti4*ti2;

       ssa = ti;
       ssb = ti3/48.0_dp;
       ssc = 3.0_dp*ti2/16.0_dp;
       ssd = 11.0_dp*ti/16.0_dp;
       sse = 1.0_dp;

       ra = coordinates(:,i);

       do nnI = 1,nrnnlist(i)

          Rb(1) = nnRx(nnI,i);
          Rb(2) = nnRy(nnI,i);
          Rb(3) = nnRz(nnI,i);
          J = nnType(nnI,i);

          rab = rb-ra  ! obs b - a !!!
          dr = sqrt(rab(1)**2+rab(2)**2+rab(3)**2)
          magr = dr
          magr2 = dr*dr

          if (dr <= coulcut .and. dr > 1e-12) then
             tj = tfact*hubbardu(spindex(j))
             dc = rab/dr
             z = abs(calpha*magr)
             numrep_erfc = erfc(z)
             ca = numrep_erfc/magr
             coulombv = coulombv + charges(j)*ca
             ca = ca + 2.0_dp*calpha*exp( -calpha2*magr2 )/sqrtpi
             force = -keconst*charges(i)*charges(j)*ca/magr
             expti = exp(-ti*magr)

             if (splist(spindex(i)) == splist(spindex(j)))then
                coulombv = coulombv - charges(j)*expti*(ssb*magr2 + ssc*magr + ssd + sse/magr)
                force = force + (keconst*charges(i)*charges(j)*expti)*((sse/magr2 - 2*ssb*magr - ssc) +&
                     ssa*(ssb*magr2 + ssc*magr + ssd + sse/magr))
             else
                tj2 = tj*tj
                tj3 = tj2*tj
                tj4 = tj2*tj2
                tj6 = tj4*tj2
                exptj = exp( -tj*magr )
                ti2mtj2 = ti2 - tj2
                tj2mti2 = -ti2mtj2
                sa = ti
                sb = tj4*ti/(2 * ti2mtj2 * ti2mtj2)
                sc = (tj6 - 3*tj4*ti2)/(ti2mtj2 * ti2mtj2 * ti2mtj2)
                sd = tj
                se = ti4*tj/(2 * tj2mti2 * tj2mti2)
                sf = (ti6 - 3*ti4*tj2)/(tj2mti2 * tj2mti2 * tj2mti2)

                coulombv = coulombv - (charges(j)*(expti*(sb - (sc/magr)) + exptj*(se - (sf/magr))))
                force = force + keconst*charges(i)*charges(j)*((expti*(sa*(sb - (sc/magr)) - (sc/magr2))) +&
                     (exptj*(sd*(se - (sf/magr)) - (sf/magr2))))

             endif

             fcoul = fcoul + dc*force

          endif

       enddo
       coul_forces_r(:,i) = fcoul
       coul_pot_r(i) = coulombv
    enddo
    !$omp end parallel do

    coul_pot_r = keconst*coul_pot_r

  end subroutine get_ewald_list_real


  !> This routine computes the real space contribution of the Ewald summation using a neighbor list.
  !! \param spindex Species index list.
  !! \param splist Element symbol for every species.
  !! \param coordinates Coordinates for every atom in the system.
  !! \param charges Charges for every atom in the system.
  !! \param hubbardu Hubbard parameter U for every species.
  !! \param volr Volume of the system cell.
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param coul_acc Coulomb accuracy.
  !! \param nnRx See neighlist_type structure.
  !! \param nnRy See neighlist_type structure.
  !! \param nnRz See neighlist_type structure.
  !! \param nrnnlist See neighlist_type structure.
  !! \param nnType See neighlist_type structure.
  !! \param coul_forces_r Coulombic forces (real space contribution)
  !! \param coul_pot_r Coulombic potential (real space contribution)
  subroutine get_ewald_list_real_dcalc(spindex,splist,coordinates,charges,hubbardu&
       ,lattice_vectors,volr,coul_acc,timeratio,nnIx,nnIy,&
       nnIz,nrnnlist,nnType&
       ,coul_forces_r,coul_pot_r)

    character(2), intent(in)             ::  splist(:)
    integer                              ::  atomi, i, j, nats
    integer                              ::  nnI
    integer, intent(in)                  ::  spindex(:)
    integer(1),allocatable, intent(in)    ::  nnIx(:,:),nnIy(:,:),nnIz(:,:)
    real(dp)                             ::  a2xa3(3), ca, calpha, calpha2
    real(dp)                             ::  coul_acc, coulcut, coulcut2, coulombv
    real(dp)                             ::  coulvol, dc(3), dr, expti, rmod
    real(dp)                             ::  exptj, fcoul(3), force, keconst
    real(dp)                             ::  magr, magr2, numrep_erfc, pi
    real(dp)                             ::  r0b(3), ra(3), rab(3), rb(3)
    real(dp)                             ::  relperm, sa, sb, sc
    real(dp)                             ::  sd, se, sf, sqrtp
    real(dp)                             ::  sqrtpi, sqrtx, ssa, ssb
    real(dp)                             ::  ssc, ssd, sse, tfact
    real(dp)                             ::  ti, ti2, ti2mtj2, ti3
    real(dp)                             ::  ti4, ti6, tj
    real(dp)                             ::  tj2, tj2mti2, tj3, tj4
    real(dp)                             ::  tj6, z, Lx, Ly, Lz
    real(dp), allocatable, intent(inout)  ::  coul_forces_r(:,:), coul_pot_r(:)
    real(dp), intent(in)                 ::  charges(:), coordinates(:,:), hubbardu(:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  timeratio
    integer, allocatable, intent(in)     ::  nrnnlist(:), nnType(:,:)
    real(dp), intent(in)                 ::  volr
    integer, allocatable                 ::  already(:)

    !Estimated ration between real & k space
    pi = 3.14159265358979323846264338327950_dp

    nats = size(charges,dim=1)

    if(.not.allocated(coul_forces_r))allocate(coul_forces_r(3,nats))
    if(.not.allocated(coul_pot_r))allocate(coul_pot_r(nats))

    coul_pot_r = 0.0_dp
    coul_forces_r = 0.0_dp

    sqrtx = sqrt(-log(coul_acc));
    calpha = sqrt(pi)*((timeratio*nats/(volr**2))**(1.0_dp/6.0_dp));
    coulcut = sqrtx/calpha;
    calpha2 = calpha*calpha;
    if (coulcut > 50.0_dp) then
       coulcut = 50.0_dp;
       calpha = sqrtx/coulcut;
    endif
    coulcut2 = coulcut*coulcut
    calpha2 = calpha*calpha

    relperm = 1.0_dp;

    !! The 14.399 factor corresponds to 1/(4*pi*epsilon0) in eV*Ang
    keconst = 14.3996437701414_dp*relperm;
    tfact  = 16.0_dp/(5.0_dp*keconst);

    sqrtpi = sqrt(pi);

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)

    !$omp parallel do default(none) private(i) &
    !$omp private(fcoul,coulombv) &
    !$omp private(ti,ti2,ti3,ti4,ti6,ssa,ssb,ssc,ssd,sse) &
    !$omp private(tj,tj2,tj3,tj4,tj6,ti2mtj2,sa,sb,sc,sd,se,sf) &
    !$omp private(ra,rb,nni,dr,rab,magr,magr2,j) &
    !$omp private(dc,z,numrep_erfc,ca,force,expti,exptj,tj2mti2,rmod) &
    !$omp shared(nats,hubbardu,spindex,coordinates,sqrtpi,keconst,Lx,Ly,Lz ) &
    !$omp shared(nrnnlist,coulcut,nnType,tfact,nnIx,nnIy,nnIz,splist) &
    !$omp shared(coul_forces_r, coul_pot_r, calpha, charges, calpha2)
    do i =1,nats

       fcoul = 0.0_dp
       coulombv = 0.0_dp

       ti = tfact*hubbardu(spindex(i));

       ti2 = ti*ti;
       ti3 = ti2*ti;
       ti4 = ti2*ti2;
       ti6 = ti4*ti2;

       ssa = ti;
       ssb = ti3/48.0_dp;
       ssc = 3.0_dp*ti2/16.0_dp;
       ssd = 11.0_dp*ti/16.0_dp;
       sse = 1.0_dp;

       ra = coordinates(:,i);

       do nni = 1,nrnnlist(i)

          j = nnType(nni,i);

          if(allocated(nnIx))then
             Rb(1) = coordinates(1,j) + nnIx(nni,i)*Lx
             Rb(2) = coordinates(2,j) + nnIy(nni,i)*Ly
             Rb(3) = coordinates(3,j) + nnIz(nni,i)*Lz
             rab = rb-ra
             rmod = prg_norm2(rab)
           
          else

             Rb(1) = coordinates(1,j)
             Rb(2) = coordinates(2,j)
             Rb(3) = coordinates(3,j)

             rab = rb-ra

             rmod = prg_norm2(rab)
             !rmod = sqrt(rab(1)**2 + rab(2)**2 + rab(3)**2)

          !   if(norm2(rab) > coulcut)then
             if(rmod > coulcut)then
                rab(1) = modulo((Rb(1) - Ra(1) + Lx/2.0_dp),Lx) - Lx/2.0_dp
                rab(2) = modulo((Rb(2) - Ra(2) + Ly/2.0_dp),Ly) - Ly/2.0_dp
                rab(3) = modulo((Rb(3) - Ra(3) + Lz/2.0_dp),Lz) - Lz/2.0_dp
             endif
          endif

          !dr = rmod
          dr = norm2(rab)
          !dr = sqrt(rab(1)**2 + rab(2)**2 + rab(3)**2)

          magr = dr
          magr2 = dr*dr

          if (dr <= coulcut .and. dr > 1e-12) then
             tj = tfact*hubbardu(spindex(j))
             dc = rab/dr
             z = abs(calpha*magr)
             numrep_erfc = erfc(z)
             ca = numrep_erfc/magr
             coulombv = coulombv + charges(j)*ca
             ca = ca + 2.0_dp*calpha*exp( -calpha2*magr2 )/sqrtpi
             force = -keconst*charges(i)*charges(j)*ca/magr
             expti = exp(-ti*magr)

             if (splist(spindex(i)) == splist(spindex(j)))then
                coulombv = coulombv - charges(j)*expti*(ssb*magr2 + ssc*magr + ssd + sse/magr)
                force = force + (keconst*charges(i)*charges(j)*expti)*((sse/magr2 - 2*ssb*magr - ssc) +&
                     ssa*(ssb*magr2 + ssc*magr + ssd + sse/magr))
             else
                tj2 = tj*tj
                tj3 = tj2*tj
                tj4 = tj2*tj2
                tj6 = tj4*tj2
                exptj = exp( -tj*magr )
                ti2mtj2 = ti2 - tj2
                tj2mti2 = -ti2mtj2
                sa = ti
                sb = tj4*ti/(2 * ti2mtj2 * ti2mtj2)
                sc = (tj6 - 3*tj4*ti2)/(ti2mtj2 * ti2mtj2 * ti2mtj2)
                sd = tj
                se = ti4*tj/(2 * tj2mti2 * tj2mti2)
                sf = (ti6 - 3*ti4*tj2)/(tj2mti2 * tj2mti2 * tj2mti2)

                coulombv = coulombv - (charges(j)*(expti*(sb - (sc/magr)) + exptj*(se - (sf/magr))))
                force = force + keconst*charges(i)*charges(j)*((expti*(sa*(sb - (sc/magr)) - (sc/magr2))) +&
                     (exptj*(sd*(se - (sf/magr)) - (sf/magr2))))

             endif
       
             fcoul = fcoul + dc*force
             
          endif
!         write(*,*)"i,fcoul",i,fcoul, dr, coulcut
       enddo

     !$omp critical
       coul_forces_r(:,i) = fcoul
       coul_pot_r(i) = coulombv
     !$omp end critical 

    enddo
    !$omp end parallel do

    coul_pot_r = keconst*coul_pot_r

  end subroutine get_ewald_list_real_dcalc

  !> This routine computes the reciprocal space contribution of the Ewald summation.
  !! \param atomi Atom index where Ewald Real will be calculated.
  !! \param spindex Species index list.
  !! \param splist Element symbol for every species.
  !! \param coordinates Coordinates for every atom in the system.
  !! \param charges Charges for every atom in the system.
  !! \param hubbardu Hubbard parameter U for every species.
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param coul_acc Coulomb accuracy.
  !! \param coul_forces_k Coulombic forces (k space contribution)
  !! \param coul_pot_k Coulombic potential (k space contribution)
  subroutine get_ewald_recip(spindex,splist,coordinates,charges,hubbardu,&
       lattice_vectors,recip_vectors,volr,coul_acc,coul_forces_k,coul_pot_k);

    character(2), intent(in)             ::  splist(:)
    integer                              ::  i, ik, l, lmax, m
    integer                              ::  mmax, mmin, n, nats
    integer                              ::  nmax, nmin
    integer                              ::  id, nt
    integer, intent(in)                  ::  spindex(:)
    real(dp)                             ::  a2xa3(3), calpha, calpha2, corrfact
    real(dp)                             ::  cossum, cossum2, coulcut, coulcut2
    real(dp)                             ::  coulombv, coulvol, dot, fcoul(3)
    real(dp)                             ::  force, fourcalpha2, k(3), k2
    real(dp)                             ::  kcutoff, kcutoff2, keconst, kepref
    real(dp)                             ::  l11, l12, l13, m21, knorm
    real(dp)                             ::  m22, m23, pi, prefactor
    real(dp)                             ::  previr, relperm, sinsum, sinsum2
    real(dp)                             ::  sqrtp, sqrtpi, sqrtx, tfact
    real(dp)                             ::  timeratio
    real(dp), allocatable                ::  coslist(:), sinlist(:)
    real(dp), allocatable, intent(inout) ::  coul_forces_k(:,:), coul_pot_k(:)
    real(dp), intent(in)                 ::  charges(:), coordinates(:,:), coul_acc, hubbardu(:)
    real(dp), intent(in)                 ::  lattice_vectors(:,:), recip_vectors(:,:), volr
    real(dp), allocatable                ::  k1_list(:),k2_list(:),k3_list(:),ksq_list(:)
    integer                              ::  Nk

    timeratio = 10.0_dp  !Estimated ration between real & k space
    pi = 3.14159265358979323846264338327950_dp

    nats = size(charges,dim=1)

    if(.not.allocated(coul_forces_k))allocate(coul_forces_k(3,nats))
    if(.not.allocated(coul_pot_k))allocate(coul_pot_k(nats))

    coul_pot_k = 0.0_dp
    coul_forces_k = 0.0_dp

    allocate(sinlist(nats))
    allocate(coslist(nats))

    sqrtx = sqrt(-log(coul_acc))
    calpha = sqrt(pi)*((timeratio*nats/(volr**2))**(1.0_dp/6.0_dp))
    coulcut = sqrtx/calpha
    calpha2 = calpha*calpha
    if (coulcut > 50.0_dp) then
       coulcut = 50.0_dp
       calpha = sqrtx/coulcut
    endif
    coulcut2 = coulcut*coulcut
    kcutoff = 2.0_dp*calpha*sqrtx
    kcutoff2 = kcutoff*kcutoff
    calpha2 = calpha*calpha
    fourcalpha2 = 4.0_dp*calpha2

    knorm = prg_norm2(recip_vectors(1,:))
    lmax = floor(kcutoff / knorm)
    knorm = prg_norm2(recip_vectors(2,:))
    mmax = floor(kcutoff / knorm)
    knorm = prg_norm2(recip_vectors(3,:))
    nmax = floor(kcutoff / knorm)

    !! The 14.399 factor corresponds to 1/(4*pi*epsilon0) in eV*Ang
    relperm = 1.0_dp
    keconst = 14.3996437701414_dp*relperm

    sqrtpi = sqrt(pi);

    coul_forces_k(:,:) = 0.0_dp
    coul_pot_k(:) = 0.0_dp
    sinlist(:) = 0.0_dp
    coslist(:) = 0.0_dp

    call get_k_lists(recip_vectors,k1_list,k2_list,k3_list,ksq_list,Nk,kcutoff)

    !$omp parallel do default(none) &
    !$omp shared(nats,charges,k1_list,k2_list,k3_list,ksq_list) &
    !$omp shared(coordinates,Nk, keconst, volr, calpha2, pi) &
    !$omp private(i,ik,dot,force,sinlist,coslist,cossum,sinsum) &
    !$omp private(cossum2,sinsum2, prefactor,kepref, previr) &
    !$omp reduction(+:coul_forces_k,coul_pot_k)
    do ik = 1, Nk
       prefactor = 8*pi*exp(-ksq_list(ik)/(4*calpha2))/(volr*ksq_list(ik))
       previr = (2/ksq_list(ik)) + (2/(4*calpha2))

       ! doing the sin and cos sums

       cossum = 0.0_dp
       sinsum = 0.0_dp

       do i = 1,nats
          dot = k1_list(ik)*coordinates(1,i) + k2_list(ik)*coordinates(2,i) &
               + k3_list(ik)*coordinates(3,i)
          ! we re-use these in the next loop...
          sinlist(i) = sin(dot)
          coslist(i) = cos(dot)
          cossum = cossum + charges(i)*coslist(i)
          sinsum = sinsum + charges(i)*sinlist(i)
       enddo

       cossum2 = cossum*cossum
       sinsum2 = sinsum*sinsum

       ! add up energy and force contributions

       kepref = keconst*prefactor

       do i = 1,nats
          coul_pot_k(i) = coul_pot_k(i) + &
               kepref*(coslist(i)*cossum + sinlist(i)*sinsum)
          force = kepref * charges(i) * &
               (sinlist(i)*cossum - coslist(i)*sinsum)

          coul_forces_k(1,i) = coul_forces_k(1,i) + force*k1_list(ik)
          coul_forces_k(2,i) = coul_forces_k(2,i) + force*k2_list(ik)
          coul_forces_k(3,i) = coul_forces_k(3,i) + force*k3_list(ik)
       enddo

       kepref = keconst*prefactor * (cossum2 + sinsum2)

    enddo
    !$omp end parallel do

    ! Point self energy
    corrfact = 2*keconst*calpha/sqrtpi;
    coul_pot_k = coul_pot_k - corrfact*charges;

  end subroutine get_ewald_recip


  !> This routine constructs the k points list, e.g., for get_ewald_recip
  !! \param recip_vectors Reciprocal lattice vectors
  !! \param k1_list list of k(1) values
  !! \param k2_list list of k(2) values
  !! \param k3_list list of k(3) values
  !! \param ksq_list list of k2 (k squared) values
  !! \param Nk number of k values in lists (allocated length is longer, due to cutoff)
  !! \param coul_acc Coulomb accuracy.
  !! \param kcutoff Cutoff value for k
  subroutine get_k_lists(recip_vectors,k1_list,k2_list,k3_list,ksq_list,Nk,kcutoff)

    integer                              ::  i, l, lmax, m
    integer                              ::  mmax, mmin, n
    integer                              ::  nmax, nmin
    real(dp), intent(in)                 ::  kcutoff
    real(dp)                             ::  kcutoff2
    real(dp)                             ::  l11, l12, l13, m21
    real(dp)                             ::  m22, m23, pi, knorm
    real(dp)                             ::  k(3),k2
    real(dp), allocatable, intent(inout) ::  k1_list(:),k2_list(:),k3_list(:),ksq_list(:)
    integer, intent(out)                 ::  Nk
    real(dp), intent(in)                 ::  recip_vectors(:,:)

    knorm =  prg_norm2(recip_vectors(1,:))
    lmax = floor(kcutoff / knorm)
    knorm = prg_norm2(recip_vectors(2,:))
    mmax = floor(kcutoff / knorm)
    knorm = prg_norm2(recip_vectors(3,:))
    nmax = floor(kcutoff / knorm)

    !! Maximum value of Nk
    Nk = (2*lmax+1)*(2*mmax+1)*(2*nmax+1)

    !! Use max value for allocation of lists
    if(.not.allocated(k1_list))allocate(k1_list(Nk))
    if(.not.allocated(k2_list))allocate(k2_list(Nk))
    if(.not.allocated(k3_list))allocate(k3_list(Nk))
    if(.not.allocated(ksq_list))allocate(ksq_list(Nk))

    kcutoff2 = kcutoff * kcutoff

    i = 1

    do l = 0,lmax

       if (l == 0) then
          mmin = 0
       else
          mmin = -mmax
       endif

       l11 = l*recip_vectors(1,1)
       l12 = l*recip_vectors(1,2)
       l13 = l*recip_vectors(1,3)

       do m = mmin,mmax

          nmin = -nmax

          if(l .eq. 0 .and. m .eq. 0) nmin = 1

          m21 = l11 + m*recip_vectors(2,1)
          m22 = l12 + m*recip_vectors(2,2)
          m23 = l13 + m*recip_vectors(2,3)
          do n = nmin,nmax
             k(1) = m21 + n*recip_vectors(3,1)
             k(2) = m22 + n*recip_vectors(3,2)
             k(3) = m23 + n*recip_vectors(3,3)
             k2 = k(1)*k(1) + k(2)*k(2) + k(3)*k(3)
             if (k2 <= kcutoff2)then
                !! Should check whether i > size(k1), etc here
                k1_list(i) = k(1)
                k2_list(i) = k(2)
                k3_list(i) = k(3)
                ksq_list(i) = k2
                i = i + 1
             end if
          end do
       end do
    end do
    Nk = i-1

  end subroutine get_k_lists

end module coulomb_latte_mod
