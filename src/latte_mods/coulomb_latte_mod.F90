!> Module to compute the coulombic contribution of the Forces.
!! \ingroup LATTE
!!
module coulomb_latte_mod

  use bml
  use prg_extras_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: low = kind(100)
  public :: get_ewald_real, get_ewald_recip, get_coulcut
  public :: get_ewald_list_real, get_ewald_list_real_dcalc
  public :: get_ewald_list_real_dcalc_vect

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
    integer,allocatable, intent(in)    ::  nnIx(:,:),nnIy(:,:),nnIz(:,:)
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

          if(rmod > coulcut)then
            rab(1) = modulo((Rb(1) - Ra(1) + Lx/2.0_dp),Lx) - Lx/2.0_dp
            rab(2) = modulo((Rb(2) - Ra(2) + Ly/2.0_dp),Ly) - Ly/2.0_dp
            rab(3) = modulo((Rb(3) - Ra(3) + Lz/2.0_dp),Lz) - Lz/2.0_dp
          endif
        endif

        dr = norm2(rab)

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
      !$omp critical
      coul_forces_r(:,i) = fcoul
      coul_pot_r(i) = coulombv
      !$omp end critical
    enddo
    !$omp end parallel do


    coul_pot_r = keconst*coul_pot_r

  end subroutine get_ewald_list_real_dcalc

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
  subroutine get_ewald_list_real_dcalc_vect(spindex,splist,coordinates,charges,hubbardu&
       ,lattice_vectors,volr,coul_acc,timeratio,nnIx,nnIy,&
       nnIz,nrnnlist,nnType&
       ,coul_forces_r,coul_pot_r)

    character(2), intent(in)             ::  splist(:)
    integer                              ::  atomi, i, j, nats
    integer                              ::  nnI
    integer, intent(in)                  ::  spindex(:)
    integer,allocatable, intent(in)    ::  nnIx(:,:),nnIy(:,:),nnIz(:,:)
    real(dp)                             ::  a2xa3(3), ca, calpha, calpha2
    real(dp)                             ::  coul_acc, coulcut, coulcut2, coulombv
    real(dp)                             ::  dcoulombv,dforce
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
    integer                              ::  maxn
    real(dp), allocatable                ::  dx_vec(:), dy_vec(:),dz_vec(:), dr_vec(:), dr2_vec(:)
    real(dp), allocatable                ::  qj_vec(:), qij_vec(:), ti_list(:), ti2_list(:), ti4_list(:), ti6_list(:)
    real(dp), allocatable                ::  tj_list(:), tj2_list(:), tj4_list(:), tj6_list(:)
    real(dp), allocatable                ::  ti2mtj2_vec(:), ti2mtj2_2_vec(:),ti2mtj2_3_vec(:) 
    real(dp), allocatable                ::  sa_vec(:), sb_vec(:), sc_vec(:)
    real(dp), allocatable                ::  sd_vec(:), se_vec(:), sf_vec(:)
    real(dp), allocatable                ::  ssa_vec(:), ssb_vec(:), ssc_vec(:), ssd_vec(:), sse_vec(:)
    real(dp), allocatable                ::  f_vec(:), f1_vec(:), f2_vec(:)
    real(dp), allocatable                ::  e_vec(:), e1_vec(:), e2_vec(:)
    real(dp), allocatable                ::  expti_vec(:), exptj_vec(:), ca_vec(:)
    logical, allocatable                 ::  neigh_mask(:), same_sp_mask(:), distance_mask(:)
    logical                              ::  use_modulo_trick

    !Estimated ration between real & k space
    pi = 3.14159265358979323846264338327950_dp

    nats = size(charges,dim=1)

    maxn = maxval(nrnnlist)

    relperm = 1.0_dp;

    !! The 14.399 factor corresponds to 1/(4*pi*epsilon0) in eV*Ang
    keconst = 14.3996437701414_dp*relperm;
    tfact  = 16.0_dp/(5.0_dp*keconst);

    sqrtpi = sqrt(pi);


    if(allocated(sa_vec))then
       if(size(sa_vec).lt.maxn)then
          deallocate(dx_vec)
          deallocate(dy_vec)
          deallocate(dz_vec)
          deallocate(dr_vec)
          deallocate(dr2_vec)
          deallocate(qj_vec)
          deallocate(qij_vec)
          deallocate(neigh_mask)
          deallocate(same_sp_mask)
          deallocate(distance_mask)
          deallocate(ti_list)
          deallocate(ti2_list)
          deallocate(ti4_list)
          deallocate(ti6_list)
          deallocate(tj_list)
          deallocate(tj2_list)
          deallocate(tj4_list)
          deallocate(tj6_list)
          deallocate(ti2mtj2_vec)
          deallocate(ti2mtj2_2_vec)
          deallocate(ti2mtj2_3_vec)
          deallocate(sa_vec)
          deallocate(sb_vec)
          deallocate(sc_vec)
          deallocate(sd_vec)
          deallocate(se_vec)
          deallocate(sf_vec)
          deallocate(ssa_vec)
          deallocate(ssb_vec)
          deallocate(ssc_vec)
          deallocate(ssd_vec)
          deallocate(sse_vec)
          deallocate(e_vec)
          deallocate(e1_vec)
          deallocate(e2_vec)
          deallocate(f_vec)
          deallocate(f1_vec)
          deallocate(f2_vec)
          deallocate(ca_vec)
          deallocate(expti_vec)
          deallocate(exptj_vec)
       endif
    endif
    
    if(.not.allocated(sa_vec))then
       allocate(dx_vec(maxn))
       allocate(dy_vec(maxn))
       allocate(dz_vec(maxn))
       allocate(dr_vec(maxn))
       allocate(dr2_vec(maxn))
       allocate(qj_vec(maxn))
       allocate(qij_vec(maxn))
       allocate(neigh_mask(maxn))
       allocate(same_sp_mask(maxn))
       allocate(distance_mask(maxn))
       allocate(tj_list(maxn))
       allocate(tj2_list(maxn))
       allocate(tj4_list(maxn))
       allocate(tj6_list(maxn))
       allocate(ti2mtj2_vec(maxn))
       allocate(ti2mtj2_2_vec(maxn))
       allocate(ti2mtj2_3_vec(maxn))
       allocate(sa_vec(maxn))
       allocate(sb_vec(maxn))
       allocate(sc_vec(maxn))
       allocate(sd_vec(maxn))
       allocate(se_vec(maxn))
       allocate(sf_vec(maxn))
       allocate(ssa_vec(maxn))
       allocate(ssb_vec(maxn))
       allocate(ssc_vec(maxn))
       allocate(ssd_vec(maxn))
       allocate(sse_vec(maxn))
       allocate(e_vec(maxn))
       allocate(e1_vec(maxn))
       allocate(e2_vec(maxn))
       allocate(f_vec(maxn))
       allocate(f1_vec(maxn))
       allocate(f2_vec(maxn))
       allocate(ca_vec(maxn))
       allocate(expti_vec(maxn))
       allocate(exptj_vec(maxn))
    endif

    if(.not.allocated(ti_list))then
       allocate(ti_list(nats))
       allocate(ti2_list(nats))
       allocate(ti4_list(nats))
       allocate(ti6_list(nats))
       !$omp parallel do default(none) private(i) &
       !$omp shared(nats,spindex,ti_list,hubbardu)
       do i=1,nats
          ti_list(i) = hubbardu(spindex(i))
       enddo
       !$omp end parallel do
    endif
    ti_list(:) = tfact * ti_list(:)
    ti2_list(:) = ti_list(:) * ti_list(:)
    ti4_list(:) = ti2_list(:) * ti2_list(:)
    ti6_list(:) = ti4_list(:) * ti2_list(:)
    
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

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)

    if(allocated(nnIx))then
       use_modulo_trick = .false.
    else
       use_modulo_trick = .true.
    endif
    
    if(.not.allocated(coul_forces_r))allocate(coul_forces_r(3,nats))
    if(.not.allocated(coul_pot_r))allocate(coul_pot_r(nats))
    
    !$omp parallel do default(none) private(i,j,nnI,tj_list,tj2_list,tj4_list,tj6_list) &
    !$omp shared(nats,nrnnlist,charges,nnType,coordinates,Lx,Ly,Lz) &
    !$omp private(qj_vec,qij_vec,dx_vec,dy_vec,dz_vec,dr_vec,dr2_vec) &
    !$omp private(ti2mtj2_vec,ti2mtj2_2_vec,ti2mtj2_3_vec) &
    !$omp private(sa_vec,sb_vec,sc_vec,sd_vec,se_vec,sf_vec) &
    !$omp private(ssa_vec,ssb_vec,ssc_vec,ssd_vec,sse_vec) &
    !$omp private(ca_vec,expti_vec,exptj_vec) &
    !$omp private(e_vec,e1_vec,e2_vec) &
    !$omp private(f_vec,f1_vec,f2_vec) &
    !$omp shared(calpha,calpha2,sqrtpi,keconst,use_modulo_trick,coulcut) &
    !$omp shared(nnIx,nnIy,nnIz) &
    !$omp private(neigh_mask,same_sp_mask,distance_mask) &
    !$omp shared(coul_forces_r,coul_pot_r) &
    !$omp shared(ti_list,ti2_list,ti4_list,ti6_list)
    do i=1,nats
       neigh_mask = .false.
       neigh_mask(1:nrnnlist(i)) = .true.
       dx_vec = 0.0_dp
       dy_vec = 0.0_dp
       dz_vec = 0.0_dp
       qj_vec = 0.0_dp
       ti2mtj2_vec = 0.0_dp
       ti2mtj2_2_vec = 0.0_dp
       ti2mtj2_3_vec = 0.0_dp
       tj_list(:) = 0.0_dp
       tj2_list(:) = 0.0_dp
       tj4_list(:) = 0.0_dp
       tj6_list(:) = 0.0_dp
       do nnI=1,nrnnlist(i)
          j = nnType(nnI,i)
          qj_vec(nnI) = charges(j)
          dx_vec(nnI) = coordinates(1,j)
          dy_vec(nnI) = coordinates(2,j)
          dz_vec(nnI) = coordinates(3,j)          
          ti2mtj2_vec(nnI) = ti2_list(j)
          tj_list(nnI) = ti_list(j)
          tj2_list(nnI) = ti2_list(j)
          tj4_list(nnI) = ti4_list(j)
          tj6_list(nnI) = ti6_list(j)
       enddo
       qij_vec(:) = qj_vec(:) * charges(i)
       if(use_modulo_trick)then
          dx_vec(:) = modulo(coordinates(1,i) - dx_vec(:) + Lx/2.0_dp,Lx) - Lx/2.0_dp
          dy_vec(:) = modulo(coordinates(2,i) - dy_vec(:) + Ly/2.0_dp,Ly) - Ly/2.0_dp
          dz_vec(:) = modulo(coordinates(3,i) - dz_vec(:) + Lz/2.0_dp,Lz) - Lz/2.0_dp
       else
          dx_vec(:) = coordinates(1,i) - dx_vec(:) - nnIx(nni,i)*Lx
          dy_vec(:) = coordinates(2,i) - dy_vec(:) - nnIy(nni,i)*Ly
          dz_vec(:) = coordinates(3,i) - dz_vec(:) - nnIz(nni,i)*Lz
       endif
       dr2_vec(:) = dx_vec(:)*dx_vec(:)+dy_vec(:)*dy_vec(:)+dz_vec(:)*dz_vec(:)
       dr_vec(:) = sqrt(dr2_vec(:))
       distance_mask(:) = (dr_vec(:) <= coulcut).and.(dr_vec(:) > 1.0D-12)
       dr2_vec(:) = merge(dr2_vec(:),1.0_dp,distance_mask(:))
       dr_vec(:) = merge(dr_vec(:),1.0_dp,distance_mask(:))
       ti2mtj2_vec(:) = ti2_list(i) - ti2mtj2_vec(:)
       same_sp_mask(:) = ti2mtj2_vec(:).eq.0.0_dp
       ti2mtj2_vec(:) = merge(1.0_dp,ti2mtj2_vec(:),same_sp_mask(:))
       ti2mtj2_2_vec(:) = ti2mtj2_vec(:)*ti2mtj2_vec(:)
       ti2mtj2_3_vec(:) = ti2mtj2_2_vec(:)*ti2mtj2_vec(:)
       ti2mtj2_2_vec(:) = 2.0_dp*ti2mtj2_2_vec(:)
       sa_vec(:) = ti_list(i)
       sb_vec(:) = tj4_list(:)*ti_list(i)/ti2mtj2_2_vec(:)
       sc_vec(:) = + (tj6_list(:) - 3.0_dp*tj4_list(:)*ti2_list(i))/ti2mtj2_3_vec(:)
       sd_vec(:) = tj_list(:)
       se_vec(:) = ti4_list(i)*tj_list(:)/ti2mtj2_2_vec(:)
       sf_vec(:) = - (ti6_list(i) - 3.0_dp*ti4_list(i)*tj2_list(:))/ti2mtj2_3_vec(:)
       ssa_vec(:) = ti_list(i)
       ssb_vec(:) = ti2_list(i)*ti_list(i)/48.0_dp
       ssc_vec(:) = 3.0_dp*ti2_list(i)/16.0_dp
       ssd_vec(:) = 11.0_dp*ti_list(i)/16.0_dp
       sse_vec(:) = 1.0_dp
       expti_vec(:) = exp(-ti_list(i)*dr_vec(:))
       exptj_vec(:) = exp(-tj_list(:)*dr_vec(:))
       e1_vec(:) = - qj_vec(:)*expti_vec(:)* &
            (ssb_vec(:)*dr2_vec(:) + &
             ssc_vec(:)*dr_vec(:) + ssd_vec(:) + sse_vec(:)/dr_vec(:))
       f1_vec(:) = keconst*qij_vec(:)*expti_vec(:)* &
            ((sse_vec(:)/dr2_vec(:) - 2.0_dp*ssb_vec(:)*dr_vec(:) - ssc_vec(:)) + &
              ssa_vec(:)* &
              (ssb_vec(:)*dr2_vec(:) + ssc_vec(:)*dr_vec(:) + ssd_vec(:) + sse_vec(:)/dr_vec(:)))
       e2_vec(:) = - qj_vec(:)* &
            (expti_vec(:)*(sb_vec(:) - sc_vec(:)/dr_vec(:)) + &
             exptj_vec(:)*(se_vec(:) - sf_vec(:)/dr_vec(:)))
       f2_vec(:) = keconst*qij_vec(:) * &
           (expti_vec(:)*(sa_vec(:)*(sb_vec(:) - sc_vec(:)/dr_vec(:)) - sc_vec(:)/dr2_vec(:)) + &
            exptj_vec(:)*(sd_vec(:)*(se_vec(:) - sf_vec(:)/dr_vec(:)) - sf_vec(:)/dr2_vec(:)))
       ca_vec(:) = erfc(abs(calpha*dr_vec(:)))/dr_vec(:)
       e_vec(:) = qj_vec(:)*ca_vec(:) + merge(e1_vec(:),e2_vec(:),same_sp_mask(:))
       ca_vec(:) = ca_vec(:) + 2.0_dp*calpha*exp(-calpha2*dr2_vec(:))/sqrtpi
       f_vec(:) = -keconst*qij_vec(:)*ca_vec(:)/dr_vec(:) + merge(f1_vec(:),f2_vec(:),same_sp_mask(:))
       coul_pot_r(i) = keconst*sum(e_vec(:),neigh_mask(:).and.distance_mask(:))
       coul_forces_r(1,i) = - sum(f_vec(:)*dx_vec(:)/dr_vec(:),neigh_mask(:).and.distance_mask(:))
       coul_forces_r(2,i) = - sum(f_vec(:)*dy_vec(:)/dr_vec(:),neigh_mask(:).and.distance_mask(:))
       coul_forces_r(3,i) = - sum(f_vec(:)*dz_vec(:)/dr_vec(:),neigh_mask(:).and.distance_mask(:))
    enddo    
    !$omp end parallel do
    
  end subroutine get_ewald_list_real_dcalc_vect

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
