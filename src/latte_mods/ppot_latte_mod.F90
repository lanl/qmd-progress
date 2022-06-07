!> A module to compute the pair potential contribution to the forces and energy.
!! @ingroup LATTE
!! \brief This module will be used to compute pair potential the forces and energy contribution
!!
module ppot_latte_mod

  use prg_openfiles_mod
  use bml
  use tbparams_latte_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: low = 8

  public :: get_PairPot_contrib, get_PairPot_contrib_int

contains

  !> This routine computes the forces and energy from the pair potentials.
  !! \param coords System coordinates.
  !! \param lattice_vectors Lattice vectors.
  !! \param spindex Index of species.
  !! \param ppot Pair potential structure.
  !! \param PairForces Pair potential forces.
  !! \param ERep Repulsive energy.
  !!
  subroutine get_PairPot_contrib(coords,lattice_vectors,spindex,ppot,PairForces,ERep)
    implicit none
    integer                              ::  i, ii, j, jj
    integer                              ::  nats
    integer                              ::  nr_shift_X, nr_shift_Y, nr_shift_Z
    integer, intent(in)                  ::  spindex(:)
    real(dp)                             ::  CUTPHI, DC(3), DPHI(3), DPOLYNOM
    real(dp)                             ::  EXPTMP, FCUT(3), FORCE, FTMP(3)
    real(dp)                             ::  FUNIV(3), LBox(3), MYR, PHI
    real(dp)                             ::  POLYNOM, PotCoef(16), R1, RCUT
    real(dp)                             ::  RCUT2, RXb, RYb, RZb
    real(dp)                             ::  Ra(3), Rb(3), UNIVPHI, VIRCUT
    real(dp)                             ::  VIRUNIV, dR2, dr, rab(3)
    real(dp), allocatable, intent(inout)  ::  PairForces(:,:)
    real(dp), intent(in)                 ::  coords(:,:), lattice_vectors(:,:)
    real(dp), intent(inout)              ::  ERep
    type(ppot_type), intent(inout)       ::  ppot(:,:)

    write(*,*)"In get_PairPot_contrib ..."

    nats = size(coords,dim=2)
    if(.not.allocated(PairForces))then
      allocate(PairForces(3,nats))
    endif

    PairForces = 0.0_dp

    write(*,*)nats

    UNIVPHI = 0;
    CUTPHI = 0;

    VIRUNIV = 0;
    VIRCUT = 0;

    LBox(1) = lattice_vectors(1,1)
    LBox(2) = lattice_vectors(2,2)
    LBox(3) = lattice_vectors(3,3)

    !$omp parallel do default(none) private(i) &
    !$omp private(FCUT,Ra,Rb,RXb,RYb,RZb,Rab,dR,dR2,DC) &
    !$omp private(POLYNOM,PHI,DPOLYNOM,DPHI,EXPTMP,FTMP,FUNIV,MYR) &
    !$omp private(FORCE,j,jj,ii) &
    !$omp private(PotCoef,R1,RCUT,RCUT2,nr_shift_X,nr_shift_Y,nr_shift_Z) &
    !$omp shared(nats,coords,spindex,ppot,lbox) &
    !$omp shared(PairForces) &
    !$omp reduction (+:UNIVPHI,CUTPHI)
    do i = 1, nats
      FUNIV = 0.0_dp
      FCUT = 0.0_dp
      Ra(1) = coords(1,i); Ra(2) = coords(2,i); Ra(3) = coords(3,i)
      ii=spindex(i)
      do j = 1,nats
        if(i.ne.j)then

          jj=spindex(j)

          PotCoef = ppot(ii,jj)%potparams;

          R1 = PotCoef(9)
          RCUT = PotCoef(10)
          RCUT2 = RCUT*RCUT;

          RXb = coords(1,j); RYb = coords(2,j); RZb = coords(3,j);
          do nr_shift_X = -1,1  ! Periodic BC shifts in X, Y and Z. Costs a lot extra!
            do nr_shift_Y = -1,1
              do nr_shift_Z = -1,1

                Rb(1) = RXb + nr_shift_X*LBox(1); ! Shifts for PBC
                Rb(2) = RYb + nr_shift_Y*LBox(2);
                Rb(3) = RZb + nr_shift_Z*LBox(3);
                Rab = Rb-Ra;  ! OBS b - a !!!
                dR2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
                dR = sqrt(dR2)
                !!                dR = norm2(Rab) ;
                !!
                !!                dR2 = dR*dR;

                if (dR < RCUT)then

                  DC = Rab/dR;

                  if (dR < R1)then

                    POLYNOM = dR*(PotCoef(2) + dR*(PotCoef(3) + dR*(PotCoef(4) + dR*PotCoef(5))));
                    PHI = PotCoef(1)*exp(POLYNOM);
                    DPOLYNOM = PotCoef(2) + dR*(2*PotCoef(3) + dR*(3*PotCoef(4) + 4*PotCoef(5)*dR));
                    DPHI = -DC*PHI*DPOLYNOM;
                    EXPTMP = PotCoef(6)*exp( PotCoef(7)*(dR - PotCoef(8)) );

                    UNIVPHI = UNIVPHI + PHI + EXPTMP;
                    FTMP = DC*PotCoef(7)*EXPTMP;
                    FUNIV = FUNIV - DPHI + FTMP;

                  else

                    MYR = dR - R1;
                    CUTPHI =  CUTPHI + PotCoef(11) + MYR*(PotCoef(12) + MYR*(PotCoef(13 ) + MYR*(PotCoef(14 ) + MYR*(PotCoef(15 ) + MYR*PotCoef(16 )))));
                    FORCE = PotCoef(12 )  + MYR*(2*PotCoef(13 ) + MYR*(3*PotCoef(14 ) + MYR*(4*PotCoef(15 ) + MYR*5*PotCoef(16 ))));
                    FCUT = FCUT + DC*FORCE;

                  endif
                endif
              enddo
            enddo
          enddo

        endif

      enddo
      PairForces(:,i) = FUNIV + FCUT;
    enddo
    ! $omp end parallel do

    ERep = 0.5*(UNIVPHI + CUTPHI);

  end subroutine get_PairPot_contrib

  !> This routine computes the forces and energy from the pair potentials.
  !! \param coords System coordinates.
  !! \param lattice_vectors Lattice vectors.
  !! \param spindex Index of species.
  !! \param ppot Pair potential structure.
  !! \param PairForces Pair potential forces.
  !! \param ERep Repulsive energy.
  !!
  subroutine get_PairPot_contrib_int(coords,lattice_vectors,nnIx,nnIy,&
       nnIz,nrnnlist,nnType,spindex,ppot,PairForces,ERep)
    implicit none
    integer                              ::  i, ii, j, jj
    integer                              ::  nats, nni
    integer                              ::  nr_shift_X, nr_shift_Y, nr_shift_Z
    integer(kind=low), intent(in)               ::  nnIx(:,:),nnIy(:,:),nnIz(:,:)
    integer, intent(in)                  ::  spindex(:)
    integer, intent(in)                  ::  nrnnlist(:), nnType(:,:)
    real(dp)                             ::  CUTPHI, DC(3), DPHI(3), DPOLYNOM
    real(dp)                             ::  EXPTMP, FCUT(3), FORCE, FTMP(3)
    real(dp)                             ::  FUNIV(3), Lx, Ly, Lz, MYR, PHI
    real(dp)                             ::  POLYNOM, PotCoef(16), R1, RCUT
    real(dp)                             ::  RCUT2, RXb, RYb, RZb
    real(dp)                             ::  Ra(3), Rb(3), UNIVPHI, VIRCUT
    real(dp)                             ::  VIRUNIV, dR2, dr, rab(3)
    real(dp), allocatable, intent(inout)  ::  PairForces(:,:)
    real(dp), intent(in)                 ::  coords(:,:), lattice_vectors(:,:)
    real(dp), intent(inout)              ::  ERep
    type(ppot_type), intent(inout)       ::  ppot(:,:)

    write(*,*)"In get_PairPot_contrib ..."

    nats = size(coords,dim=2)
    if(.not.allocated(PairForces))then
      allocate(PairForces(3,nats))
    endif

    PairForces = 0.0_dp

    write(*,*)nats

    UNIVPHI = 0;
    CUTPHI = 0;

    VIRUNIV = 0;
    VIRCUT = 0;

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)

    !$omp parallel do default(none) private(i) &
    !$omp private(FCUT,Ra,Rb,RXb,RYb,RZb,Rab,dR,dR2,DC) &
    !$omp private(POLYNOM,PHI,DPOLYNOM,DPHI,EXPTMP,FTMP,FUNIV,MYR) &
    !$omp private(FORCE,j,jj,ii,nni) &
    !$omp private(PotCoef,R1,RCUT,RCUT2,nr_shift_X,nr_shift_Y,nr_shift_Z) &
    !$omp shared(nats,coords,spindex,ppot,Lx,Ly,Lz) &
    !$omp shared(PairForces,nnIx,nnIy,nnIz,nrnnlist,nnType) &
    !$omp reduction (+:UNIVPHI,CUTPHI)
    do i = 1, nats
      FUNIV = 0.0_dp
      FCUT = 0.0_dp
      Ra(1) = coords(1,i); Ra(2) = coords(2,i); Ra(3) = coords(3,i)
      ii=spindex(i)

      do nni = 1,nrnnlist(i)
        j = nnType(nni,i)

        if(i.ne.j)then

          jj=spindex(j)

          PotCoef = ppot(ii,jj)%potparams;

          R1 = PotCoef(9)
          RCUT = PotCoef(10)
          RCUT2 = RCUT*RCUT;

          Rb(1) = coords(1,j)
          Rb(2) = coords(2,j)
          Rb(3) = coords(3,j)

          !  Rb(1) = Rb(1) + nnIx(nni,i)*Lx; ! Shifts for PBC
          !  Rb(2) = Rb(2) + nnIy(nni,i)*Ly;
          !  Rb(3) = Rb(3) + nnIz(nni,i)*Lz;


          rab(1) = modulo((Rb(1) - Ra(1) + Lx/2.0_dp),Lx) - Lx/2.0_dp
          rab(2) = modulo((Rb(2) - Ra(2) + Ly/2.0_dp),Ly) - Ly/2.0_dp
          rab(3) = modulo((Rb(3) - Ra(3) + Lz/2.0_dp),Lz) - Lz/2.0_dp


          dR2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
          dR = sqrt(dR2)

          if (dR < RCUT)then

            DC = Rab/dR;

            if (dR < R1)then

              POLYNOM = dR*(PotCoef(2) + dR*(PotCoef(3) + dR*(PotCoef(4) + dR*PotCoef(5))));
              PHI = PotCoef(1)*exp(POLYNOM);
              DPOLYNOM = PotCoef(2) + dR*(2*PotCoef(3) + dR*(3*PotCoef(4) + 4*PotCoef(5)*dR));
              DPHI = -DC*PHI*DPOLYNOM;
              EXPTMP = PotCoef(6)*exp( PotCoef(7)*(dR - PotCoef(8)) );

              UNIVPHI = UNIVPHI + PHI + EXPTMP;
              FTMP = DC*PotCoef(7)*EXPTMP;
              FUNIV = FUNIV - DPHI + FTMP;

            else

              MYR = dR - R1;
              CUTPHI =  CUTPHI + PotCoef(11) + MYR*(PotCoef(12) + MYR*(PotCoef(13 ) + MYR*(PotCoef(14 ) + MYR*(PotCoef(15 ) + MYR*PotCoef(16 )))));
              FORCE = PotCoef(12 )  + MYR*(2*PotCoef(13 ) + MYR*(3*PotCoef(14 ) + MYR*(4*PotCoef(15 ) + MYR*5*PotCoef(16 ))));
              FCUT = FCUT + DC*FORCE;

            endif
          endif

        endif

      enddo
      PairForces(:,i) = FUNIV + FCUT;
    enddo
    ! $omp end parallel do

    ERep = 0.5*(UNIVPHI + CUTPHI);

  end subroutine get_PairPot_contrib_int

end module ppot_latte_mod
