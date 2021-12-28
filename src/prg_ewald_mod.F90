! Ewald sum routines for kernel calculation
module prg_ewald_mod

  use bml
  use prg_timer_mod
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: Ewald_Real_Space_Single
  public :: Ewald_Real_Space_Single_latte
  public :: Ewald_Real_Space
  public :: Ewald_Real_Space_latte
  public :: Ewald_Real_Space_Test
  public :: Ewald_Real_Space_Matrix_latte
  public :: Ewald_k_space_latte_single
  public :: Ewald_k_space_latte
  public :: Ewald_k_space_Test
  public :: Ewald_k_space_Matrix_latte

contains

  !> Find Coulomb potential on site I from single charge at site J
  subroutine  Ewald_Real_Space_Single_latte(COULOMBV,I,RXYZ,Box,Nr_elem, &
       DELTAQ,J,U,Element_Pointer,Nr_atoms,COULACC,HDIM,Max_Nr_Neigh)

    implicit none

    integer,    parameter     :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, Nr_elem, HDIM, Max_Nr_Neigh, I, J, Element_Pointer(Nr_atoms)
    real(PREC), intent(in)    :: COULACC, DELTAQ(Nr_atoms)
    real(PREC)                :: TFACT, RELPERM, KECONST
    real(PREC), intent(in)    :: RXYZ(3,Nr_atoms), Box(3,3)
    real(PREC), intent(in)    :: U(Nr_elem)
    real(PREC)                :: COULCUT, COULCUT2
    real(PREC), intent(out)   :: COULOMBV
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: NUMREP_ERFC, CA, FORCE, EXPTI, EXPTJ
    real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
    real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2

    integer                   :: K, ccnt,l,m,n

    COULVOL = Box(1,1)*Box(2,2)*Box(3,3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0
    COULCUT = 12.D0
    CALPHA = SQRTX/COULCUT
    COULCUT2 = COULCUT*COULCUT
    CALPHA2 = CALPHA*CALPHA

    RELPERM = ONE
    KECONST = 14.3996437701414D0*RELPERM
    TFACT  = 16.0D0/(5.0D0*KECONST)

    COULOMBV = ZERO

    TI = TFACT*U(Element_Pointer(I))
    TI2 = TI*TI
    TI3 = TI2*TI
    TI4 = TI2*TI2
    TI6 = TI4*TI2

    SSA = TI
    SSB = TI3/48.D0
    SSC = 3.D0*TI2/16.D0
    SSD = 11.D0*TI/16.D0
    SSE = 1.D0

    Ra(1) = RXYZ(1,I)
    Ra(2) = RXYZ(2,I)
    Ra(3) = RXYZ(3,I)

    do k = -1,1
      do m = -1,1
        do l = -1,1

          Rb(1) = RXYZ(1,J)+k*box(1,1)
          Rb(2) = RXYZ(2,J)+m*box(2,2)
          Rb(3) = RXYZ(3,J)+l*box(3,3)
          Rab = Rb-Ra  ! OBS b - a !!!
          dR = norm2(Rab)
          MAGR = dR
          MAGR2 = dR*dR

          if ((dR <= COULCUT).and.(dR > 1e-12)) then

            TJ = TFACT*U(Element_Pointer(J))
            DC = Rab/dR

            Z = abs(CALPHA*MAGR)
            NUMREP_ERFC = erfc(Z)

            CA = NUMREP_ERFC/MAGR
            COULOMBV = COULOMBV + DELTAQ(J)*CA
            ccnt = ccnt + 1
            CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
            EXPTI = exp(-TI*MAGR )

            if (Element_Pointer(I).eq.Element_Pointer(J)) then
              COULOMBV = COULOMBV - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
              ccnt = ccnt + 1
            else
              TJ2 = TJ*TJ
              TJ3 = TJ2*TJ
              TJ4 = TJ2*TJ2
              TJ6 = TJ4*TJ2
              EXPTJ = exp( -TJ*MAGR )
              TI2MTJ2 = TI2 - TJ2
              TJ2MTI2 = -TI2MTJ2
              SA = TI
              SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
              SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
              SD = TJ
              SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
              SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

              COULOMBV = COULOMBV - (DELTAQ(J)*(EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
            endif
          endif
        enddo
      enddo
    enddo

    COULOMBV = KECONST*COULOMBV

  end subroutine Ewald_Real_Space_Single_latte

  subroutine  Ewald_Real_Space_Single(COULOMBV,FCOUL,I,RX,RY,RZ,LBox, &
       DELTAQ,J,U,Element_Type,Nr_atoms,COULACC,TIMERATIO,HDIM,Max_Nr_Neigh)


    implicit none

    integer,    parameter     :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, HDIM, Max_Nr_Neigh, I, J
    real(PREC), intent(in)    :: COULACC, TIMERATIO,DELTAQ(Nr_atoms)
    real(PREC)                :: TFACT, RELPERM, KECONST
    real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
    real(PREC), intent(in)    :: U(Nr_atoms)
    real(PREC)                :: COULCUT, COULCUT2
    character(10), intent(in) :: Element_Type(Nr_atoms)
    real(PREC), intent(out)   :: COULOMBV, FCOUL(3)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: NUMREP_ERFC, CA, FORCE, EXPTI, EXPTJ
    real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
    real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2

    integer                   :: K, ccnt,l,m,n

    COULVOL = LBox(1)*LBox(2)*LBox(3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0
    COULCUT = 12.D0
    CALPHA = SQRTX/COULCUT
    COULCUT2 = COULCUT*COULCUT
    CALPHA2 = CALPHA*CALPHA

    RELPERM = ONE
    KECONST = 14.3996437701414D0*RELPERM
    TFACT  = 16.0D0/(5.0D0*KECONST)

    FCOUL = ZERO
    COULOMBV = ZERO

    TI = TFACT*U(I)
    TI2 = TI*TI
    TI3 = TI2*TI
    TI4 = TI2*TI2
    TI6 = TI4*TI2

    SSA = TI
    SSB = TI3/48.D0
    SSC = 3.D0*TI2/16.D0
    SSD = 11.D0*TI/16.D0
    SSE = 1.D0

    Ra(1) = RX(I)
    Ra(2) = RY(I)
    Ra(3) = RZ(I)

    do k = -1,1
      do m = -1,1
        do l = -1,1

          Rb(1) = RX(J)+k*Lbox(1)
          Rb(2) = RY(J)+m*Lbox(2)
          Rb(3) = RZ(J)+l*Lbox(3)
          Rab = Rb-Ra  ! OBS b - a !!!
          dR = norm2(Rab)
          MAGR = dR
          MAGR2 = dR*dR

          if ((dR <= COULCUT).and.(dR > 1e-12)) then

            TJ = TFACT*U(J)
            DC = Rab/dR

            Z = abs(CALPHA*MAGR)
            NUMREP_ERFC = erfc(Z)

            CA = NUMREP_ERFC/MAGR
            COULOMBV = COULOMBV + DELTAQ(J)*CA
            ccnt = ccnt + 1
            CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
            FORCE = -KECONST*DELTAQ(I)*DELTAQ(J)*CA/MAGR
            EXPTI = exp(-TI*MAGR )

            if (Element_Type(I).eq.Element_Type(J)) then
              COULOMBV = COULOMBV - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
              ccnt = ccnt + 1
              FORCE = FORCE + (KECONST*DELTAQ(I)*DELTAQ(J)*EXPTI)*((SSE/MAGR2 - TWO*SSB*MAGR - SSC) &
                   + SSA*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR))
            else
              TJ2 = TJ*TJ
              TJ3 = TJ2*TJ
              TJ4 = TJ2*TJ2
              TJ6 = TJ4*TJ2
              EXPTJ = exp( -TJ*MAGR )
              TI2MTJ2 = TI2 - TJ2
              TJ2MTI2 = -TI2MTJ2
              SA = TI
              SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
              SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
              SD = TJ
              SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
              SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

              COULOMBV = COULOMBV - (DELTAQ(J)*(EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
              FORCE = FORCE + KECONST*DELTAQ(I)*DELTAQ(J)*((EXPTI*(SA*(SB - (SC/MAGR)) - (SC/MAGR2))) &
                   + (EXPTJ*(SD*(SE - (SF/MAGR)) - (SF/MAGR2))))
            endif

            FCOUL(1) = FCOUL(1) + DC(1)*FORCE
            FCOUL(2) = FCOUL(2) + DC(2)*FORCE
            FCOUL(3) = FCOUL(3) + DC(3)*FORCE
          endif
        enddo
      enddo
    enddo

    COULOMBV = KECONST*COULOMBV

  end subroutine Ewald_Real_Space_Single

  subroutine  Ewald_Real_Space_Matrix_latte(E,RXYZ,Box,U,Element_Pointer,Nr_atoms,COULACC,nebcoul,totnebcoul,HDIM,Max_Nr_Neigh,Nr_Elem)

    implicit none

    integer,    parameter     :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, HDIM, Max_Nr_Neigh, Nr_Elem
    real(PREC), intent(in)    :: COULACC
    real(PREC), intent(out)   :: E(Nr_atoms,Nr_atoms)
    real(PREC)                :: TFACT, RELPERM, KECONST
    real(PREC), intent(in)    :: RXYZ(3,Nr_atoms), Box(3,3)
    real(PREC), intent(in)    :: U(Nr_elem)
    real(PREC)                :: COULCUT, COULCUT2
    integer, intent(in)       :: Element_Pointer(Nr_atoms)
    integer,    intent(in)    :: totnebcoul(Nr_atoms), nebcoul(4,Max_Nr_Neigh,Nr_atoms)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: NUMREP_ERFC, CA, FORCE, EXPTI, EXPTJ
    real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
    real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2
    integer                   :: I,J,K, ccnt, newj, PBCI,PBCJ,PBCK

    COULVOL = Box(1,1)*Box(2,2)*Box(3,3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0
    COULCUT = 12.D0
    CALPHA = SQRTX/COULCUT
    COULCUT2 = COULCUT*COULCUT
    CALPHA2 = CALPHA*CALPHA

    RELPERM = ONE
    KECONST = 14.3996437701414D0*RELPERM
    TFACT  = 16.0D0/(5.0D0*KECONST)

    E = 0.0

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(E,U,Element_Pointer,RXYZ,totnebcoul,nebcoul,coulcut,calpha)

    do I = 1,Nr_atoms

      TI = TFACT*U(Element_Pointer(I))
      TI2 = TI*TI
      TI3 = TI2*TI
      TI4 = TI2*TI2
      TI6 = TI4*TI2

      SSA = TI
      SSB = TI3/48.D0
      SSC = 3.D0*TI2/16.D0
      SSD = 11.D0*TI/16.D0
      SSE = 1.D0

      Ra(1) = RXYZ(1,I)
      Ra(2) = RXYZ(2,I)
      Ra(3) = RXYZ(3,I)

      do newj = 1,totnebcoul(I)
        J = NEBCOUL(1, NEWJ, I)
        PBCI = NEBCOUL(2, NEWJ, I)
        PBCJ = NEBCOUL(3, NEWJ, I)
        PBCK = NEBCOUL(4, NEWJ, I)
        Rb(1) = RXYZ(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
             REAL(PBCK)*BOX(3,1)

        Rb(2) = RXYZ(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
             REAL(PBCK)*BOX(3,2)

        Rb(3) = RXYZ(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
             REAL(PBCK)*BOX(3,3)
        Rab = Rb-Ra  ! OBS b - a !!!
        dR = norm2(Rab)
        MAGR = dR
        MAGR2 = dR*dR

        if ((dR <= COULCUT).and.(dR > 1e-12)) then

          TJ = TFACT*U(Element_Pointer(J))
          DC = Rab/dR

          Z = abs(CALPHA*MAGR)
          NUMREP_ERFC = erfc(Z)

          CA = NUMREP_ERFC/MAGR
          E(I,J) = E(I,J) + CA
          ccnt = ccnt + 1
          CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
          EXPTI = exp(-TI*MAGR )

          if (Element_Pointer(I).eq.Element_Pointer(J)) then
            E(I,J) = E(I,J) - EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
            ccnt = ccnt + 1
          else
            TJ2 = TJ*TJ
            TJ3 = TJ2*TJ
            TJ4 = TJ2*TJ2
            TJ6 = TJ4*TJ2
            EXPTJ = exp( -TJ*MAGR )
            TI2MTJ2 = TI2 - TJ2
            TJ2MTI2 = -TI2MTJ2
            SA = TI
            SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
            SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
            SD = TJ
            SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
            SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

            E(I,J) = E(I,J) - ((EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
          endif

        endif
      enddo
    enddo
    !$OMP END PARALLEL DO

    E = KECONST*E

  end subroutine Ewald_Real_Space_Matrix_latte

  subroutine  Ewald_Real_Space_latte(COULOMBV,I,RXYZ,Box, &
       DELTAQ,U,Element_Pointer,Nr_atoms,COULACC,nebcoul,totnebcoul,HDIM,Max_Nr_Neigh,Nr_Elem)

    implicit none

    integer,    parameter     :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, HDIM, Max_Nr_Neigh, I, Nr_Elem
    real(PREC), intent(in)    :: COULACC
    real(PREC)                :: TFACT, RELPERM, KECONST
    real(PREC), intent(in)    :: RXYZ(3,Nr_atoms), Box(3,3), DELTAQ(Nr_atoms)
    real(PREC), intent(in)    :: U(Nr_elem)
    real(PREC)                :: COULCUT, COULCUT2
    integer, intent(in)       :: Element_Pointer(Nr_atoms)
    integer,    intent(in)    :: totnebcoul(Nr_atoms), nebcoul(4,Max_Nr_Neigh,Nr_atoms)
    real(PREC), intent(out)   :: COULOMBV
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: NUMREP_ERFC, CA, FORCE, EXPTI, EXPTJ
    real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
    real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2
    integer                   :: J,K, ccnt, newj, PBCI,PBCJ,PBCK

    COULVOL = Box(1,1)*Box(2,2)*Box(3,3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0
    COULCUT = 12.D0
    CALPHA = SQRTX/COULCUT
    COULCUT2 = COULCUT*COULCUT
    CALPHA2 = CALPHA*CALPHA

    RELPERM = ONE
    KECONST = 14.3996437701414D0*RELPERM
    TFACT  = 16.0D0/(5.0D0*KECONST)

    COULOMBV = ZERO

    TI = TFACT*U(Element_Pointer(I))
    TI2 = TI*TI
    TI3 = TI2*TI
    TI4 = TI2*TI2
    TI6 = TI4*TI2

    SSA = TI
    SSB = TI3/48.D0
    SSC = 3.D0*TI2/16.D0
    SSD = 11.D0*TI/16.D0
    SSE = 1.D0

    Ra(1) = RXYZ(1,I)
    Ra(2) = RXYZ(2,I)
    Ra(3) = RXYZ(3,I)

    do newj = 1,totnebcoul(I)
      J = NEBCOUL(1, NEWJ, I)
      PBCI = NEBCOUL(2, NEWJ, I)
      PBCJ = NEBCOUL(3, NEWJ, I)
      PBCK = NEBCOUL(4, NEWJ, I)
      Rb(1) = RXYZ(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
           REAL(PBCK)*BOX(3,1)

      Rb(2) = RXYZ(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
           REAL(PBCK)*BOX(3,2)

      Rb(3) = RXYZ(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
           REAL(PBCK)*BOX(3,3)
      Rab = Rb-Ra  ! OBS b - a !!!
      dR = norm2(Rab)
      MAGR = dR
      MAGR2 = dR*dR

      if ((dR <= COULCUT).and.(dR > 1e-12)) then

        TJ = TFACT*U(Element_Pointer(J))
        DC = Rab/dR

        Z = abs(CALPHA*MAGR)
        NUMREP_ERFC = erfc(Z)

        CA = NUMREP_ERFC/MAGR
        COULOMBV = COULOMBV + DELTAQ(J)*CA
        ccnt = ccnt + 1
        CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
        EXPTI = exp(-TI*MAGR )

        if (Element_Pointer(I).eq.Element_Pointer(J)) then
          COULOMBV = COULOMBV - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
          ccnt = ccnt + 1
        else
          TJ2 = TJ*TJ
          TJ3 = TJ2*TJ
          TJ4 = TJ2*TJ2
          TJ6 = TJ4*TJ2
          EXPTJ = exp( -TJ*MAGR )
          TI2MTJ2 = TI2 - TJ2
          TJ2MTI2 = -TI2MTJ2
          SA = TI
          SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
          SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
          SD = TJ
          SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
          SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

          COULOMBV = COULOMBV - (DELTAQ(J)*(EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
        endif

      endif
    enddo
    COULOMBV = KECONST*COULOMBV

  end subroutine Ewald_Real_Space_latte

  subroutine  Ewald_Real_Space_Test(COULOMBV,I,RX,RY,RZ,LBox, &
       DELTAQ,U,Element_Type,Nr_atoms,COULACC,nnRx,nnRy,nnRz,nrnnlist,nnType,Max_Nr_Neigh)

    implicit none

    integer,    parameter     :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, Max_Nr_Neigh, I
    real(PREC), intent(in)    :: COULACC
    real(PREC)                :: TFACT, RELPERM, KECONST
    real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3), DELTAQ(Nr_atoms)
    real(PREC), intent(in)    :: U(Nr_atoms)
    real(PREC)                :: COULCUT, COULCUT2
    character(10), intent(in) :: Element_Type(Nr_atoms)
    integer,    intent(in)    :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)    :: nnRx(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)    :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(out)   :: COULOMBV
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: NUMREP_ERFC, CA, FORCE, EXPTI, EXPTJ
    real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
    real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2
    integer                   :: J,K, ccnt, nnI

    COULVOL = LBox(1)*LBox(2)*LBox(3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0
    COULCUT = 12.D0
    CALPHA = SQRTX/COULCUT
    COULCUT2 = COULCUT*COULCUT
    CALPHA2 = CALPHA*CALPHA

    RELPERM = ONE
    KECONST = 14.3996437701414D0*RELPERM
    TFACT  = 16.0D0/(5.0D0*KECONST)

    COULOMBV = ZERO

    TI = TFACT*U(I)
    TI2 = TI*TI
    TI3 = TI2*TI
    TI4 = TI2*TI2
    TI6 = TI4*TI2

    SSA = TI
    SSB = TI3/48.D0
    SSC = 3.D0*TI2/16.D0
    SSD = 11.D0*TI/16.D0
    SSE = 1.D0

    Ra(1) = RX(I)
    Ra(2) = RY(I)
    Ra(3) = RZ(I)

    !    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nnI,J,Rb,Rab,dR,MAGR,MAGR2,TJ,DC,Z,NUMREP_ERFC,CA) &
    !    !$OMP REDUCTION(+:COULOMBV)
    do nnI = 1,nrnnlist(I)
      Rb(1) = nnRx(I,nnI)
      Rb(2) = nnRy(I,nnI)
      Rb(3) = nnRz(I,nnI)
      J = nnType(I,nnI)
      Rab = Rb-Ra  ! OBS b - a !!!
      dR = norm2(Rab)
      MAGR = dR
      MAGR2 = dR*dR

      if ((dR <= COULCUT).and.(dR > 1e-12)) then

        TJ = TFACT*U(J)
        DC = Rab/dR

        ! Not Using Numerical Recipes ERFC
        Z = abs(CALPHA*MAGR)
        NUMREP_ERFC = erfc(Z)

        CA = NUMREP_ERFC/MAGR
        COULOMBV = COULOMBV + DELTAQ(J)*CA
        ccnt = ccnt + 1
        !TEST(ccnt) = DELTAQ(J)*CA
        CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
        EXPTI = exp(-TI*MAGR )

        if (Element_Type(I).eq.Element_Type(J)) then
          COULOMBV = COULOMBV - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
          ccnt = ccnt + 1
          !TEST(ccnt) = - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
        else
          TJ2 = TJ*TJ
          TJ3 = TJ2*TJ
          TJ4 = TJ2*TJ2
          TJ6 = TJ4*TJ2
          EXPTJ = exp( -TJ*MAGR )
          TI2MTJ2 = TI2 - TJ2
          TJ2MTI2 = -TI2MTJ2
          SA = TI
          SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
          SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
          SD = TJ
          SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
          SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

          COULOMBV = COULOMBV - (DELTAQ(J)*(EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
        endif

      endif
    enddo
    !    !$OMP END PARALLEL DO
    COULOMBV = KECONST*COULOMBV

  end subroutine Ewald_Real_Space_Test

  subroutine  Ewald_Real_Space(COULOMBV,FCOUL,I,RX,RY,RZ,LBox, &
       DELTAQ,U,Element_Type,Nr_atoms,COULACC,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)

    implicit none

    integer,    parameter     :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, HDIM, Max_Nr_Neigh, I
    real(PREC), intent(in)    :: COULACC, TIMERATIO
    real(PREC)                :: TFACT, RELPERM, KECONST
    real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3), DELTAQ(Nr_atoms)
    real(PREC), intent(in)    :: U(Nr_atoms)
    real(PREC)                :: COULCUT, COULCUT2
    character(10), intent(in) :: Element_Type(Nr_atoms)
    integer,    intent(in)    :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)    :: nnRx(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(in)    :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
    real(PREC), intent(out)   :: COULOMBV, FCOUL(3)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: NUMREP_ERFC, CA, FORCE, EXPTI, EXPTJ
    real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
    real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2
    integer                   :: J,K, ccnt, nnI

    COULVOL = LBox(1)*LBox(2)*LBox(3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0
    COULCUT = 12.D0
    CALPHA = SQRTX/COULCUT
    COULCUT2 = COULCUT*COULCUT
    CALPHA2 = CALPHA*CALPHA

    RELPERM = ONE
    KECONST = 14.3996437701414D0*RELPERM
    TFACT  = 16.0D0/(5.0D0*KECONST)

    FCOUL = ZERO
    COULOMBV = ZERO

    TI = TFACT*U(I)
    TI2 = TI*TI
    TI3 = TI2*TI
    TI4 = TI2*TI2
    TI6 = TI4*TI2

    SSA = TI
    SSB = TI3/48.D0
    SSC = 3.D0*TI2/16.D0
    SSD = 11.D0*TI/16.D0
    SSE = 1.D0

    Ra(1) = RX(I)
    Ra(2) = RY(I)
    Ra(3) = RZ(I)

    do nnI = 1,nrnnlist(I)
      Rb(1) = nnRx(I,nnI)
      Rb(2) = nnRy(I,nnI)
      Rb(3) = nnRz(I,nnI)
      J = nnType(I,nnI)
      Rab = Rb-Ra  ! OBS b - a !!!
      dR = norm2(Rab)
      MAGR = dR
      MAGR2 = dR*dR

      if ((dR <= COULCUT).and.(dR > 1e-12)) then

        TJ = TFACT*U(J)
        DC = Rab/dR

        ! Not Using Numerical Recipes ERFC
        Z = abs(CALPHA*MAGR)
        NUMREP_ERFC = erfc(Z)

        CA = NUMREP_ERFC/MAGR
        COULOMBV = COULOMBV + DELTAQ(J)*CA
        ccnt = ccnt + 1
        !TEST(ccnt) = DELTAQ(J)*CA
        CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
        FORCE = -KECONST*DELTAQ(I)*DELTAQ(J)*CA/MAGR
        EXPTI = exp(-TI*MAGR )

        if (Element_Type(I).eq.Element_Type(J)) then
          COULOMBV = COULOMBV - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
          ccnt = ccnt + 1
          !TEST(ccnt) = - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
          FORCE = FORCE + (KECONST*DELTAQ(I)*DELTAQ(J)*EXPTI)*((SSE/MAGR2 - TWO*SSB*MAGR - SSC) &
               + SSA*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR))
        else
          TJ2 = TJ*TJ
          TJ3 = TJ2*TJ
          TJ4 = TJ2*TJ2
          TJ6 = TJ4*TJ2
          EXPTJ = exp( -TJ*MAGR )
          TI2MTJ2 = TI2 - TJ2
          TJ2MTI2 = -TI2MTJ2
          SA = TI
          SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
          SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
          SD = TJ
          SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
          SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

          COULOMBV = COULOMBV - (DELTAQ(J)*(EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
          FORCE = FORCE + KECONST*DELTAQ(I)*DELTAQ(J)*((EXPTI*(SA*(SB - (SC/MAGR)) - (SC/MAGR2))) &
               + (EXPTJ*(SD*(SE - (SF/MAGR)) - (SF/MAGR2))))
        endif

        FCOUL(1) = FCOUL(1) + DC(1)*FORCE
        FCOUL(2) = FCOUL(2) + DC(2)*FORCE
        FCOUL(3) = FCOUL(3) + DC(3)*FORCE
      endif
    enddo
    COULOMBV = KECONST*COULOMBV

  end subroutine Ewald_Real_Space

  subroutine Ewald_k_Space_latte(COULOMBV,RXYZ,Box,DELTAQ,Nr_atoms,COULACC,Max_Nr_Neigh)

    implicit none

    integer, parameter        :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0, FOUR = 4.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0, EIGHT = 8.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, Max_Nr_Neigh
    real(PREC), intent(in)    :: COULACC
    real(PREC)                :: KECONST, TFACT, RELPERM
    real(PREC), intent(in)    :: RXYZ(3,Nr_atoms), Box(3,3), DELTAQ(Nr_atoms)
    real(PREC)                :: COULCUT, COULCUT2
    real(PREC), intent(out)   :: COULOMBV(Nr_atoms)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: CORRFACT,FOURCALPHA2, FORCE
    real(PREC)                :: RECIPVECS(3,3),SINLIST(Nr_atoms),COSLIST(Nr_Atoms)
    real(PREC)                :: K(3),L11,L12,L13,M21,M22,M23,K2,KCUTOFF,KCUTOFF2,PREFACTOR
    real(PREC)                :: PREVIR, COSSUM,SINSUM,DOT, KEPREF, COSSUM2, SINSUM2

    integer                   :: I,J,L,M,N, ccnt, nnI, LMAX,MMAX,NMAX,NMIN,MMIN

    COULVOL = Box(1,1)*Box(2,2)*Box(3,3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0

    COULCUT = 12.0D0
    CALPHA = SQRTX/COULCUT

    COULCUT2 = COULCUT*COULCUT
    KCUTOFF = TWO*CALPHA*SQRTX
    KCUTOFF2 = KCUTOFF*KCUTOFF
    CALPHA2 = CALPHA*CALPHA
    FOURCALPHA2 = FOUR*CALPHA2

    RECIPVECS = ZERO
    RECIPVECS(1,1) = TWO*pi/Box(1,1)
    RECIPVECS(2,2) = TWO*pi/Box(2,2)
    RECIPVECS(3,3) = TWO*pi/Box(3,3)
    LMAX = floor(KCUTOFF / sqrt(RECIPVECS(1,1)*RECIPVECS(1,1)))
    MMAX = floor(KCUTOFF / sqrt(RECIPVECS(2,2)*RECIPVECS(2,2)))
    NMAX = floor(KCUTOFF / sqrt(RECIPVECS(3,3)*RECIPVECS(3,3)))

    RELPERM = 1.D0
    KECONST = 14.3996437701414D0*RELPERM

    COULOMBV = ZERO
    SINLIST = ZERO
    COSLIST = ZERO

    do L = 0,LMAX

      if (L.eq.0) then
        MMIN = 0
      else
        MMIN = -MMAX
      endif

      L11 = L*RECIPVECS(1,1)
      L12 = L*RECIPVECS(1,2)
      L13 = L*RECIPVECS(1,3)

      do M = MMIN,MMAX

        NMIN = -NMAX
        if ((L==0).and.(M==0)) then
          NMIN = 1
        endif

        M21 = L11 + M*RECIPVECS(2,1)
        M22 = L12 + M*RECIPVECS(2,2)
        M23 = L13 + M*RECIPVECS(2,3)

        do N = NMIN,NMAX
          K(1) = M21 + N*RECIPVECS(3,1)
          K(2) = M22 + N*RECIPVECS(3,2)
          K(3) = M23 + N*RECIPVECS(3,3)
          K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)
          if (K2.le.KCUTOFF2) then
            PREFACTOR = EIGHT*pi*exp(-K2/(4.D0*CALPHA2))/(COULVOL*K2)
            PREVIR = (2.D0/K2) + (2.D0/(4.D0*CALPHA2));

            COSSUM = 0.D0
            SINSUM = 0.D0

            ! Doing the sin and cos sums
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,DOT) &
            !$OMP REDUCTION(+:COSSUM) &
            !$OMP REDUCTION(+:SINSUM)
            do I = 1,Nr_atoms
              DOT = K(1)*RXYZ(1,I) + K(2)*RXYZ(2,I) + K(3)*RXYZ(3,I)
              ! We re-use these in the next loop...
              SINLIST(I) = sin(DOT)
              COSLIST(I) = cos(DOT)
              COSSUM = COSSUM + DELTAQ(I)*COSLIST(I)
              SINSUM = SINSUM + DELTAQ(I)*SINLIST(I)
            enddo
            !$OMP END PARALLEL DO
            COSSUM2 = COSSUM*COSSUM
            SINSUM2 = SINSUM*SINSUM

            ! Add up energy and force contributions

            KEPREF = KECONST*PREFACTOR
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
            do I = 1,Nr_atoms
              COULOMBV(I) = COULOMBV(I) + KEPREF*(COSLIST(I)*COSSUM + SINLIST(I)*SINSUM)
            enddo
            !$OMP END PARALLEL DO

            KEPREF = KEPREF*(COSSUM2 + SINSUM2)
          endif
        enddo
      enddo
    enddo

    ! Point self energy
    CORRFACT = 2.D0*KECONST*CALPHA/SQRTPI;
    COULOMBV = COULOMBV - CORRFACT*DELTAQ;

  end subroutine Ewald_k_Space_latte

  subroutine Ewald_k_Space_Matrix_latte(E,RXYZ,Box,Nr_atoms,COULACC,Max_Nr_Neigh,nebcoul,totnebcoul)

    implicit none

    integer, parameter        :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0, FOUR = 4.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0, EIGHT = 8.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, Max_Nr_Neigh,nebcoul(4,Max_Nr_Neigh,Nr_atoms),totnebcoul(Nr_atoms)
    real(PREC), intent(in)    :: COULACC
    real(PREC)                :: KECONST, TFACT, RELPERM
    real(PREC), intent(in)    :: RXYZ(3,Nr_atoms), Box(3,3)
    real(PREC)                :: COULCUT, COULCUT2
    real(PREC), intent(out)   :: E(Nr_atoms,Nr_atoms)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: CORRFACT,FOURCALPHA2, FORCE
    real(PREC)                :: RECIPVECS(3,3),SINLIST(Nr_atoms),COSLIST(Nr_Atoms)
    real(PREC)                :: K(3),L11,L12,L13,M21,M22,M23,K2,KCUTOFF,KCUTOFF2,PREFACTOR
    real(PREC)                :: PREVIR, COSSUM,SINSUM,DOT, KEPREF, COSSUM2, SINSUM2

    integer                   :: I,J,L,M,N, newj, ccnt, nnI, LMAX,MMAX,NMAX,NMIN,MMIN

    COULVOL = Box(1,1)*Box(2,2)*Box(3,3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0

    COULCUT = 12.0D0
    CALPHA = SQRTX/COULCUT

    COULCUT2 = COULCUT*COULCUT
    KCUTOFF = TWO*CALPHA*SQRTX
    KCUTOFF2 = KCUTOFF*KCUTOFF
    CALPHA2 = CALPHA*CALPHA
    FOURCALPHA2 = FOUR*CALPHA2

    RECIPVECS = ZERO
    RECIPVECS(1,1) = TWO*pi/Box(1,1)
    RECIPVECS(2,2) = TWO*pi/Box(2,2)
    RECIPVECS(3,3) = TWO*pi/Box(3,3)
    LMAX = floor(KCUTOFF / sqrt(RECIPVECS(1,1)*RECIPVECS(1,1)))
    MMAX = floor(KCUTOFF / sqrt(RECIPVECS(2,2)*RECIPVECS(2,2)))
    NMAX = floor(KCUTOFF / sqrt(RECIPVECS(3,3)*RECIPVECS(3,3)))

    RELPERM = 1.D0
    KECONST = 14.3996437701414D0*RELPERM

    !COULOMBV = ZERO
    SINLIST = ZERO
    COSLIST = ZERO

    E = 0.0

    do L = 0,LMAX

      if (L.eq.0) then
        MMIN = 0
      else
        MMIN = -MMAX
      endif

      L11 = L*RECIPVECS(1,1)
      L12 = L*RECIPVECS(1,2)
      L13 = L*RECIPVECS(1,3)

      do M = MMIN,MMAX

        NMIN = -NMAX
        if ((L==0).and.(M==0)) then
          NMIN = 1
        endif

        M21 = L11 + M*RECIPVECS(2,1)
        M22 = L12 + M*RECIPVECS(2,2)
        M23 = L13 + M*RECIPVECS(2,3)

        do N = NMIN,NMAX
          K(1) = M21 + N*RECIPVECS(3,1)
          K(2) = M22 + N*RECIPVECS(3,2)
          K(3) = M23 + N*RECIPVECS(3,3)
          K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)
          if (K2.le.KCUTOFF2) then
            PREFACTOR = EIGHT*pi*exp(-K2/(4.D0*CALPHA2))/(COULVOL*K2)
            PREVIR = (2.D0/K2) + (2.D0/(4.D0*CALPHA2));


            ! Doing the sin and cos sums
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,DOT)
            do I = 1,Nr_atoms
              DOT = K(1)*RXYZ(1,I) + K(2)*RXYZ(2,I) + K(3)*RXYZ(3,I)
              ! We re-use these in the next loop...
              SINLIST(I) = sin(DOT)
              COSLIST(I) = cos(DOT)
              !   COSSUM = COSSUM + DELTAQ(I)*COSLIST(I)
              !   SINSUM = SINSUM + DELTAQ(I)*SINLIST(I)
            enddo
            !$OMP END PARALLEL DO

            ! Add up energy and force contributions

            KEPREF = KECONST*PREFACTOR
            CORRFACT = 2.D0*KECONST*CALPHA/SQRTPI;
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,NEWJ)
            do I = 1,Nr_atoms
              do newj = 1,totnebcoul(I)
                J = NEBCOUL(1, NEWJ, I)
                E(I,J) = E(I,J) + KEPREF*(COSLIST(I)*COSLIST(J)+SINLIST(I)*SINLIST(J))
                !    COULOMBV(I) = COULOMBV(I) + KEPREF*(COSLIST(I)*COSSUM + SINLIST(I)*SINSUM)
              enddo
              E(I,I) = E(I,I) - CORRFACT
            enddo
            !$OMP END PARALLEL DO

            !KEPREF = KEPREF*(COSSUM2 + SINSUM2)
          endif
        enddo
      enddo
    enddo

    ! Point self energy
    !CORRFACT = 2.D0*KECONST*CALPHA/SQRTPI;
    !COULOMBV = COULOMBV - CORRFACT*DELTAQ;

  end subroutine Ewald_k_Space_Matrix_latte

  subroutine Ewald_k_Space_latte_single(COULOMBV,J,RXYZ,Box,DELTAQ,Nr_atoms,COULACC)

    implicit none

    integer, parameter        :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0, FOUR = 4.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0, EIGHT = 8.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, J
    real(PREC), intent(in)    :: COULACC
    real(PREC)                :: KECONST, TFACT, RELPERM
    real(PREC), intent(in)    :: RXYZ(3,Nr_atoms), Box(3,3), DELTAQ(Nr_atoms)
    real(PREC)                :: COULCUT, COULCUT2
    real(PREC), intent(out)   :: COULOMBV(Nr_atoms)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2
    real(PREC)                :: CORRFACT,FOURCALPHA2
    real(PREC)                :: RECIPVECS(3,3)
    real(PREC)                :: K(3),L11,L12,L13,M21,M22,M23,K2,KCUTOFF,KCUTOFF2,PREFACTOR
    real(PREC)                :: IDOT, JDOT, COSJDOT, SINJDOT, KEPREF

    integer                   :: I,L,M,N, ccnt, nnI, LMAX,MMAX,NMAX,NMIN,MMIN

    COULVOL = Box(1,1)*Box(2,2)*Box(3,3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0

    COULCUT = 12.0D0
    CALPHA = SQRTX/COULCUT

    COULCUT2 = COULCUT*COULCUT
    KCUTOFF = TWO*CALPHA*SQRTX
    KCUTOFF2 = KCUTOFF*KCUTOFF
    CALPHA2 = CALPHA*CALPHA
    FOURCALPHA2 = FOUR*CALPHA2

    RECIPVECS = ZERO
    RECIPVECS(1,1) = TWO*pi/Box(1,1)
    RECIPVECS(2,2) = TWO*pi/Box(2,2)
    RECIPVECS(3,3) = TWO*pi/Box(3,3)
    LMAX = floor(KCUTOFF / sqrt(RECIPVECS(1,1)*RECIPVECS(1,1)))
    MMAX = floor(KCUTOFF / sqrt(RECIPVECS(2,2)*RECIPVECS(2,2)))
    NMAX = floor(KCUTOFF / sqrt(RECIPVECS(3,3)*RECIPVECS(3,3)))

    RELPERM = 1.D0
    KECONST = 14.3996437701414D0*RELPERM

    COULOMBV = ZERO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,IDOT,JDOT,COSJDOT,SINJDOT,L,M,N,MMIN,MMAX,NMIN,NMAX,L11,M22,K,K2,PREFACTOR,KEPREF)
    do L = 0,LMAX

      if (L.eq.0) then
        MMIN = 0
      else
        MMIN = -MMAX
      endif

      L11 = L*RECIPVECS(1,1)

      do M = MMIN,MMAX

        NMIN = -NMAX
        if ((L==0).and.(M==0)) then
          NMIN = 1
        endif

        M22 = M*RECIPVECS(2,2)

        do N = NMIN,NMAX
          K(1) = L11
          K(2) = M22
          K(3) = N*RECIPVECS(3,3)
          K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)

          PREFACTOR = EIGHT*pi*exp(-K2/(4.D0*CALPHA2))/(COULVOL*K2)
          KEPREF = KECONST*PREFACTOR
          JDOT = K(1)*RXYZ(1,J) + K(2)*RXYZ(2,J) + K(3)*RXYZ(3,J)
          SINJDOT = sin(JDOT)
          COSJDOT = cos(JDOT)
          do I = 1,Nr_atoms
            IDOT = K(1)*RXYZ(1,I) + K(2)*RXYZ(2,I) + K(3)*RXYZ(3,I)
            COULOMBV(I) = COULOMBV(I) + KEPREF*DELTAQ(J)*(COSJDOT*cos(IDOT)+SINJDOT*sin(IDOT))
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! Point self energy
    CORRFACT = 2.D0*KECONST*CALPHA/SQRTPI;
    COULOMBV = COULOMBV - CORRFACT*DELTAQ;

  end subroutine Ewald_k_Space_latte_single

  subroutine Ewald_k_Space_Test(COULOMBV,RX,RY,RZ,LBox,DELTAQ,Nr_atoms,COULACC,Max_Nr_Neigh)
    !
    implicit none
    !
    integer, parameter        :: PREC = 8
    real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0, FOUR = 4.D0
    real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0, EIGHT = 8.D0
    real(PREC), parameter     :: pi = 3.14159265358979323846264D0
    real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
    integer,    intent(in)    :: Nr_atoms, Max_Nr_Neigh
    real(PREC), intent(in)    :: COULACC
    real(PREC)                :: KECONST, TFACT, RELPERM
    real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3), DELTAQ(Nr_atoms)
    real(PREC)                :: COULCUT, COULCUT2
    real(PREC), intent(out)   :: COULOMBV(Nr_atoms)
    real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
    real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
    real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
    real(PREC)                :: CORRFACT,FOURCALPHA2, FORCE
    real(PREC)                :: RECIPVECS(3,3),SINLIST(Nr_atoms),COSLIST(Nr_Atoms)
    real(PREC)                :: K(3),L11,L12,L13,M21,M22,M23,K2,KCUTOFF,KCUTOFF2,PREFACTOR
    real(PREC)                :: PREVIR, COSSUM,SINSUM,DOT, KEPREF, COSSUM2, SINSUM2

    integer                   :: I,J,L,M,N, ccnt, nnI, LMAX,MMAX,NMAX,NMIN,MMIN
    !
    COULVOL = LBox(1)*LBox(2)*LBox(3)
    SQRTX = sqrt(-log(COULACC))

    ccnt = 0

    COULCUT = 12.0D0
    CALPHA = SQRTX/COULCUT

    COULCUT2 = COULCUT*COULCUT
    KCUTOFF = TWO*CALPHA*SQRTX
    KCUTOFF2 = KCUTOFF*KCUTOFF
    CALPHA2 = CALPHA*CALPHA
    FOURCALPHA2 = FOUR*CALPHA2

    RECIPVECS = ZERO
    RECIPVECS(1,1) = TWO*pi/LBox(1)
    RECIPVECS(2,2) = TWO*pi/LBox(2)
    RECIPVECS(3,3) = TWO*pi/LBox(3)
    LMAX = floor(KCUTOFF / sqrt(RECIPVECS(1,1)*RECIPVECS(1,1)))
    MMAX = floor(KCUTOFF / sqrt(RECIPVECS(2,2)*RECIPVECS(2,2)))
    NMAX = floor(KCUTOFF / sqrt(RECIPVECS(3,3)*RECIPVECS(3,3)))

    RELPERM = 1.D0
    KECONST = 14.3996437701414D0*RELPERM

    COULOMBV = ZERO
    SINLIST = ZERO
    COSLIST = ZERO

    do L = 0,LMAX

      if (L.eq.0) then
        MMIN = 0
      else
        MMIN = -MMAX
      endif

      L11 = L*RECIPVECS(1,1)
      L12 = L*RECIPVECS(1,2)
      L13 = L*RECIPVECS(1,3)

      do M = MMIN,MMAX

        NMIN = -NMAX
        if ((L==0).and.(M==0)) then
          NMIN = 1
        endif

        M21 = L11 + M*RECIPVECS(2,1)
        M22 = L12 + M*RECIPVECS(2,2)
        M23 = L13 + M*RECIPVECS(2,3)

        do N = NMIN,NMAX
          K(1) = M21 + N*RECIPVECS(3,1)
          K(2) = M22 + N*RECIPVECS(3,2)
          K(3) = M23 + N*RECIPVECS(3,3)
          K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)
          if (K2.le.KCUTOFF2) then
            PREFACTOR = EIGHT*pi*exp(-K2/(4.D0*CALPHA2))/(COULVOL*K2)
            PREVIR = (2.D0/K2) + (2.D0/(4.D0*CALPHA2));

            COSSUM = 0.D0
            SINSUM = 0.D0

            ! Doing the sin and cos sums
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,DOT) &
            !$OMP REDUCTION(+:COSSUM) &
            !$OMP REDUCTION(+:SINSUM)
            do I = 1,Nr_atoms
              DOT = K(1)*RX(I) + K(2)*RY(I) + K(3)*RZ(I)
              ! We re-use these in the next loop...
              SINLIST(I) = sin(DOT)
              COSLIST(I) = cos(DOT)
              COSSUM = COSSUM + DELTAQ(I)*COSLIST(I)
              SINSUM = SINSUM + DELTAQ(I)*SINLIST(I)
            enddo
            !$OMP END PARALLEL DO
            COSSUM2 = COSSUM*COSSUM
            SINSUM2 = SINSUM*SINSUM

            ! Add up energy and force contributions
            !
            KEPREF = KECONST*PREFACTOR
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
            do I = 1,Nr_atoms
              COULOMBV(I) = COULOMBV(I) + KEPREF*(COSLIST(I)*COSSUM + SINLIST(I)*SINSUM)
            enddo
            !$OMP END PARALLEL DO

            KEPREF = KEPREF*(COSSUM2 + SINSUM2)
          endif
        enddo
      enddo
    enddo

    ! Point self energy
    CORRFACT = 2.D0*KECONST*CALPHA/SQRTPI;
    COULOMBV = COULOMBV - CORRFACT*DELTAQ;

  end subroutine Ewald_k_Space_Test

end module prg_ewald_mod
