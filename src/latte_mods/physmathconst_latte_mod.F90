!> Physical and Mathematical constants.
!! Here we include some unit transformations that might be 
!! used in any physics/chemistry code develompmet.
!! \ingroup LATTE
!!
module PhysMathConst_latte_mod

  implicit none

  private 

  integer, parameter :: dp = kind(1.0d0)

  !> Physical and Mathematical constants.
  !!
  real(dp),public, parameter :: kb = 8.61734215D-5  ![eV/K]
  real(dp),public, parameter :: pi = 3.14159265358979323846264338327950_dp
  real(dp),public, parameter :: hbar = 0.65821188926_dp ![eV fs]
  real(dp),public, parameter :: e0 = 5.526349954E-3_dp ![eV^-1 Ang^-1] 
  real(dp),public, parameter :: Ang2overeVtoAU = 96.54 
  real(dp),public, parameter :: Ang2overeVtoAng3 = 1.0_dp/(4.0_dp*pi*e0)
   real(dp),public, parameter :: AUtoAng3 = 0.529_dp**3
  
  !> Unit transformations.
  !!
  real(dp),public, parameter :: mA2efs = 160.21_dp  
  real(dp),public, parameter :: kcalmol2ev = 1.0_dp/23.06_dp
  real(dp),public, parameter :: hartree2ev = 27.2116_dp
  real(dp),public, parameter :: G0 = 77473.39_dp
  real(dp),public, parameter :: ev2cmm1=8065.73
  real(dp),public, parameter :: eA2Debye=4.80245_dp

end module PhysMathConst_latte_mod
