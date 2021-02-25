program algo 

integer, allocatable :: caca(:)

allocate(caca(3))
caca = 0
write(*,*)"caca",caca
call gpmdcov_algo(caca)

contains
  subroutine gpmdcov_algo(cacaint)

   implicit none
   integer, allocatable :: cacaint(:)

  cacaint = 0
   write(*,*)"caca", cacaint
 end subroutine gpmdcov_algo
end
