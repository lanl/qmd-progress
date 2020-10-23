program main

  use prg_sp2_tensorcore_mod
  
  real(4), allocatable :: H(:)
  real(8), allocatable :: D(:)
  real(4) :: eps, bndfil, idemtol
  integer :: minsp2iter, maxsp2iter, verbose
  character(3) :: sp2conv
  
  eps = 1.0e-16
  N = 100 
  Nocc = 1 
  bndfil = real(Nocc,8)/real(N,8)
  sp2conv = "rel"
  idemtol = 1.0e-16 
  verbose = 1

  allocate(H(N*N))

  !Build H 
  write(*,*)"Constructing H ..."
  do i=1,N
    do j=1,N
      H(i*N + j) = exp(-0.5*abs(i-j)*sin(real((i+1),8)))
      H(j*N + i) = H(i*N + j)
    enddo
  enddo
  allocate(D(N*N))

! call algo()
  write(*,*)"Entering prg_sp2_tensorcore ..."
  call prg_sp2_tensorcore(N,H,D,eps,bndfil,minsp2iter,&
       & maxsp2iter,sp2conv,idemtol,verbose)

end 
