program main

  use prg_sp2_tensorcore_mod
  
  real(4), allocatable :: H(:)
  real(8), allocatable :: D(:)
  real(4) :: eps, bndfil, idemtol
  integer :: minsp2iter, maxsp2iter, verbose
  character(3) :: sp2conv
  
  eps = 1.0e-16
  N = 10
  Nocc = 1 
  bndfil = real(Nocc,8)/real(N,8)
  sp2conv = "rel"
  idemtol = 1.0e-16 
  verbose = 1

  allocate(H(N*N))

  !Build H 
  write(*,*)"Constructing H ..."
  do i=1,N
    do j=i,N
      H((i-1)*N + j) = exp(-0.5*abs(real(i-j,4)))*sin(real(i,4))
      H((j-1)*N + i) = H((i-1)*N + j)
    enddo
  write(*,*)(H((i-1)*N + j),j=1,N)
  enddo
  allocate(D(N*N))

  D = 0.0D0
  
  write(*,*)"Entering prg_sp2_tensorcore ..."
  call prg_sp2_tensorcore_f(N,H,D,eps,bndfil,minsp2iter,&
       & maxsp2iter,sp2conv,idemtol,verbose)

  write(*,*)"Writing Density Matrix ..."
!  do i=1,N
!    write(*,*)(D((i-1)*N + j),j=1,N)
!  enddo
  


end 
