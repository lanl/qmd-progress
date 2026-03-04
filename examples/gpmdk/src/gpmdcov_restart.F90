  !> Dump GPMD
  subroutine gpmdcov_restart()
    use gpmdcov_vars
    use prg_parallel_mod

    integer :: array_size(1), nats(1)
    
    if(myRank.eq.1)then
       write(*,*)"Restarting ..."
       
       open(1,file='restart.dmp',form="unformatted",access="sequential",status="old")
       read(1)nats(1)
    endif
    
    call prg_BcastIntParallel(nats,1,1)
    sy%nats = nats(1)
    
    if(.not.allocated(sy%symbol))allocate(sy%symbol(sy%nats))
    if(.not.allocated(sy%atomic_number))allocate(sy%atomic_number(sy%nats))
    if(.not.allocated(sy%coordinate))allocate(sy%coordinate(3,sy%nats))
    if(.not.allocated(sy%velocity))allocate(sy%velocity(3,sy%nats))
    if(.not.allocated(sy%force))allocate(sy%force(3,sy%nats))
    if(.not.allocated(sy%net_charge))allocate(sy%net_charge(sy%nats))
    if(.not.allocated(sy%mass))allocate(sy%mass(sy%nats))
    if(.not.allocated(sy%spindex))allocate(sy%spindex(sy%nsp))
    if(.not.allocated(sy%lattice_vector))allocate(sy%lattice_vector(3,3))
    if(.not.allocated(sy%spatnum))allocate(sy%spatnum(sy%nsp))
    if(.not.allocated(sy%spmass))allocate(sy%spmass(sy%nsp))
    if(.not.allocated(n))allocate(n(sy%nats))
    if(.not.allocated(n_0))allocate(n_0(sy%nats))
    if(.not.allocated(n_1))allocate(n_1(sy%nats))
    if(.not.allocated(n_2))allocate(n_2(sy%nats))
    if(.not.allocated(n_3))allocate(n_3(sy%nats))
    if(.not.allocated(n_4))allocate(n_4(sy%nats))
    if(.not.allocated(n_5))allocate(n_5(sy%nats))

    if(myRank.eq.1)then
       read(1)sy%symbol
       read(1)sy%atomic_number
       read(1)sy%coordinate
       read(1)sy%velocity
       
       read(1)sy%force
       read(1)sy%net_charge
       read(1)sy%mass
       read(1)sy%spindex
       read(1)sy%lattice_vector
       read(1)sy%spatnum
       read(1)sy%spmass

       
       read(1)mdstep
       read(1)n
       read(1)n_0
       read(1)n_1
       read(1)n_2
       read(1)n_3
       read(1)n_4
       read(1)n_5
       close(1)

    endif

#define BCAST_CHAR(x) \
    if(myRank.eq.1)then; \
       array_size(1) = size(x); \
    endif; \
    call prg_BcastIntParallel(array_size,1,1); \
    call prg_BcastParallel(x,array_size(1),1)

#define BCAST_INT(x) \
    if(myRank.eq.1)then; \
       array_size(1) = size(x); \
    endif; \
    call prg_BcastIntParallel(array_size,1,1); \
    call prg_BcastIntParallel(x,array_size(1),1)

#define BCAST_REAL(x) \
    if(myRank.eq.1)then; \
       array_size(1) = size(x); \
    endif; \
    call prg_BcastIntParallel(array_size,1,1); \
    call prg_BcastRealParallel(x,array_size(1),1)

    BCAST_CHAR(sy%symbol)
    BCAST_INT(sy%atomic_number)
    BCAST_REAL(sy%coordinate)
    BCAST_REAL(sy%velocity)
    BCAST_REAL(sy%force)
    BCAST_REAL(sy%net_charge)
    BCAST_REAL(sy%mass)
    BCAST_INT(sy%spindex)
    BCAST_REAL(sy%lattice_vector)
    BCAST_INT(sy%spatnum)
    BCAST_REAL(sy%spmass)
    BCAST_REAL(n)
    BCAST_REAL(n_0)
    BCAST_REAL(n_1)
    BCAST_REAL(n_2)
    BCAST_REAL(n_3)
    BCAST_REAL(n_4)
    BCAST_REAL(n_5)
    
    mdstep = 0
           
  end subroutine gpmdcov_restart

