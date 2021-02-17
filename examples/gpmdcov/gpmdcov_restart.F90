  !> Dump GPMD
  subroutine gpmdcov_restart()
    use gpmdcov_vars

    write(*,*)"Restarting ..."

    open(1,file='restart.dmp',form="unformatted",access="sequential",status="old")
    read(1)sy%nats

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

    if(.not.allocated(n))allocate(n(sy%nats))
    if(.not.allocated(n_0))allocate(n_0(sy%nats))
    if(.not.allocated(n_1))allocate(n_1(sy%nats))
    if(.not.allocated(n_2))allocate(n_2(sy%nats))
    if(.not.allocated(n_3))allocate(n_3(sy%nats))
    if(.not.allocated(n_4))allocate(n_4(sy%nats))
    if(.not.allocated(n_5))allocate(n_5(sy%nats))

    read(1)mdstep
    read(1)n
    read(1)n_0
    read(1)n_1
    read(1)n_2
    read(1)n_3
    read(1)n_4
    read(1)n_5

    mdstep = 0

    close(1)

  end subroutine gpmdcov_restart

