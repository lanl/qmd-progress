
  !> Dump GPMD
  subroutine gpmdcov_dump()
    
    use gpmdcov_vars

    if(myRank == 1)then
      write(*,*)"Dumping resetart file ..."

      open(1,file='restart.dmp',form="unformatted",access="sequential")
      write(1)sy%nats
      write(1)sy%symbol
      write(1)sy%atomic_number
      write(1)sy%coordinate
      write(1)sy%velocity
      write(1)sy%force
      write(1)sy%net_charge
      write(1)sy%mass
      write(1)sy%spindex
      write(1)sy%lattice_vector
      write(1)sy%spatnum
      write(1)sy%spmass

      !Dump xlbo
      write(1)mdstep
      write(1)n
      write(1)n_0
      write(1)n_1
      write(1)n_2
      write(1)n_3
      write(1)n_4
      write(1)n_5

      close(1)
    endif

  end subroutine gpmdcov_dump

