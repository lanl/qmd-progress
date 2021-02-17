
  !> To write output file or perform some analysis
  !!
  subroutine gpmdcov_writeout()
    use gpmdcov_vars

    if(allocated(sy%resindex))deallocate(sy%resindex)
    allocate(sy%resindex(sy%nats))
    sy%resindex(sy%nats)=-100

#ifdef DO_MPI
    !      !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
      !       do ipt = 1,gpat%TotalParts
#endif
      write(filename,*)ipt
      auxchar = adjustl(trim(filename))
      filename = "part_"//auxchar
      call prg_write_system(syprt(ipt),filename,"pdb")
      do j=1,gpat%sgraph(ipt)%llsize
        sy%resindex(gpat%sgraph(ipt)%core_halo_index(j)+1) = ipt
      enddo
    enddo
    call prg_write_system(sy,"system_parts","pdb")

    !> Writting the extension of the graph as a resindex
    !     if(.not.allocated(row1))allocate(row1(sy%nats))
    !     call bml_get_row(g_bml,1,row1)
    !     write(*,*)"N attached to molecule i =", sum(row1(:))
    !
    !     row1 = 0
    !     do i=1,nl%nrnnstruct(20000)
    !       row1(nl%nnStruct(20000,i)) = 1
    !     enddo
    !
    !     sy%resindex = row1

    !     if(.not.allocated(auxcharge1))allocate(auxcharge1(sy%nats))
    !     if(.not.allocated(sy%net_charge))allocate(sy%net_charge(sy%nats))

    !     auxcharge1=sy%net_charge
    !     sy%net_charge = 0
    !     do j=1,gpat%sgraph(1)%llsize
    !       sy%net_charge(gpat%sgraph(1)%core_halo_index(j)+1) = 1.0_dp
    !     enddo
    !     do j=gpat%sgraph(1)%llsize+1,gpat%sgraph(1)%lsize
    !       sy%net_charge(gpat%sgraph(1)%core_halo_index(j)+1) = -1.0_dp
    !     enddo

    !     call prg_write_system(sy,"connections_of_atom_i","pdb")
    !     call prg_write_trajectory(sy,mdstep,1,0.1_dp,"traj_parts","xyz")
    !     sy%net_charge = auxcharge1

    deallocate(sy%resindex)
    !     deallocate(row1)

  end subroutine gpmdcov_writeout

