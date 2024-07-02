module gpmdcov_reshuffle_mod

  contains 

  !> Reshuffle the parts
  !! partsInEachRank(getNRanks()) stores the number of partitions assigned to rank i
  !! reshuffle(j,i) assigns partition reshuffle(j,i) to rank i for j=1,partsInEachRank(getNRanks())
  !!
  subroutine gpmdcov_reshuffle()
    use gpmdcov_vars
    use gpmdcov_writeout_mod

    integer :: maxnparts, np

    maxnparts = 0
    do i=1,numRanks
      np = gpat%localPartMax(i)-gpat%localPartMin(i)+1
      maxnparts = max(maxnparts,np)
    enddo

    if(allocated(reshuffle))then
      deallocate(reshuffle)
      deallocate(partsInEachRank)
    endif

    allocate(reshuffle(maxnparts,numRanks))
    allocate(partsInEachRank(numRanks))

    reshuffle = 0
    icount = 0
    partsInEachRank = 0

    do j=1,maxnparts
      do i=1,getNRanks()
        np = gpat%localPartMax(i)-gpat%localPartMin(i)+1
        if(np > partsInEachRank(i))then
          icount = icount + 1
          partsInEachRank(i) = partsInEachRank(i) + 1
          reshuffle(partsInEachRank(i),i) = icount
          if(icount == nparts) exit
        endif
      enddo
      if(icount == nparts) exit
      do i=getNRanks(),1,-1
        np = gpat%localPartMax(i)-gpat%localPartMin(i)+1
        if(np > partsInEachRank(i))then
          icount = icount + 1
          partsInEachRank(i) = partsInEachRank(i) + 1
          reshuffle(partsInEachRank(i),i) = icount
          if(icount == nparts) exit
        endif
      enddo
      if(icount == nparts) exit
    enddo

    costperrankmax = 0.0d0
    costperrankmin = 1.0d+10

    do i=1,getNRanks()
      costperrank = 0.0d0
      do j=1,maxnparts
        !if(reshuffle(j,i)>0)write(*,*)i,j,reshuffle(j,i),gpat%sgraph(reshuffle(j,i))%lsize
        costperrank = costperrank + real((gpat%sgraph(reshuffle(j,i))%lsize)**3)
      enddo
      costperrankmax = max(costperrank,costperrankmax)
      costperrankmin = min(costperrank,costperrankmin)
    enddo
    !call gpmdcov_msII("gpmdcov_reshuffle","Measure of workload assymetry per rank &
    !       &"//to_string( (costperrankmax - costperrankmin)/costperrankmin ),lt%verbose,myRank)

  end subroutine gpmdcov_reshuffle

end module gpmdcov_reshuffle_mod
