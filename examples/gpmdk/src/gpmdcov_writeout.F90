module gpmdcov_writeout_mod

 use prg_extras_mod
 use prg_system_mod
 use prg_graph_mod
 use bml
 
 public :: gpmdcov_message, gpmdcov_msI, gpmdcov_msII
 public :: gpmdcov_msIII, gpmdcov_msMem, gpmdcov_msInt, gpmdcov_msVectInt, gpmdcov_msVectRel
 public :: gpmdcov_writepartitionout, gpmdcov_color_message, gpmdcov_error_message
 public :: gpmdcov_warning_message, gpmdcov_info_message, gpmdcov_status_message

contains

  !> Write a message with a predefined color
  subroutine gpmdcov_color_message(premsg,routine,message,color,verbose,verbTol,rank)
    implicit none
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: premsg
    character(*), intent(in) :: message
    character(*), intent(in) :: routine
    character(*), intent(in) :: color
    character(10) :: imsg, fmsg


    fmsg = '[0m'
    if(adjustl(trim(color)) == "red")then 
        imsg = '[31m'
    elseif(adjustl(trim(color)) == "green")then 
        imsg = '[32m'
    elseif(adjustl(trim(color)) == "yellow")then 
        imsg = '[33m'
    elseif(adjustl(trim(color)) == "blue")then 
        imsg = '[34m'
    else
        STOP 'Color not available'
    endif


    if (verbose >= verbTol) then
      if (present(rank)) then
        if(rank == 1)then
          write(*,*)""
          print*,achar(27),imsg,adjustl(trim(premsg))," at routine ",&
                  &adjustl(trim(routine)),"; ",adjustl(trim(message)),achar(27),fmsg
          write(*,*)""
        endif
      else
        write(*,*)""
        print*,achar(27),imsg,adjustl(trim(premsg))," at routine ",&
                &adjustl(trim(routine)),"; ",adjustl(trim(message)),achar(27),fmsg
        write(*,*)""
      endif
    endif

  end subroutine gpmdcov_color_message


  !> Write an error message
  !!
  subroutine gpmdcov_error_message(routine,message,verbose,verbTol,rank)
    implicit none
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_color_message(">>> ERROR!!!",routine,message,"red",verbose,verbTol,rank)    

  end subroutine gpmdcov_error_message 


  !> Write an warning message
  !!
  subroutine gpmdcov_warning_message(routine,message,verbose,verbTol,rank)
    implicit none
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_color_message("WARNING!!!",routine,message,"yellow",verbose,verbTol,rank)

  end subroutine gpmdcov_warning_message


  !> Write an information message
  !!
  subroutine gpmdcov_info_message(routine,message,verbose,verbTol,rank)
    implicit none
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_color_message(">>> INFO:",routine,message,"blue",verbose,verbTol,rank)

  end subroutine gpmdcov_info_message


  !> Write an status message
  !!
  subroutine gpmdcov_status_message(routine,message,verbose,verbTol,rank)
    implicit none
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_color_message(">>> STATUS:",routine,message,"green",verbose,verbTol,rank)

  end subroutine gpmdcov_status_message


  !> To write output file or perform some analysis
  !!
  subroutine gpmdcov_writepartitionout(sy,syprt,gpat,reshuffle,partsInEachRank,myRank)
  use prg_partition_mod
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    type(system_type)     ::  sy
    type(system_type), allocatable, intent(inout)    ::  syprt(:)
    type(graph_partitioning_t), intent(inout)        ::  gpat
    integer :: ipt,iptt,j
    integer, intent(in) :: myRank
    integer, allocatable, intent(in) :: partsInEachRank(:),reshuffle(:,:)
    character(25) :: filename
    character(20) :: auxchar

    if(allocated(sy%resindex))deallocate(sy%resindex)
    allocate(sy%resindex(sy%nats))
    sy%resindex=0

#ifdef DO_MPI
   !      !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
    do iptt=1,partsInEachRank(myRank)
       ipt=reshuffle(iptt,myRank)
#else
    do ipt=1,gpat%TotalParts
#endif

      write(filename,*)ipt
      auxchar = adjustl(trim(filename))
      filename = "part_"//auxchar
      if(.not.allocated(syprt(ipt)%net_charge)) then
        allocate(syprt(ipt)%net_charge(gpat%sgraph(ipt)%llsize))
        syprt(ipt)%net_charge = 0.0_dp
      endif

      call prg_write_system(syprt(ipt),filename,"pdb")

    enddo

    if(myRank == 1)then 
      do ipt = 1,gpat%TotalParts
        do j=1,gpat%sgraph(ipt)%llsize
          sy%resindex(gpat%sgraph(ipt)%core_halo_index(j)+1) = ipt
        enddo
      enddo
      if(.not.allocated(sy%net_charge))then 
        allocate(sy%net_charge(sy%nats))
        sy%net_charge = 0.0_dp
      endif
      if(myRank == 1) call prg_write_system(sy,"system_parts","pdb")
    endif

    deallocate(sy%resindex)

  end subroutine gpmdcov_writepartitionout

  subroutine gpmdcov_message(routine,message,verbose,verbTol,rank)
    implicit none 
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    if (verbose >= verbTol) then
      if (present(rank)) then
        if(rank == 1)then
          write(*,*)""
          write(*,*)"At routine ",routine,"; ",message
          write(*,*)""
        endif
      else
        write(*,*)""
        write(*,*)"At routine ",routine,"; ",message
        write(*,*)""
      endif
    endif

  end subroutine gpmdcov_message

  subroutine gpmdcov_msI(routine,message,verbose,rank)
    implicit none 
    integer, intent(in) :: verbose
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_message(routine,message,verbose,1,rank)

  end subroutine gpmdcov_msI

  subroutine gpmdcov_msII(routine,message,verbose,rank)
    implicit none 
    integer, intent(in) :: verbose
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_message(routine,message,verbose,2,rank)

  end subroutine gpmdcov_msII

  subroutine gpmdcov_msIII(routine,message,verbose,rank)
    implicit none 
    integer, intent(in) :: verbose
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(*), intent(in) :: routine

    call gpmdcov_message(routine,message,verbose,3,rank)

  end subroutine gpmdcov_msIII

  subroutine gpmdcov_msMem(routine,message,verbose,rank)
    implicit none 
    integer, intent(in) :: verbose
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    character(200) :: message1
    character(*), intent(in) :: routine

   if(verbose >= 4 .and. rank == 1)then
        message1 = "At "//trim(adjustl(routine))//" "//trim(adjustl(message))
        call prg_get_mem(routine,message1) 
   endif

  end subroutine gpmdcov_msMem


 subroutine gpmdcov_msRel(message,rel,verbose,verbTol,rank)
    implicit none 
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    real(8), intent(in) :: rel

    if (verbose >= verbTol) then
      if (present(rank)) then
        if(rank == 1)then
         write(*,*)""
         write(*,*)message," ",rel
         write(*,*)""
        endif
      else
       write(*,*)""
       write(*,*)message," ",rel
       write(*,*)""
      endif
    endif

  end subroutine gpmdcov_msRel

 subroutine gpmdcov_msInt(message,inte,verbose,verbTol,rank)
    implicit none 
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    character(*), intent(in) :: message
    integer, intent(in) :: inte

    if (verbose >= verbTol) then
      if (present(rank)) then
        if(rank == 1)then
         write(*,*)""
         write(*,*)message," ",inte
         write(*,*)""
        endif
      else
       write(*,*)""
       write(*,*)message," ",inte
       write(*,*)""
      endif
    endif

  end subroutine gpmdcov_msInt

 subroutine gpmdcov_msVectRel(message,rel,verbose,verbTol,rank)
    implicit none 
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    integer :: i
    character(*), intent(in) :: message
    real(8), allocatable, intent(in) :: rel(:)

    if (verbose >= verbTol) then
      if (present(rank)) then
        if(rank == 1)then
         write(*,*)""
         write(*,*)message
         do i = 1,size(rel,dim=1)
           write(*,*)rel(i)
         enddo
         write(*,*)""
        endif
      else
         write(*,*)message
         do i = 1,size(rel,dim=1)
           write(*,*)rel(i)
         enddo
         write(*,*)""
      endif
    endif

  end subroutine gpmdcov_msVectRel

 subroutine gpmdcov_msVectInt(message,inte,verbose,verbTol,rank)
    implicit none 
    integer, intent(in) :: verbose, verbTol
    integer, optional, intent(in) :: rank
    integer :: i
    character(*), intent(in) :: message
    integer, allocatable, intent(in) :: inte(:)

    if (verbose >= verbTol) then
      if (present(rank)) then
        if(rank == 1)then
         write(*,*)""
         write(*,*)message
         do i = 1,size(inte,dim=1)
           write(*,*)inte(i)
         enddo
         write(*,*)""
        endif
      else
         write(*,*)message
         do i = 1,size(inte,dim=1)
           write(*,*)inte(i)
         enddo
         write(*,*)""
      endif
    endif

  end subroutine gpmdcov_msVectInt

end module gpmdcov_writeout_mod
