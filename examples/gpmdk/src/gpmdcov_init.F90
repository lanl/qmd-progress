module gpmdcov_Init_mod

contains

  !> Initialize the program variables and parse input files.
  !!
  subroutine gpmdcov_Init(lib_on)

    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    use gpmdcov_neighbor_mod
    use gpmdcov_parser_mod
    use gpmdcov_allocation_mod
    use gpmdcov_prepareMD_mod, only : gpmdcov_addVelocity

    integer :: ix, iy, iz, numreps
    type(system_type) :: syaux
    real(dp) :: minEdge

    logical, intent(in) :: lib_on

    !> Start progress
    !  This will initialize the MPI ranks
    if(.not. lib_init) call prg_progress_init()
    
    !> Get MPI rank
#ifdef DO_MPI
    myRank = getMyRank() + 1
    numRanks = getNRanks()
#else
    myRank = 1
    numRanks = 1
#endif

    call gpmdcov_msI("gpmdcov_Init","GPMD started ...",lt%verbose,myRank)

    !> Get the input file from argumets.
    if(lib_on)then
      inputfile = "gpmd_input.in"
    else
      call gpmdcov_msI("gpmdcov_Init","Reading inputfile ...",lt%verbose,myRank)
      call getarg(1,inputfile)
    endif

    if(.not. lib_init)then

    !> Parsing input file. This file contains all the variables needed to
    !  run the scf including the sp2 (solver) variables. lt is "latte_type" structure
    !  containing all the variables.
    call parse_latte(lt,inputfile)

    !> Parsing SP2 input paramenters. This will read the variables in the input file.
    !  sp2 is the "sp2data_type".
    call prg_parse_sp2(sp2,inputfile)

    !> Parsing GSP2 input paramenters. This will read the variables in the input file.
    !  gsp2 is the "gsp2data_type". These parameters control all the graph
    !  related variables incluind partition methos, number of parts, etc.
    call prg_parse_gsp2(gsp2,inputfile)

    if((gsp2%nx .ne. 0) .and. (gsp2%ny .ne. 0) .and. (gsp2%nz .ne. 0))then
       nparts = gsp2%nx*gsp2%ny*gsp2%nz
       if(nparts .ne. gsp2%partition_count)then 
               call gpmdcov_msI("gpmdcov_Init","!!!ERROR: If PartitionCountX/Y/Z are set, PartitionCount &
               &should be set to PartitionCountx * PartitionCountY * PartitionCountZ",lt%verbose,Myrank)
               stop
       endif
    endif


    !> Parsing Extended Lagrangian input paramenters. This will read all the
    ! variables related to the Extended lagrangian method.
    !  xl is the "xlbo_type".
    call prg_parse_xlbo(xl,inputfile)

    !> Parsing Z sparse propagation. This will read all the variables related to
    ! the propagation and construction of the inverse overlap.
    call prg_parse_zsp(zsp,inputfile)

    !> Parsing specific variales for the gpmd code
    call gpmdcov_parse(trim(adjustl(inputfile)),gpmdt)

    !> Parse variables for the kernel method
    call gpmdcov_parseKernel(kernel,inputfile)
    endif 


    !> Parsing system coordinates. This reads the coords.pdb file to get the position of every
    !  atom in the system. sy is the "system_type" structure containing all the variables.

    ! Reading the system
    if(gpmdt%replicatex + gpmdt%replicatey + gpmdt%replicatez > 0)then
      if(.not. lib_on) call prg_parse_system(syaux,lt%coordsfile)
      call prg_replicate_system(syaux,sy,gpmdt%replicatex,gpmdt%replicatey,gpmdt%replicatez)
      call gpmdcov_reallocate_intVect(sy%atomic_number,sy%nats)
      call gpmdcov_reallocate_realVect(sy%mass,sy%nats)
      call gpmdcov_reallocate_charVect(sy%splist,syaux%nsp)
      call gpmdcov_reallocate_realVect(sy%spmass,syaux%nsp)
      call gpmdcov_reallocate_intVect(sy%spatnum,syaux%nsp)
      sy%nsp = syaux%nsp
      sy%splist = syaux%splist
      sy%spmass = syaux%spmass
      sy%spatnum = syaux%spatnum

      if(gpmdt%htod)then
         numh = 0
         do i = 1,sy%nats
            if(syaux%atomic_number(i) .eq. 1)then
               syaux%mass(i) = 2.013553212745 ! Change H mass to D mass
               numh = numh + 1
            endif
         enddo
         call gpmdcov_msI("gpmdcov_init","Changing H mass to D mass for "&
                              &//to_string(numh)//" H atoms",lt%verbose,myRank)
      endif
      
      numreps = 0
      do ix = 1,gpmdt%replicatex
        do iy = 1,gpmdt%replicatey
          do iz = 1,gpmdt%replicatez
            numreps = numreps + 1
            sy%atomic_number(1 + (numreps - 1)*syaux%nats:numreps*syaux%nats) = syaux%atomic_number(:)
            sy%mass(1 + (numreps - 1)*syaux%nats:numreps*syaux%nats) = syaux%mass(:)
          enddo
        enddo
      enddo

      call prg_destroy_system(syaux)
      if (myRank  ==  1) then
        if(lt%verbose >= 2) call prg_write_system(sy,adjustl(trim(lt%jobname))//"_postReplication","pdb")
      endif
      
    else
    
       if(.not. lib_on) call prg_parse_system(sy,lt%coordsfile)

       if(gpmdt%htod)then
          numh = 0
          do i = 1,sy%nats
             if(sy%atomic_number(i) .eq. 1)then
                sy%mass(i) = 2.013553212745 ! Change H mass to D mass
                numh = numh + 1
             endif
          enddo
          call gpmdcov_msI("gpmdcov_init","Changing H mass to D mass for "&
               &//to_string(numh)//" H atoms",lt%verbose,myRank)
       endif
    endif

    call gpmdcov_msI("gpmdcov_init","Number of atoms = "//to_string(sy%nats),lt%verbose,myRank)


    !> This variable sets the maximun number of non-zeros per row. I will be
    ! used when sparse matrices are used. Typically sparse matrices will only be
    ! used for storing the graph and neighbor list.
    if(gsp2%mdim > 0)then
      mdim = gsp2%mdim
    else
      mdim = sy%nats
    endif

    if(.not.allocated(origin)) allocate(origin(3))
    !> Center sytem inside the box and fold it by the lattice_vectors. This is
    ! done only for visualization purposes.
    !origin = 0.0_dp
    if(gpmdt%trfl)then
      call prg_translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin)
    endif
    !origin = 0.0_dp

    if(gpmdt%restartfromdump)then
       call gpmdcov_msI("gpmdcov_init","Restarting from dump file ",lt%verbose,myRank)
       call gpmdcov_restart()
    endif
    
    !> This will construct a pdb file with the system that was read by the code.
    ! It should be used to verify that the program indeed is reading the system
    ! of interest.
    if (myRank  ==  1) then
      if(lt%verbose >= 2) call prg_write_system(sy,adjustl(trim(lt%jobname))//"_centered","pdb")
    endif

    !> Get the Coulombic cut off.
    call get_coulcut(lt%coul_acc,lt%timeratio,sy%nats,sy%lattice_vector,coulcut)

    !> Building the neighbor list.
    call gpmdcov_msMem("gpmdcov_init","Before build_nlist",lt%verbose,myRank)
    mls_i = mls()

    minEdge = min(norm2(sy%lattice_vector(1,:)),norm2(sy%lattice_vector(2,:)),norm2(sy%lattice_vector(3,:)))
    if(coulcut < 2*minEdge) nlistSparse = .true.
    nlistSparse = .false.
    if(nlistSparse)then 
      call gpmdcov_msI("gpmdcov_init", "Doing Linear Scaling Neighbor list construction... ",lt%verbose,myRank)
      call gpmdcov_build_nlist_sparse_v2(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose,myRank,numRanks)
    else
      call gpmdcov_msI("gpmdcov_init", "Doing Full Neighbor list construction... ",lt%verbose,myRank)
      call gpmdcov_build_nlist_full(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose,myRank,numRanks)
    endif
    !call gpmdcov_build_nlist_sparse(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose,myRank,numRanks)
    !call build_nlist_int(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)
   ! LBox(1) = sy%lattice_vector(1,1)
   ! LBox(2) = sy%lattice_vector(2,2)
   ! LBox(3) = sy%lattice_vector(3,3)
    
   ! call gpmd_nearestneighborlist(nl%nrnnlist,nl%nndist,nl%nnRx,nl%nnRy,nl%nnRz,nl%nnType, &
   !               &sy%coordinate(1,:),sy%coordinate(2,:),sy%coordinate(3,:),LBox,coulcut,sy%nats,200)
    call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_nlist "&
                     &//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)
    call gpmdcov_msMem("gpmdcov_init","After build_nlist",lt%verbose,myRank)

    !> LATTE Hamiltonian parameter
    call load_latteTBparams(tb,sy%splist,lt%parampath)

    !> Get the reciprocal vectors
    call prg_get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)

    !> Bond integrals parameters for LATTE Hamiltonian.
    call gpmdcov_msMem("gpmdcov","Before load_bintTBparamsH",lt%verbose,myRank)
    call load_bintTBparamsH(sy%splist,tb%onsite_energ,&
         typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,lt%parampath)
    call write_bintTBparamsH(typeA,typeB,&
         intKind,intPairsH,intPairsS,adjustl(trim(lt%jobname))//"_mybondints.nonorth")

    !> Load Pair potentials for LATTE TB.
    call load_PairPotTBparams(lt%parampath,sy%splist,ppot,lt%verbose)
    !write(*,*)"ppot",ppot(1,1)%potparams(:)
    !stop

    !> Allocate bounds vector.
    allocate(gbnd(2))

    !> mdstep needs to be initialized.
    mdstep = 0

    !> This is just to get the number of total orbitals
    call get_hindex(sy%spindex,tb%norbi,sy%estr%hindex,norb,lt%verbose)
    sy%estr%norbs = norb
    call gpmdcov_msRel("Total Number of Orbitals:",real(sy%estr%norbs,dp),lt%verbose,myRank)

    !> Get total occupation
    !  WARNING: This could change depending on the TB method being used.
    
    sy%estr%nel = sum(element_numel(sy%atomic_number(:)))
    bndfilTotal = sy%estr%nel/(2.0_dp*norb)
    call gpmdcov_msRel("Total Number of Electrons:",real(sy%estr%nel,dp),lt%verbose,myRank)

    call gpmdcov_msMem("gpmdcov","After gpmd_Init",lt%verbose,myRank)

  end subroutine gpmdcov_Init

end module gpmdcov_Init_mod
