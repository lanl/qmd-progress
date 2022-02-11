module gpmdcov_Init_mod

contains

  !> Initialize the program variables and parse input files.
  !!
  subroutine gpmdcov_Init()

    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    integer :: thisJ

    !> Start progress
    !  This will initialize the MPI ranks
    call prg_progress_init()

    !> Get MPI rank
    myRank = getMyRank() + 1

    call gpmdcov_msI("gpmdcov_Init","GPMD started ...",lt%verbose,myRank)

    !> Get the input file from argumets.
    call gpmdcov_msI("gpmdcov_Init","Reading inputfile ...",lt%verbose,myRank)
    call getarg(1,inputfile)   

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

    !> Parsing Extended Lagrangian input paramenters. This will read all the
    ! variables related to the Extended lagrangian method.
    !  xl is the "xlbo_type".
    call prg_parse_xlbo(xl,inputfile)

    !> Parsing Z sparse propagation. This will read all the variables related to
    ! the propagation and construction of the inverse overlap. 
    call prg_parse_zsp(zsp,inputfile)

    !> Parsing system coordinates. This reads the coords.pdb file to get the position of every
    !  atom in the system. sy is the "system_type" structure containing all the variables.
    call prg_parse_system(sy,lt%coordsfile)

    !> This variable sets the maximun number of non-zeros per row. I will be
    ! used when sparse matrices are used. Typically sparse matrices will only be
    ! used for storing the graph and neighbor list. 
    if(gsp2%mdim > 0)then
      mdim = gsp2%mdim
    else
      mdim = sy%nats
    endif

    !> Center sytem inside the box and fold it by the lattice_vectors. This is
    ! done only for visualization purposes. 
    call prg_translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin)

    !> This will construct a pdb file with the system that was read by the code.
    ! It should be used to verify that the program indeed is reading the system
    ! of interest.
    if (myRank  ==  1) then
      if(lt%verbose >= 2) call prg_write_system(sy,adjustl(trim(lt%jobname))//"_centered","pdb")
    endif

    !> Get the Coulombic cut off.
    call get_coulcut(lt%coul_acc,lt%timeratio,sy%nats,sy%lattice_vector,coulcut)
   
    !> Building the neighbor list.
    call gpmdcov_msMem("gpmdcov","Before build_nlist",lt%verbose,myRank)
    call build_nlist_int(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)
    call gpmdcov_msMem("gpmdcov","After build_nlist",lt%verbose,myRank)

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
    call load_PairPotTBparams(lt%parampath,sy%splist,ppot)

    !> Allocate bounds vector.
    allocate(gbnd(2))

    !> mdstep needs to be initialized.
    mdstep = 0

    !> This is just to get the number of total orbitals
    call get_hindex(sy%spindex,tb%norbi,sy%estr%hindex,norb,lt%verbose)
    sy%estr%norbs = norb
    call gpmdcov_msRel("Total Number of Orbitals:",real(sy%estr%norbs,dp),lt%verbose,myRank)

    !> Get occupation based on last shell population.
    !  WARNING: This could change depending on the TB method being used.
    sy%estr%nel = sum(element_numel(sy%atomic_number(:)),&
         & size(sy%atomic_number,dim=1))
    bndfilTotal = sy%estr%nel/(2.0_dp*norb)
    call gpmdcov_msRel("Total Number of Electrons:",real(sy%estr%nel,dp),lt%verbose,myRank)
    
    !> Parse variables for the kernel method
    call gpmdcov_parseKernel(kernel,inputfile)

    call gpmdcov_msMem("gpmdcov","After gpmd_Init",lt%verbose,myRank)

  end subroutine gpmdcov_Init

end module gpmdcov_Init_mod
