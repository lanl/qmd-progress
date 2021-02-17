
  !> Initialize the program variables and parse input files.
  !!
  subroutine gpmdcov_Init()
    use gpmdcov_vars

    !> Start progress
    call prg_progress_init()

    !> Get MPI rank
    myRank = getMyRank() + 1

    if (printRank()  ==  1) then
      write(*,*)"" ; write(*,*) "GPMD started ..."; write(*,*)""
    endif

    !> Get the input file from argumets.
    call getarg(1, inputfile)
    if (printRank()  ==  1) then
      write(*,*)""; write(*,*)"Reading ",inputfile,"..."; write(*,*)""
    endif

    !> Parsing input file. This file contains all the variables needed to
    !  run the scf including the sp2 (solver) variables. lt is "latte_type" structure
    !  containing all the variables.
    !  file://~/progress/build/doc/html/structlatteparser__latte__mod_1_1latte__type.html
    call parse_latte(lt,inputfile)

    !> Parsing SP2 input paramenters. This will read the variables in the input file.
    !  sp2 is the "sp2data_type".
    call prg_parse_sp2(sp2,inputfile)

    !> Parsing GSP2 input paramenters. This will read the variables in the input file.
    !  gsp2 is the "gsp2data_type".
    call prg_parse_gsp2(gsp2,inputfile)

    !> Parsing Extended Lagrangian input paramenters. This will read the variables in the input file.
    !  xl is the "xlbo_type".
    call prg_parse_xlbo(xl,inputfile)

    !> Parsing Z sparse propagation.
    !  note: no need to pass a structure.
    ! call genZSPsolver%parse(zsp,inputfile)
    call prg_parse_zsp(zsp,inputfile)

    !> Parsing system coordinates. This reads the coords.pdb file to get the position of every
    !  atom in the system. sy is the "system_type" structure containing all the variables.
    !  file://~/progress/build/doc/html/classsystem__latte__mod.html
    call prg_parse_system(sy,lt%coordsfile)

    if(gsp2%mdim > 0)then
      mdim = gsp2%mdim
    else
      mdim = sy%nats
    endif

    !> Center sytem inside the box and fold it by the lattice_vectors.
    call prg_translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin)

    if (myRank  ==  1) then
      if(lt%verbose >= 2) call prg_write_system(sy,adjustl(trim(lt%jobname))//"_centered","pdb")
    endif

    !> Get the Coulombic cut off.
    call get_coulcut(lt%coul_acc,lt%timeratio,sy%nats,sy%lattice_vector,coulcut)

    !  if(lt%nlisteach > 1 .and. &
    !    min(sy%lattice_vector(1,1),sy%lattice_vector(2,2),sy%lattice_vector(3,3))/2.0_dp < coulcut)then
    !    write(*,*)"STOP: Make NlisEach= 1 under LATTE{} in order to continue ..."
    !    stop
    !  endif

    if(lt%restart) call gpmdcov_restart()

    !> Building the neighbor list.
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","Before build_nlist")
    call build_nlist_int(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","After build_nlist")

    !> Get Huckel hamiltonian. Computes the Extended Huckel Hamiltonian from the
    !  atom coordinates. The main inputs are the huckelTBparams and the system coordinate (sy%coordinate)
    !  The main outputs are Hamiltonian (ham_bml) and Overlap (over_bml) matrices.
    !  call get_hshuckel(ham_bml,over_bml,sy%coordinate,sy%spindex,sy%spatnum,&
    !     "../../huckelTBparams",lt%bml_type,lt%mdim,lt%threshold&
    !     ,tb%nsp,tb%splist,tb%basis,tb%numel,tb%onsite_energ,&
    !     tb%norbi,tb%hubbardu)

    !> LATTE Hamiltonian parameter
    call load_latteTBparams(tb,sy%splist,lt%parampath)

    !> Get the reciprocal vectors
    call prg_get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)

    !> Bond integrals parameters for LATTE Hamiltonian.
    call load_bintTBparamsH(sy%splist,tb%onsite_energ,&
         typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,lt%parampath)
    call write_bintTBparamsH(typeA,typeB,&
         intKind,intPairsH,intPairsS,adjustl(trim(lt%jobname))//"_mybondints.nonorth")

    !> Load Pair potentials for LATTE TB.
    call load_PairPotTBparams(lt%parampath,sy%splist,ppot)

    !> Allocate bounds vector.
    allocate(gbnd(2))

    !> mdstep needs to be prg_initialized.
    mdstep = 0

    !> This is just to get the number of total orbitals
    call get_hindex(sy%spindex,tb%norbi,sy%estr%hindex,norb)
    sy%estr%norbs = norb

    !> Get occupation based on last shell population.
    !  WARNING: This could change depending on the TB method being used.
    sy%estr%nel = sum(element_numel(sy%atomic_number(:)),&
         & size(sy%atomic_number,dim=1))
    bndfilTotal = sy%estr%nel/(2.0_dp*norb)
    write(*,*)"Total Number of Electrons:",sy%estr%nel


    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","After gpmd_Init")

  end subroutine gpmdcov_Init

