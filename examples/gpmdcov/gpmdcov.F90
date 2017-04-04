!> High-level program to perform GRAPH-BASED QMD.
!!
!! \ingroup PROGRAMS
!! \brief  High-level program to perform GRAPH-BASED QMD with a DFTB Hamiltonian 
!!  using a full parallel Graph-based approach.
!! \note To test this program with the \verbatim run_test \endverbatim script.
!!
!! \author C. F. A. Negre
!! (cnegre@lanl.gov)
!!
!! \author S. Mniszewski
!! (smn@lanl.gov)
!! 
!! \author A. M. N. Niklasson
!! (amn@lanl.gov)
!!
!! \author M. E. Wall
!! (mewall@lanl.gov)
!!
!! Verbose levels: 
!!   >= 0- Only print tags.
!!   >= 1- Print useful scalar physical data (e.g., Total Energy)
!!   >= 2- Print single file output data (e.g., 0-SCF Charges) 
!!   >= 3- Write out trajectory data.
!!   >= 4- Write out physically meaningful 1D arrays (e.g., Charges for the species)
!!   >= 5- Write out physically meaningful 2D arrays (e.g., H matrix)
!!
program gpmd

  !BML lib.
  use bml
  !Progress and LATTE lib modules.
  use prg_progress_mod
  use prg_system_mod
  use prg_ptable_mod
  use latteparser_latte_mod
  use huckel_latte_mod
  use tbparams_latte_mod 
  use ham_latte_mod
  use coulomb_latte_mod
  use prg_charges_mod
  use prg_initmatrices_mod
  use prg_genz_mod
  use prg_nonortho_mod
  use prg_pulaymixer_mod
  use prg_dos_mod
  use prg_densitymatrix_mod
  use neighborlist_latte_mod
  use ppot_latte_mod  
  use hsderivative_latte_mod
  use slaterkosterforce_latte_mod
  use prg_pulaycomponent_mod
  use nonorthocoulombforces_latte_mod
  ! Graph partitioning modules
  use prg_parallel_mod
  use prg_timer_mod
  use prg_graphsp2parser_mod
  use prg_sp2parser_mod
  use prg_sp2_mod
  use prg_graph_mod
  use prg_subgraphLoop_mod
  use prg_homolumo_mod
  use prg_xlbo_mod
  use md_latte_mod
  use prg_partition_mod
  use prg_extras_mod 

  implicit none     

  integer, parameter                ::  dp = kind(1.0d0)
  integer                           ::  seed = 1
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable         ::  intKind(:)
  character(50)                     ::  inputfile, filename
  character(2)                      ::  auxchar
  integer                           ::  mdstep, Nr_SCF_It, i, icount, ierr
  integer                           ::  j, nel, norb, pp(100), nnodes, iii
  integer                           ::  nparts, niter=500, npat, ipt
  integer                           ::  ii, jj, iscf, norb_core
  integer                           ::  mdim
  integer, allocatable              ::  hindex(:,:), hnode(:), vectorint(:)
  integer, allocatable              ::  xadj(:), adjncy(:), CH_count(:)  
  integer, allocatable              ::  part(:), core_count(:), Halo_count(:,:)  
  real(dp)                          ::  C0, C1, C2, C3
  real(dp)                          ::  C4, C5, ECoul, EKIN
  real(dp)                          ::  EPOT, ERep, Energy, Etot
  real(dp)                          ::  F2V, KE2T, MVV2KE, M_prg_init
  real(dp)                          ::  TRRHOH, Temp, Time, alpha
  real(dp)                          ::  bndfil, cc, coulcut, dt
  real(dp)                          ::  dx, egap, ehomo, elumo
  real(dp)                          ::  kappa, scferror, traceMult, vv(100)
  real(dp)                          ::  sumCubes, maxCH, Ef, smooth_maxCH, pnorm=6
  real(dp)                          ::  dvdw, d, mls_i, Efstep
  real(dp), allocatable             ::  FPUL(:,:), FSCOUL(:,:), FTOT(:,:), PairForces(:,:)
  real(dp), allocatable             ::  SKForce(:,:), VX(:), VY(:), VZ(:), collectedforce(:,:)
  real(dp), allocatable             ::  charges_old(:), coul_forces(:,:), coul_forces_k(:,:), coul_forces_r(:,:)
  real(dp), allocatable             ::  coul_pot_k(:), coul_pot_r(:), dqin(:,:), dqout(:,:)
  real(dp), allocatable             ::  eigenvals(:), gbnd(:), n(:), n_0(:)
  real(dp), allocatable             ::  n_1(:), n_2(:), n_3(:), n_4(:), acceprat(:)
  real(dp), allocatable             ::  n_5(:), onsitesH(:,:), onsitesS(:,:), rhoat(:)
  real(dp), allocatable             ::  origin(:), row(:), row1(:), auxcharge(:), auxcharge1(:)
  real(dp), allocatable             ::  g_dense(:,:),tch
  type(bml_matrix_t)                ::  aux_bml, dH0x_bml, dH0y_bml, dH0z_bml
  type(bml_matrix_t)                ::  dSx_bml, dSy_bml, dSz_bml, eigenvects
  type(bml_matrix_t)                ::  g_bml, ham0_bml, ham_bml
  type(bml_matrix_t)                ::  over_bml, rho_bml, rhoat_bml
  type(bml_matrix_t)                ::  rhoh_bml, zmat_bml, gch_bml
  type(bml_matrix_t)                ::  copy_g_bml, gcov_bml, aux1_bml
  type(graph_partitioning_t)        ::  gp
  type(graph_partitioning_t)        ::  gpat
  type(gsp2data_type)               ::  gsp2
  type(intpairs_type), allocatable  ::  intPairsH(:,:), intPairsS(:,:)
  type(latte_type)                  ::  lt
  type(neighlist_type)              ::  nl
  type(ppot_type), allocatable      ::  ppot(:,:)
  type(sp2data_type)                ::  sp2
  type(system_type)                 ::  sy
  type(system_type), allocatable    ::  syprt(:)  
  type(tbparams_type)               ::  tb
  type(xlbo_type)                   ::  xl
  logical                           ::  first_part = .true.
  logical                           ::  converged = .false.

  type(bml_matrix_t)                :: ZK1_bml, ZK2_bml, ZK3_bml 
  type(bml_matrix_t)                :: ZK4_bml, ZK5_bml, ZK6_bml  
  integer                           :: igenz !Couter to keep track of the times zmat is computed.
  logical                           :: tZSP, restart
  type(genZSPinp)                   :: zsp

  integer, allocatable              :: xadj_cov(:), adjncy_cov(:), CH_count_cov(:)  
  integer, allocatable              :: part_cov(:), core_count_cov(:), Halo_count_cov(:,:)  
  integer                           :: vsize(2)
  integer                           :: nparts_cov, myRank


  !!!!!!!!!!!!!!!!!!!!!!!!
  !> Main program driver
  !!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the program variables and parse input files. 
  call gpmd_Init()

  !We give a first guess of the Fermi level.
  !For water system, this is -2.5 eV.  
  Ef =  lt%efermi
   
  !> Initial partition of the system based on the covalency graph. 
  !! This will need to be replace by a first SP2 algorithm to compute a 
  !! first density matrix.
  call gpmd_Part()
    
  !> Initialize partitions.   
  call gpmd_InitParts() 
! stop

  !> Comput first charges. 
  
  call gpmd_FirstCharges()
    
  !> First SCF loop up to maxscf.
  call gpmd_DM_Min(lt%maxscf,sy%net_charge,.true.)

  !> First calculation of energies and forces.
  call gpmd_EnergAndForces(sy%net_charge)

  !> Setup the Molecular Dynamics (MD) calculation.
  call gpmd_PrepareMD()

  !> Perform the MD simulation.
  call gpmd_MDloop()

  !> Finalize the program.
  call gpmd_Finalize()

contains


  !> Initialize the program variables and parse input files. 
  !!
  subroutine gpmd_Init()
    implicit none

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

    if(lt%nlisteach > 1 .and. &
      min(sy%lattice_vector(1,1),sy%lattice_vector(2,2),sy%lattice_vector(3,3))/2.0_dp < coulcut)then 
      write(*,*)"STOP: Make NlisEach= 1 under LATTE{} in order to continue ..."
      stop
    endif
    
    if(lt%restart) call gpmd_restart()
        
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
    
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","After gpmd_Init")

  end subroutine gpmd_Init

  
  !> Partition by systems
  !!
  subroutine gpmd_Part
    implicit none
    integer, allocatable :: graph_h(:,:)
    integer, allocatable :: graph_p(:,:)
    real(dp)             :: mls_ii
            
    if(mdstep < 1)then 
      if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "Before prg_get_covgraph")
      write(*,*)"MPI rank",myRank, "in prg_get_covgraph .."
      call prg_get_covgraph(sy,nl%nnStructMindist,nl%nnStruct,nl%nrnnstruct&
        ,gsp2%bml_type,gsp2%covgfact,g_bml,gsp2%mdim,lt%verbose)
!       call bml_write_matrix(g_bml,"g_bml")
      write(*,*)"MPI rank",myRank, "done with prg_get_covgraph .."   
      if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "After prg_get_covgraph")  
    else 
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Ander's way of graph construction.   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"In prg_get_covgraph_h .."
      mls_ii = mls()
      call prg_get_covgraph_h(sy,nl%nnStructMindist,nl%nnStruct,nl%nrnnstruct,gsp2%nlgcut,graph_h,gsp2%mdim,lt%verbose)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_get_covgraph_h = ", mls()-mls_ii
    
#ifdef DO_MPI
      do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif
        call prg_collect_graph_p(syprt(ipt)%estr%orho,gpat%sgraph(ipt)%llsize,sy%nats,syprt(ipt)%estr%hindex,&
          gpat%sgraph(ipt)%core_halo_index,graph_p,gsp2%gthreshold,gsp2%mdim,lt%verbose)
          
        call bml_deallocate(syprt(ipt)%estr%orho)
      enddo             
        
      mls_i = mls()


#ifdef DO_MPI
      if (getNRanks() > 1) then    
        call prg_sumIntReduceN(graph_p, mdim*sy%nats)        
      endif
#endif
    
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_sumIntReduceN for graph", mls() - mls_i
                    
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"In prg_merge_graph .."    
      mls_ii = mls()
      call prg_merge_graph(graph_p,graph_h)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_merge_graph = ", mls()-mls_ii
    
      deallocate(graph_h)!
     
      !Transform graph into bml format.
      if(bml_get_N(gcov_bml).GT.0) call bml_deallocate(g_bml)
!       call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,mdim,g_bml,lt%bml_dmode)  
      call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,mdim,g_bml)        
              
      if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","Before prg_graph2bml")          
      write(*,*)"MPI rank",myRank, "in prg_graph2bml .."          
      mls_ii = mls()
      call prg_graph2bml(graph_p,gsp2%bml_type,g_bml)       
    
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_graph2bml = ", mls()-mls_ii
      if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","After prg_graph2bml")  
        
    endif
    
    if(allocated(syprt))then 
      do ipt=1,gpat%TotalParts   
        call prg_destroy_subsystems(syprt(ipt),lt%verbose)
      enddo
      deallocate(syprt)
    endif  

    if (myRank  ==  1 .and. lt%verbose >= 5) then
      call bml_print_matrix("gcov",g_bml,0,15,0,15) 
    endif
      
    if(mod(mdstep,gsp2%parteach)==0.or.mdstep == 0 .or.mdstep == 1)then 
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"In graph_part .."
      mls_ii = mls()
      call gpmd_graphpart()
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_graphpart = ", mls()-mls_ii
      write(*,*)"MPI rank",myRank, "done with graph_part .."    
    endif
    
    !To partition by molecule.
    ! write(*,*)"part by mol"
    ! call prg_molpartition(sy,nparts_cov,nl%nnStructMindist,nl%nnStruct,nl%nrnnstruct,"O ",gpat)

    mls_ii = mls()    
    do i=1,gpat%TotalParts
      call bml_matrix2submatrix_index(g_bml,&
        gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
        gpat%sgraph(i)%core_halo_index, &
        vsize,.true.)
      gpat%sgraph(i)%lsize = vsize(1)
      gpat%sgraph(i)%llsize = vsize(2)
    enddo
    if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for bml_matrix2submatrix_index = ", mls()-mls_ii
    
    if(allocated(syprt))deallocate(syprt)    
    allocate(syprt(gpat%TotalParts))    

    !> For every partition get the partial CH systems.
    if(myRank == 1)then
        write(*,*)""; write(*,*)"Getting CH subsystems ..."; write(*,*)""
    endif
        
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","Before prg_get_subsystem")  
    mls_ii = mls()
#ifdef DO_MPI
    do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
    do ipt = 1,gpat%TotalParts
#endif
      call prg_get_subsystem(sy,gpat%sgraph(ipt)%lsize,gpat%sgraph(ipt)%core_halo_index,syprt(ipt))            
    enddo
    if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_get_subsystem = ", mls()-mls_ii 
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","After prg_get_subsystem")  
 
 
    !To analyze partitions with VMD.
    if(myRank == 1)then 
      if(mod(mdstep,20) == 0.or.mdstep == 2)then
        call gpmd_writeout()
      endif
    endif    
    
  end subroutine gpmd_Part


  !>  Initialize the partition. 
  !!
  subroutine gpmd_InitParts

      if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov","Before gpmd_InitParts")  
#ifdef DO_MPI
      do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif

      if(myRank == 1) then 
          write(*,*)""
          write(*,*)"         #######################"
          write(*,*)"           Initializing partition =",ipt
          write(*,*)"         #######################"
          write(*,*)""
        if(lt%verbose >= 1)then           
          write(*,*)""
          write(*,*)"Number of atoms in the core =      ",gpat%sgraph(ipt)%llsize
          write(*,*)"Number of atoms in the cores+halo =",gpat%sgraph(ipt)%lsize
          write(*,*)""
        endif  
      endif  

      !> Get the mapping of the Hamiltonian index with the atom index 
      if(allocated(syprt(ipt)%estr%hindex))deallocate(syprt(ipt)%estr%hindex)    
      allocate(syprt(ipt)%estr%hindex(2,syprt(ipt)%nats))    

      call get_hindex(syprt(ipt)%spindex,tb%norbi,syprt(ipt)%estr%hindex,norb)                              
      syprt(ipt)%estr%norbs = norb

      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%ham0)    
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%over)    

      call get_hsmat(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%coordinate,&
        syprt(ipt)%lattice_vector,syprt(ipt)%spindex,&
        tb%norbi,syprt(ipt)%estr%hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)

      if (myRank  ==  1 .and. lt%verbose >= 5)then 
          write(*,*)"H0 and S for part:"
          call bml_print_matrix("H0",syprt(ipt)%estr%ham0,0,6,0,6)
          call bml_print_matrix("S",syprt(ipt)%estr%over,0,6,0,6)
      endif

      !> Get occupation based on last shell population. 
      !  WARNING: This could change depending on the TB method being used.
      nel = sum(element_numel(syprt(ipt)%atomic_number(:)),&
        size(syprt(ipt)%atomic_number,dim=1))
      bndfil = nel/(2.0_dp*norb)

      !> Initialize the density matrix (rho_bml) and inverse overlap factor (zmat_bml).
      if(bml_get_N(syprt(ipt)%estr%zmat).GT.0) call bml_deallocate(syprt(ipt)%estr%zmat)
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%zmat)      

      !> Get the Inverse square root overlap matrix.
      if(lt%verbose >= 1) call prg_timer_start(dyn_timer,"Build Z for part")
      call gpmd_buildz(syprt(ipt)%estr%over,syprt(ipt)%estr%zmat)
      if(lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      if(myRank == 1 .and. lt%verbose >= 5)then 
        write(*,*)"Z matrix for part:"
        call bml_print_matrix("Z",syprt(ipt)%estr%zmat,0,6,0,6)
      endif

    enddo
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "After gpmd_InitParts")  
    
  end subroutine gpmd_InitParts


  !> First Charge computation.
  !! \brief Here we compute the first "non-scf charges".
  !!  
  subroutine gpmd_FirstCharges()
    implicit none

    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "Before gpmd_FirstCharges")  
    if(.not.allocated(sy%net_charge))allocate(sy%net_charge(sy%nats))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over parts
      !!!!!!!!!!!!!!!!!!!!!!!!!!!    
      sy%net_charge = 0.0_dp

#ifdef DO_MPI
      do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif

      norb = syprt(ipt)%estr%norbs

      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)                 
      
      !> Initialize the orthogonal versions of ham and rho.      
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)                 
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)      

      !> Orthogonalize ham.
      if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(ortho_timer)
      call prg_orthogonalize(syprt(ipt)%estr%ham0,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
        lt%threshold,lt%bml_type,lt%verbose)
      if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(ortho_timer)

      call gpmd_rhosolver(syprt(ipt)%estr%oham,syprt(ipt)%estr%orho)

      !> Deprg_orthogonalize rho.       
      if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(deortho_timer)
      call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
        lt%threshold,lt%bml_type,lt%verbose)
      if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(deortho_timer)

      !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing 
      !! the charges.
      call prg_get_charges(syprt(ipt)%estr%rho, syprt(ipt)%estr%over, syprt(ipt)%estr%hindex, syprt(ipt)%net_charge, tb%numel,& 
        syprt(ipt)%spindex, lt%mdim, lt%threshold)

      if(lt%verbose.GE.4 .and. myRank  ==  1)then 
        write(*,*)""
        write(*,*)"Total charge of part",ipt,"=",sum(syprt(ipt)%net_charge(:),size(syprt(ipt)%net_charge,dim=1))
        write(*,*)""
        write(*,*)""; write(*,*)"Part charges:"
        do j=1,syprt(ipt)%nats
          write(*,*)syprt(ipt)%symbol(j),syprt(ipt)%net_charge(j)
        enddo
      endif

      do j=1,gpat%sgraph(ipt)%llsize          
        jj = gpat%sgraph(ipt)%core_halo_index(j)+1
        sy%net_charge(jj) = syprt(ipt)%net_charge(j)
      enddo
       
      call bml_deallocate(syprt(ipt)%estr%oham)
      call bml_deallocate(syprt(ipt)%estr%orho)
      call bml_deallocate(syprt(ipt)%estr%rho)        
                    
    enddo 
    ! End of loop over parts.
    
    mls_i = mls()
    
#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
      call prg_sumRealReduceN(sy%net_charge, sy%nats)
    endif
#endif

    if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"MPI rank finished prg_sumIntReduceN for qs", mls() - mls_i
    
    !> Gather charges from all the parts.
    if(.not.allocated(charges_old))allocate(charges_old(sy%nats))

    if (printRank()  ==  1) then
      if(lt%verbose >= 2)then 
        write(*,*)""
        write(*,*)"Total charge of the system=",sum(sy%net_charge(:),size(sy%net_charge,dim=1))
        write(*,*)""                              
        if(lt%verbose >= 5) call prg_write_system(sy,"charged_system","pdb")            
        write(*,*)""; write(*,*)"Full System charges:"
        do j=1,sy%nats
          write(*,*)j,sy%symbol(j),sy%net_charge(j)
        enddo
      endif
    endif  

    sy%net_charge(:) = sy%net_charge(:)- sum(sy%net_charge(:))/real(sy%nats)
    
    charges_old = sy%net_charge
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "After gpmd_FirstCharges")  
  end subroutine gpmd_FirstCharges


  !>  SCF loop
  !!
  subroutine gpmd_DM_Min(Nr_SCF,nguess,mix)
    implicit none
    integer :: Nr_SCF
    real(dp), allocatable :: nguess(:)
    logical :: mix
    real(dp) :: tch, tch1

    converged = .false.
    charges_old = nguess
    
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "Before gpmd_DM_Min")  
    
    ! Beginning of the SCF loop.
    if(.not.allocated(auxcharge))allocate(auxcharge(sy%nats))
    
    do iscf=1,Nr_SCF    

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This is done for the whole system 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (myRank == 1) write(*,*)"SCF iter", iscf

      !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
      if (myRank == 1) write(*,*)"In real Coul ..."

      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Real coul")
!       call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate&
!         ,nguess,tb%hubbardu,sy%lattice_vector,&
!         sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
!         nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
        ,nguess,tb%hubbardu,sy%lattice_vector,&
        sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
         nl%nnIz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);        
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      if (myRank  ==  1) write(*,*)"In recip Coul ..."    

      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Recip coul")
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
        ,nguess,tb%hubbardu,sy%lattice_vector,&
        sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);  
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      if(iscf == Nr_SCF) converged = .true. 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over parts
      !!!!!!!!!!!!!!!!!!!!!!!!!!!

      auxcharge = 0.0_dp
      mls_i = mls()
#ifdef DO_MPI
      do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif

        norb = syprt(ipt)%estr%norbs

        if(.not.allocated(syprt(ipt)%estr%coul_pot_k))then 
          allocate(syprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
          allocate(syprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))
        endif  

          syprt(ipt)%estr%coul_pot_k = 0.0_dp
          syprt(ipt)%estr%coul_pot_r = 0.0_dp
          syprt(ipt)%net_charge = 0.0_dp         
         
        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize          
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          syprt(ipt)%estr%coul_pot_k(j) = coul_pot_k(jj)
          syprt(ipt)%estr%coul_pot_r(j) = coul_pot_r(jj)              
          syprt(ipt)%net_charge(j) = nguess(jj)              
        enddo
        
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%ham)    
        call bml_copy_new(syprt(ipt)%estr%ham0,syprt(ipt)%estr%ham)

        !> Get the scf hamiltonian. The outputs is ham_bml.
        if (myRank == 1) write(*,*)"In prg_get_hscf ..."
        call prg_get_hscf(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%estr%ham,syprt(ipt)%spindex,&
          syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
          syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)
          
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)                         
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)      

        !> Orthogonalize ham.
        if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(ortho_timer)
        call prg_orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
          lt%threshold,lt%bml_type,lt%verbose)
        if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(ortho_timer)

        !> Now sove for the desity matrix.
        call gpmd_rhosolver(syprt(ipt)%estr%oham,syprt(ipt)%estr%orho)

!         call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml)                 
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)                 
        
        !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
        if(printRank() == 1 .and. lt%verbose >= 1) call prg_timer_start(deortho_timer)
        call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
          lt%threshold,lt%bml_type,lt%verbose)          
        if(printRank() == 1 .and. lt%verbose >= 1) call prg_timer_stop(deortho_timer)

        !> Get the system charges from rho
        call prg_get_charges(syprt(ipt)%estr%rho,syprt(ipt)%estr%over,syprt(ipt)%estr%hindex,syprt(ipt)%net_charge,tb%numel,&
          syprt(ipt)%spindex,lt%mdim,lt%threshold)

        if(lt%verbose >= 3) write(*,*)"Total charge of the part =", printRank(),&                
          sum(syprt(ipt)%net_charge(:)),size(syprt(ipt)%net_charge,dim=1)

        if(lt%verbose >= 3) write(*,*)"Total charge of the core part =", printRank() ,&                
          sum(syprt(ipt)%net_charge(1:gpat%sgraph(ipt)%llsize))
 
        do ii=1,gpat%sgraph(ipt)%llsize          
          j = gpat%sgraph(ipt)%core_halo_index(ii)+1
          auxcharge(j) = syprt(ipt)%net_charge(ii)
        enddo
        
      enddo  

      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for get qs of all parts", mls() - mls_i
      
      mls_i = mls()
      
#ifdef DO_MPI
      if (getNRanks() > 1) then
        call prg_sumRealReduceN(auxcharge(:), sy%nats)
      endif
#endif
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"MPI rank finished prg_sumRealReduceN for qs", mls() - mls_i
      
      nguess = auxcharge
      
!       tch = sum(nguess)            
!       nguess(:) = nguess(:) - tch/real(sy%nats)        
      
      if(mix)then 
        call prg_qmixer(nguess,charges_old,dqin,&
        dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,lt%verbose)
!         call prg_linearmixer(nguess,charges_old,scferror,lt%mixcoeff,lt%verbose)
      endif

      charges_old = nguess
          
      !Actualize the Fermi level.    
      tch = sum(nguess)       
      
      if(.not.allocated(acceprat))then 
        allocate(acceprat(2))
        acceprat = 0
        Efstep = 0.1_dp
      endif     

!       acceprat(5) = acceprat(4)
!       acceprat(4) = acceprat(3)      
!       acceprat(3) = acceprat(2)
      acceprat(2) = acceprat(1)
      acceprat(1) = sign(1.0_dp,tch)
      
!       if(mod(iscf,5)==0)then 
      if(acceprat(2)*acceprat(1) < 0)then 
        Efstep = Efstep*0.8_dp      
      else
        Efstep = Efstep*1.01_dp      
      endif        
      
      if(Nr_SCF.gt.10)then 
        if(iscf.gt.10)Ef = Ef -  sign(1.0_dp,tch)*min(tch**2,Efstep)
      else
        Ef = Ef -  sign(1.0_dp,tch)*min(tch**2,Efstep)
      endif  
!       endif
!       Ef = Ef -  sign(1.0_dp,tch1)*Efstep
      !Normalize charges to tch
      nguess(:) = nguess(:) - tch/real(sy%nats)
      
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Ef,q",Ef, tch, Efstep
            
      if(converged)then ! To do a last extra step.
        exit
      else  
        if(scferror < lt%scftol .and. iscf > 1) then 
          if (myRank  ==  1) then     
            write(*,*)""; write(*,*)"SCF converged within",iscf,"steps ..."
            write(*,*)"SCF error =",scferror
          endif
          converged = .true.
        endif 
      endif  
      
    enddo
    !> End of SCF loop.

    deallocate(auxcharge)

    if(myRank == 1 .and. lt%verbose >= 2)then
      write(*,*)""
      write(*,*)"Total charge of the system (After SCF) =",sum(nguess(:),size(nguess,dim=1))
      write(*,*)""                              
      call prg_write_system(sy,"charged_system","pdb")
      write(*,*)""; write(*,*)"System charges (After SCF):"
      do j=1,sy%nats
        write(*,*)sy%symbol(j),nguess(j)
      enddo
    endif  

    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "After gpmd_DM_Min")  
        
  end subroutine gpmd_DM_Min


  subroutine gpmd_EnergAndForces(charges)
    Implicit none
    real(dp), intent(in) :: charges(:)
    real(dp), allocatable :: ebandvector(:)
    
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "Before gpmd_EnergAndForces")
      
    if(.not.allocated(coul_forces)) allocate(coul_forces(3,sy%nats))    
    if(.not.allocated(FPUL))allocate(FPUL(3,sy%nats))    
    if(.not.allocated(FSCOUL))allocate(FSCOUL(3,sy%nats))    
    if(.not.allocated(SKForce))allocate(SKForce(3,sy%nats))
    if(.not.allocated(collectedforce))allocate(collectedforce(3,sy%nats))      
    if(.not.allocated(ebandvector))allocate(ebandvector(gpat%TotalParts))    
    
    FPUL = 0.0_dp
    FSCOUL = 0.0_dp
    SKForce = 0.0_dp
    collectedforce = 0.0_dp
    ebandvector = 0.0_dp
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Loop over all the parts  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DO_MPI
      do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif

      !Distribute the charges back to the parts.
      do j=1,gpat%sgraph(ipt)%lsize          
        jj = gpat%sgraph(ipt)%core_halo_index(j)+1
        syprt(ipt)%net_charge(j) = charges(jj)
      enddo
  
      norb = syprt(ipt)%estr%norbs

      norb_core = syprt(ipt)%estr%hindex(2,gpat%sgraph(ipt)%llsize) 

      if(bml_get_N(aux_bml).gt.0)then 
        call bml_deallocate(aux_bml)
        call bml_deallocate(aux1_bml)  
        deallocate(row)    
      endif

      allocate(row(norb))

      !> Get Electronic energy
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,nOrb,nOrb,aux1_bml)   
      call bml_copy_new(syprt(ipt)%estr%rho,aux_bml)


      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,nOrb,nOrb,rhoat_bml)
      call prg_build_atomic_density(rhoat_bml,tb%numel,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,norb,&
        lt%bml_type)

      
      call bml_add_deprecated(1.0_dp,aux_bml,-1.0_dp,rhoat_bml,lt%threshold)
      call bml_multiply(aux_bml,syprt(ipt)%estr%ham,aux1_bml,1.0d0, 0.0d0,lt%threshold)
      row=0.0_dp
      call bml_deallocate(rhoat_bml)
      call bml_get_diagonal(aux1_bml,row)

      TRRHOH = 0.0_dp
      do i=1,norb_core
        TRRHOH= TRRHOH+ row(i)
      enddo

      if(printRank() == 1 .and. lt%verbose >= 5) write(*,*)"Energy Band for part = ",ipt,"= ", TRRHOH   

      call bml_deallocate(aux_bml)
      call bml_deallocate(aux1_bml)
      call bml_deallocate(syprt(ipt)%estr%oham)      
      deallocate(row)    

      syprt(ipt)%estr%eband = TRRHOH
      ebandvector(ipt) = TRRHOH

      dx = 0.0001_dp;

      call get_dH(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,intPairsH,onsitesH,syprt(ipt)%symbol,&
        syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
        lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)
  
      call get_dS(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,intPairsS,onsitesS,syprt(ipt)%symbol,&
        syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
        lt%threshold, dSx_bml,dSy_bml,dSz_bml)

      if(printRank() == 1 .and. lt%verbose >= 10)then
        call bml_print_matrix("dH0x_bml",dH0x_bml,0,10,0,10)   
        call bml_print_matrix("dH0y_bml",dH0y_bml,0,10,0,10)   
        call bml_print_matrix("dH0z_bml",dH0z_bml,0,10,0,10)     

        call bml_print_matrix("dSx_bml",dSx_bml,0,10,0,10)   
        call bml_print_matrix("dSy_bml",dSy_bml,0,10,0,10)   
        call bml_print_matrix("dSz_bml",dSz_bml,0,10,0,10)     
      endif

      call get_skforce(syprt(ipt)%nats,syprt(ipt)%estr%rho,dH0x_bml,dH0y_bml,&
        dH0z_bml,syprt(ipt)%estr%hindex,syprt(ipt)%estr%SKForce,lt%threshold)
 
      call prg_get_pulayforce(syprt(ipt)%nats,syprt(ipt)%estr%zmat,syprt(ipt)%estr%ham,syprt(ipt)%estr%rho,&
        dSx_bml,dSy_bml,dSz_bml,syprt(ipt)%estr%hindex,syprt(ipt)%estr%FPUL,lt%threshold)

      call get_nonortho_coul_forces(syprt(ipt)%nats, norb, dSx_bml,dSy_bml,dSz_bml,&  
        syprt(ipt)%estr%hindex,syprt(ipt)%spindex,syprt(ipt)%estr%rho,syprt(ipt)%net_charge,syprt(ipt)%estr%coul_pot_r,&
        syprt(ipt)%estr%coul_pot_k,tb%hubbardu,syprt(ipt)%estr%FSCOUL,lt%threshold)

      do i=1,gpat%sgraph(ipt)%llsize
        FPUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FPUL(:,i)
        FSCOUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FSCOUL(:,i)
        SKForce(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%SKForce(:,i)
      enddo
                  
      call bml_deallocate(dSx_bml)
      call bml_deallocate(dSy_bml)
      call bml_deallocate(dSz_bml)

      call bml_deallocate(dH0x_bml)
      call bml_deallocate(dH0y_bml)
      call bml_deallocate(dH0z_bml)
      
      call bml_deallocate(syprt(ipt)%estr%rho)
      call bml_deallocate(syprt(ipt)%estr%ham)
      call bml_deallocate(syprt(ipt)%estr%ham0)      
      call bml_deallocate(syprt(ipt)%estr%over)
      call bml_deallocate(syprt(ipt)%estr%zmat)            
                        
    enddo
     
    collectedforce = FPUL + FSCOUL + SKForce   

    mls_i = mls()
          
#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
!       call prg_sumRealReduceN(collectedforce(1,:), sy%nats)
!       call prg_sumRealReduceN(collectedforce(2,:), sy%nats)
!       call prg_sumRealReduceN(collectedforce(3,:), sy%nats)      
      call prg_sumRealReduceN(collectedforce, sy%nats*3)
      call prg_sumRealReduceN(ebandvector, gpat%TotalParts)     
    endif
#endif  
      
    if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"MPI rank finished prg_sumRealReduceN for Forces", mls() - mls_i
    
    coul_forces =  coul_forces_r + coul_forces_k   

    !> Get Repulsive energy and forces
    if(lt%verbose >= 1) call prg_timer_start(dyn_timer,"Get pair pot")
!     call get_PairPot_contrib(sy%coordinate,sy%lattice_vector,sy%spindex,ppot,PairForces,ERep)
    call get_PairPot_contrib_int(sy%coordinate,sy%lattice_vector,nl%nnIx,nl%nnIy,&
         nl%nnIz,nl%nrnnlist,nl%nnType,sy%spindex,ppot,PairForces,ERep)    
    if(lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)
    write(*,*)"Energy Repulsive = ", ERep

    !> Get Coulombic energy
    ECoul = 0.0; 
    do i = 1,sy%nats
      ECoul = ECoul + charges(i)*(tb%hubbardu(sy%spindex(i))*charges(i) + coul_pot_r(i) + coul_pot_k(i) );
    enddo

    Etot = sum(ebandvector(:)) - 0.5_dp*ECoul  + ERep    

    if(myRank == 1 .and. lt%verbose >= 1)then         
      write(*,*)"Energy Coulomb = ", ECoul
      write(*,*)"Energy Band =", sum(ebandvector(:))
      write(*,*)"Energy Electronic =", Etot
   endif

    EPOT = Etot;

    if(.not.allocated(sy%force))allocate(sy%force(3,sy%nats))

    sy%force =  collectedforce + coul_forces + PairForces    

    if(myRank == 1 .and. lt%verbose >= 2)then 
      write(*,*)""; write(*,*)"FPUL + FSCOUL + SKForce"
      do i = 1,sy%nats
        write(*,*)"Coul Force",i,collectedforce(1,i),collectedforce(2,i),collectedforce(3,i)
      enddo

      write(*,*)""; write(*,*)"Coulomb forces:"
      do i = 1,sy%nats
        write(*,*)"Coul Force",i,coul_forces(1,i),coul_forces(2,i),coul_forces(3,i)
      enddo

      write(*,*)""; write(*,*)"Repulsive forces:"
      do i = 1,sy%nats
        write(*,*)i,PairForces(1,i),PairForces(2,i),PairForces(3,i)
      enddo

      write(*,*)""; write(*,*) "Total forces"
      do i = 1,sy%nats
        write(*,*)i,sy%force(1,i),sy%force(2,i),sy%force(3,i)
      enddo
    endif
    
    deallocate(ebandvector)
    
    if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "After gpmd_EnergAndForces")
    
  end subroutine gpmd_EnergAndForces

  !>  Preparing for MD
  !!
  subroutine gpmd_prepareMD()
    Implicit none

    !> Initialize velocities
    if(.not.allocated(sy%velocity))then 
      allocate(sy%velocity(3,sy%nats))
      sy%velocity(1,:) = 0.0_dp; sy%velocity(2,:) = 0.0_dp; sy%velocity(3,:) = 0.0_dp;    ! Initialize velocities
    endif
      
    !> Kinetic energy in eV (MVV2KE: unit conversion)
    MVV2KE = 166.0538782_dp/1.602176487_dp

    KE2T = 1.0_dp/0.000086173435_dp    

  end subroutine gpmd_prepareMD

  !>  Main MD loop
  !!  This routine performs the MD loops up to "ls%mdsteps"
  !! 
  subroutine gpmd_MDloop()
    implicit none
    real(dp) :: mls_ii

    do mdstep = 1,lt%mdsteps
    
    mls_ii = mls()
    
    if(myRank == 1)then 
      write(*,*)""
      write(*,*)"         #######################"
      write(*,*)"           MDStep =",mdstep
      write(*,*)"         #######################"
      write(*,*)""
    endif
    
      !> Get Kinetic energy
      EKIN = 0.0_dp
      do i=1,sy%nats
        EKIN = EKIN + sy%mass(i)*(sy%velocity(1,i)**2+sy%velocity(2,i)**2+sy%velocity(3,i)**2)   
      enddo    
      EKIN = 0.5_dp*MVV2KE*EKIN

      !! Statistical temperature in Kelvin
      Temp = (2.0_dp/3.0_dp)*KE2T*EKIN/real(sy%nats,dp);        
      !! Total Energy in eV
      Energy = EKIN + EPOT;        
      !! Time in fs
      Time = mdstep*lt%timestep;  

      if(myRank == 1)then 
        write(*,*)"Time [fs] = ",Time
        write(*,*)"Energy Kinetic [eV] = ",EKIN
        write(*,*)"Energy Potential [eV] = ",EPOT
        write(*,*)"Energy Total [eV] = ",Energy
        write(*,*)"Temperature [K] = ",Temp
      endif
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Preliminars", mls() - mls_ii
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul1", mls() - mls_ii
      
      !> First 1/2 of Leapfrog step
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Half Verlet")
      call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))        
      if(lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      if(myRank == 1 .and. lt%verbose.GE.5)then
        write(*,*)"Velocities"
        do i = 1,sy%nats
          write(*,*)i,sy%velocity(1,i),sy%velocity(2,i),sy%velocity(3,i)
        enddo
      endif
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul2", mls() - mls_ii
      
      !> Update positions
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Update positions")
      call updatecoords(origin,sy%lattice_vector,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:),sy%coordinate)
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
        call prg_sumRealReduceN(sy%coordinate(1,:), sy%nats)
        call prg_sumRealReduceN(sy%coordinate(2,:), sy%nats)
        call prg_sumRealReduceN(sy%coordinate(3,:), sy%nats)      
        sy%coordinate = sy%coordinate/real(getNRanks(),dp)            
    endif
#endif  
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul3", mls() - mls_ii
      !> Update neighbor list (Actialized every nlisteach times steps)
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Build Nlist")
      if(mod(mdstep,lt%nlisteach) == 0 .or. mdstep == 0 .or. mdstep == 1)then 
        call destroy_nlist(nl)
        call build_nlist_int(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)        
      endif
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      !> Repartition.
      ! This builds the new graph.
      mls_i = mls()
      call gpmd_Part()
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_Part", mls() - mls_i
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul4", mls() - mls_ii
      !> Reprg_initialize parts.
      mls_i = mls()
      call gpmd_InitParts() 
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_InitParts", mls() - mls_i
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul5", mls() - mls_ii
      
      mls_i = mls()
      call prg_xlbo_nint(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_xlbo_nint", mls() - mls_i
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul5", mls() - mls_ii
                        
      mls_i = mls()
      Nr_SCF_It = xl%maxscfiter;

      !> Use SCF the first M_prg_init MD steps
      if(mdstep < xl%minit)then 
        Nr_SCF_It = xl%maxscfInitIter
      else
        Nr_SCF_It = xl%maxscfiter
      endif  

      !> SCF loop 
      if(Nr_SCF_It.ne.0)call gpmd_DM_Min(Nr_SCF_It,n,.true.)
      
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul6", mls() - mls_ii
      
      sy%net_charge = n
      
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_DM_Min_1", mls() - mls_i      

      mls_i = mls()
      write(*,*)"Aditional DM construction ..."
      call gpmd_DM_Min(1,sy%net_charge,.false.)      
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_DM_Min_2", mls() - mls_i      
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul7", mls() - mls_ii
      
      mls_i = mls()
      call gpmd_EnergAndForces(n)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_EnergAndForces", mls() - mls_i
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul8", mls() - mls_ii
      
      mls_i = mls()
      !> Adjust forces for the linearized XLBOMD functional
      call prg_xlbo_fcoulupdate(Coul_Forces,sy%net_charge,n)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul9", mls() - mls_ii
      
      !> Total XLBOMD force
!       sy%force = SKForce + PairForces + FPUL + Coul_Forces + FSCOUL;       
      sy%force = collectedforce + PairForces + Coul_Forces 
      
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul10", mls() - mls_ii
      !> Integrate second 1/2 of leapfrog step
      call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))        

      if(lt%verbose >= 3 .and. myRank == 1)then 
          call prg_write_trajectory(sy,mdstep,5,lt%timestep,"trajectory","pdb")
      endif
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Rest", mls() - mls_i
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul11", mls() - mls_ii
      
    if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for MD iter", mls() - mls_ii
      
    ! Save MD state each 120 steps  
    if(mod(mdstep,150) == 0)call gpmd_dump()

    enddo
    ! End of MD loop.
    

  end subroutine gpmd_MDloop


  !> Solver for computing the density matrix 
  subroutine gpmd_RhoSolver(orthoh_bml,orthop_bml)
    implicit none
    type(bml_matrix_t), intent(in) :: orthoh_bml
    type(bml_matrix_t), intent(inout) :: orthop_bml
    
    if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"In solver ..."
    if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(dyn_timer,"Solver")

    if(lt%method.EQ."GSP2")then 
      call prg_timer_start(graphsp2_timer)
      call prg_subgraphSP2Loop(orthoh_bml, g_bml, orthop_bml, gp, lt%threshold)
      call prg_timer_stop(graphsp2_timer)
      ! call prg_sp2_alg1_seq(orthoh_bml,orthop_bml,lt%threshold, gp%pp, gp%maxIter, gp%vv)
    elseif(lt%method.EQ."SP2")then 
      call prg_sp2_alg2(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
        ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
    elseif(lt%method.EQ."Diag")then  
      ! call build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
      ! call prg_build_density_T(orthoh_bml,orthop_bml,lt%threshold,bndfil, 0.1_dp, Ef)
      call prg_build_density_T_Fermi(orthoh_bml,orthop_bml,lt%threshold, 0.1_dp, Ef)      
      if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"ipt =",ipt,"Ef =",Ef
    else
      stop"No valid Method in LATTE parameters"
    endif

    if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(dyn_timer,1)

! #ifdef DO_MPI_BLOCK
!     call prg_allGatherParallel(orthop_bml)
! #endif

    if(lt%verbose >= 2 .and. myRank == 1)then
      call bml_print_matrix("orthop_bml",orthop_bml,0,6,0,6)
    endif  

  end subroutine gpmd_RhoSolver


  !> Graph partitioning subroutine
  subroutine gpmd_graphpart
    implicit none
    integer :: tnnz
  
!     if(.not.allocated(xadj))then     
!     !> Symmetrize and Threshold the Matrix
!     if(gsp2%mdim > 0)then 
!       mdim = gsp2%mdim
!     else  
!       mdim = sy%nats
!     endif    
    
!     call bml_write_matrix(g_bml,"g_bml_bef")
!     call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,mdim,copy_g_bml)
!     call bml_threshold(g_bml, gsp2%gthreshold)
!     call bml_transpose(g_bml, copy_g_bml)
!     call bml_add_deprecated(0.5_dp,g_bml,0.5_dp,copy_g_bml,0.0_dp)
!     call bml_threshold(g_bml, gsp2%gthreshold)
!     call bml_deallocate(copy_g_bml)
!     call bml_write_matrix(g_bml,"g_bml_aft")
! #ifdef DO_MPI_BLOCK
!     call prg_allGatherParallel(g_bml)
! #endif

!     endif
    
    !> Create graph partitioning - Use Block or METIS or METIS+SA or METIS+KL
    if(lt%verbose >= 1) call prg_timer_start(part_timer)

    !> Block partitioning
    if (gsp2%partition_type == "Block") then
      !> Partition by orbital or atom
      if (gsp2%graph_element == "Orbital") then
        call prg_equalPartition(gpat, gsp2%nodesPerPart, sy%nats)
      else
        call prg_equalGroupPartition(gpat, hindex, nnodes, gsp2%nodesPerPart, sy%nats)
      endif

      !> METIS, METIS+SA, or METIS+KL partitioning
    else
      if(.not.allocated(xadj))then       
      
        allocate(xadj(sy%nats+1))
        
        tnnz = 0  !This must be done at the bml level
        do i=1,sy%nats
          tnnz = tnnz + bml_get_row_bandwidth(g_bml,i)
        enddo
        
!       allocate(adjncy(sy%nats*sy%nats)) !Old way 
        allocate(adjncy(tnnz+1))        
        
        call bml_adjacency(g_bml, xadj, adjncy, 1)        
        
!         call bml_deallocate(g_bml)
        
!          call prg_sortadj(xadj, adjncy)  !Not needed anymore since the g_bml is sorted.
        
        nnodes = sy%nats
      else  
        nnodes = sy%nats
      endif  

#ifdef DO_GRAPHLIB
      nparts = gsp2%partition_count

      if (first_part) then
        allocate(part(nnodes))
        allocate(core_count(nparts))
      endif

      ! if(allocated(CH_count)) deallocate(CH_count)
      ! if(allocated(Halo_count)) deallocate(Halo_count) 
      allocate(CH_count(nparts))
      allocate(Halo_count(nparts, nnodes))

      !> For METIS, if no refinement of first time, do full partitioning
      if (gsp2%partition_refinement == 'None' .or. first_part) then

        if(allocated(gpat%nnodesInPart))then 
          call prg_destroyGraphPartitioning(gpat)
        endif  
        !> Which METIS partitioning
        select case(gsp2%partition_type)
        case("METIS")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
            CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)                
        case("METIS+SA")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
            CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm) 
          call prg_simAnnealing(gpat, xadj, adjncy, part, core_count, CH_count, &
            Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
        case("METIS+KL")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
            CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm) 
          call prg_KernLin2(gpat, xadj, adjncy, part, core_count, CH_count, &
            Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
        case default
          write(*,*)"No METIS partitioning specified"
          stop ""
        end select

        first_part = .false.

        !> Not first time, do refinement
      else
        if (.not. first_part) then

          !> Which refinement
          select case(gsp2%partition_refinement)
          case("SA")
            call prg_simAnnealing(gpat, xadj, adjncy, part, core_count, CH_count, &
              Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
          case("KL")
            call prg_KernLin2(gpat, xadj, adjncy, part, core_count, CH_count, &
              Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
          case default
            write(*,*)"No refinement specified"
            stop ""
          end select

        endif
      endif

      deallocate(xadj)
      deallocate(adjncy)
      deallocate(CH_count)
      deallocate(Halo_count)

#endif
    endif

    if(lt%verbose >= 1) call prg_timer_stop(part_timer)

  end subroutine gpmd_graphpart

  subroutine gpmd_buildz(overin_bml,zmatout_bml)
    implicit none
    type(bml_matrix_t), intent(inout) :: zmatout_bml
    type(bml_matrix_t), intent(inout) :: overin_bml

    igenz = igenz + 1

    if(lt%zmat == "ZSP")then !Congruence transformation.

      call prg_buildzsparse(overin_bml,zmatout_bml,igenz,lt%mdim,&
        lt%bml_type, zk1_bml,zk2_bml,zk3_bml&
        ,zk4_bml,zk5_bml,zk6_bml,zsp%nfirst,zsp%nrefi,zsp%nreff,&
        zsp%numthresi,zsp%numthresf,zsp%integration,zsp%verbose)

    else               

      !Build Z matrix using diagonalization (usual method).
      call prg_buildzdiag(overin_bml,zmatout_bml,lt%threshold,lt%mdim,lt%bml_type)

    endif

  end subroutine gpmd_buildz

  !> To write output file or perform some analysis
  !!
  subroutine gpmd_writeout()
    implicit none
    
    if(allocated(sy%resindex))deallocate(sy%resindex) 
    allocate(sy%resindex(sy%nats))     
    sy%resindex(sy%nats)=-100
    

#ifdef DO_MPI
    do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
#else
    do ipt = 1,gpat%TotalParts
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

  end subroutine gpmd_writeout
  
  !> Finalize gpmd
  subroutine gpmd_Finalize()
    implicit none
    
    deallocate(gbnd)

    call bml_deallocate(ham_bml)
    call bml_deallocate(ham0_bml)
!     call bml_deallocate(orthoh_bml)
    call bml_deallocate(g_bml)
    call bml_deallocate(rho_bml)
    call bml_deallocate(over_bml)
    call bml_deallocate(zmat_bml)

    ! Progress is done
    call prg_progress_shutdown()

  end subroutine gpmd_Finalize

  !> Dump GPMD
  subroutine gpmd_dump()
    implicit none    
    
    if(myRank == 1)then     
    write(*,*)"Restarting ..."

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
    
  end subroutine gpmd_dump

  !> Dump GPMD
  subroutine gpmd_restart()
    implicit none    
        
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

  end subroutine gpmd_restart
      
end program gpmd
