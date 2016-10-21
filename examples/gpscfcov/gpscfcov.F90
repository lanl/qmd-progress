!> High-level program to perform QMD with a DFTB Hamiltonian.
!!
!! \ingroup PROGRAMS
!!
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
program gpmd

  !BML lib.
  use bml
  !Progress and LATTE lib modules.
  use progress_mod
  use system_mod
  use ptable_mod
  use latteparser_latte_mod
  use huckel_latte_mod
  use tbparams_latte_mod 
  use ham_latte_mod
  use coulomb_latte_mod
  use charges_mod
  use initmatrices_mod
  use genz_mod
  use nonortho_mod
  use pulaymixer_mod
  use dos_mod
  use density_mod
  use neighborlist_latte_mod
  use ppot_latte_mod  
  use hsderivative_latte_mod
  use slaterkosterforce_latte_mod
  use pulaycomponent_mod
  use nonorthocoulombforces_latte_mod
  ! Graph partitioning modules
  use parallel_mod
  use timer_mod
  use graph_sp2parser_mod
  use sp2parser_mod
  use sp2_mod
  use graph_mod
  use subgraphLoop_mod
  use homolumo_mod
  use xlbo_mod
  use md_latte_mod
  use partition_mod

  implicit none     

  integer, parameter                ::  dp = kind(1.0d0)
  integer                           ::  seed = 1
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable         ::  intKind(:)
  character(50)                     ::  inputfile, filename
  character(2)                      ::  auxchar
  integer                           ::  MD_step, Nr_SCF_It, i, icount
  integer                           ::  j, nel, norb, pp(100), nnodes
  integer                           ::  nparts, niter=500, npat, ipt
  integer                           ::  ii, jj, iscf
  integer, allocatable              ::  hindex(:,:), hnode(:)
  integer, allocatable              ::  xadj(:), adjncy(:), CH_count(:)  
  integer, allocatable              ::  part(:), core_count(:), Halo_count(:,:)  
  real(dp)                          ::  C0, C1, C2, C3
  real(dp)                          ::  C4, C5, ECoul, EKIN
  real(dp)                          ::  EPOT, ERep, Energy, Etot
  real(dp)                          ::  F2V, KE2T, MVV2KE, M_init
  real(dp)                          ::  TRRHOH, Temp, Time, alpha
  real(dp)                          ::  bndfil, cc, coulcut, dt
  real(dp)                          ::  dx, egap, ehomo, elumo
  real(dp)                          ::  kappa, scferror, traceMult, vv(100)
  real(dp)                          ::  sumCubes, maxCH, Ef, smooth_maxCH, pnorm=6
  real(dp)                          ::  dvdw, d
  real(dp), allocatable             ::  FPUL(:,:), FSCOUL(:,:), FTOT(:,:), PairForces(:,:)
  real(dp), allocatable             ::  SKForce(:,:), VX(:), VY(:), VZ(:)
  real(dp), allocatable             ::  charges_old(:), coul_forces(:,:), coul_forces_k(:,:), coul_forces_r(:,:)
  real(dp), allocatable             ::  coul_pot_k(:), coul_pot_r(:), dqin(:,:), dqout(:,:)
  real(dp), allocatable             ::  eigenvals(:), gbnd(:), n(:), n_0(:)
  real(dp), allocatable             ::  n_1(:), n_2(:), n_3(:), n_4(:)
  real(dp), allocatable             ::  n_5(:), onsitesH(:,:), onsitesS(:,:), rhoat(:)
  real(dp), allocatable             ::  origin(:)
  type(bml_matrix_t)                ::  aux_bml, dH0x_bml, dH0y_bml, dH0z_bml
  type(bml_matrix_t)                ::  dSx_bml, dSy_bml, dSz_bml, eigenvects
  type(bml_matrix_t)                ::  g_bml, ham0_bml, ham_bml, orthoh_bml
  type(bml_matrix_t)                ::  orthop_bml, over_bml, rho_bml, rhoat_bml
  type(bml_matrix_t)                ::  rhoh_bml, zmat_bml
  type(bml_matrix_t)                ::  copy_g_bml, gcov_bml
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
  logical                           :: tZSP
  type(genZSPinp)                   :: zsp

  integer, allocatable              :: xadj_cov(:), adjncy_cov(:), CH_count_cov(:)  
  integer, allocatable              :: part_cov(:), core_count_cov(:), Halo_count_cov(:,:)  
  integer                           :: vsize(2)
  integer                           :: nparts_cov


  !!!!!!!!!!!!!!!!!!!!!!!!
  !> Main program driver
  !!!!!!!!!!!!!!!!!!!!!!!!

  call gpmd_Init()

  call gpmd_Part()

  ! Loop over all the CH parts
  Ef = -2.5_dp
  do ipt = 1,gpat%TotalParts

    call gpmd_InitPart() 

    call gpmd_FirstCharges()

  enddo

  call gpmd_Gather()

  call gpmd_SCFLoop(lt%maxscf,sy%net_charge)

  call gpmd_Finalize()

contains

  !> initialize gpmd
  !!
  subroutine gpmd_Init()
    implicit none

    !> Start progress
    call progress_init()

    if (printRank() .eq. 1) then
      write(*,*) "GPMD started ..."
    endif

    call getarg(1, inputfile)
    write(*,*)""
    write(*,*)"Reading ",inputfile,"..." 
    write(*,*)""

    !> Parsing input file. This file contains all the variables needed to 
    !  run the scf including the sp2 (solver) variables. lt is "latte_type" structure 
    !  containing all the variables. 
    !  file://~/progress/build/doc/html/structlatteparser__latte__mod_1_1latte__type.html
    call parse_latte(lt,inputfile) 

    !> Parsing SP2 input paramenters. This will read the variables in the input file.
    !  sp2 is the "sp2data_type".     
    call parse_sp2(sp2,inputfile) 

    !> Parsing GSP2 input paramenters. This will read the variables in the input file. 
    !  gsp2 is the "gsp2data_type". 
    call parse_gsp2(gsp2,inputfile) 

    !> Parsing Extended Lagrangian input paramenters. This will read the variables in the input file. 
    !  xl is the "xlbo_type". 
    call parse_xlbo(xl,inputfile)

    !> Parsing Z sparse propagation. 
    !  note: no need to pass a structure. 
    ! call genZSPsolver%parse(zsp,inputfile)
    call parse_zsp(zsp,inputfile)

    !> Parsing system coordinates. This reads the coords.pdb file to get the position of every 
    !  atom in the system. sy is the "system_type" structure containing all the variables.
    !  file://~/progress/build/doc/html/classsystem__latte__mod.html
    call parse_system(sy,lt%coordsfile) 

    !> Building the neighbor list.
    call get_coulcut(lt%coul_acc,lt%timeratio,sy%nats,sy%lattice_vector,coulcut)
    call build_nlist(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)

    !> Get Huckel hamiltonian. Computes the Extended Huckel Hamiltonian from the 
    !  atom coordinates. The main inputs are the huckelTBparams and the system coordinate (sy%coordinate)
    !  The main outputs are Hamiltonian (ham_bml) and Overlap (over_bml) matrices.
    !   call get_hshuckel(ham_bml,over_bml,sy%coordinate,sy%spindex,sy%spatnum,&
    !     "../../huckelTBparams",lt%bml_type,lt%mdim,lt%threshold&
    !     ,tb%nsp,tb%splist,tb%basis,tb%numel,tb%onsite_energ,&
    !     tb%norbi,tb%hubbardu)    

    !> LATTE Hamiltonian parameter 
    call load_latteTBparams(tb,sy%splist,lt%parampath) 

    !> Get the reciprocal vectors
    call get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)

    call load_bintTBparamsH(sy%splist,tb%onsite_energ,&
      typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,lt%parampath) 
    call write_bintTBparamsH(typeA,typeB,&
      intKind,intPairsH,intPairsS,"mybondints.nonorth")  

    !> Load Pair potentials for latte TB.
    call load_PairPotTBparams(lt%parampath,sy%splist,ppot)

    !> Allocate bounds vector.
    allocate(gbnd(2))

  end subroutine gpmd_Init


  !> Partition by systems
  !!
  subroutine gpmd_Part
    implicit none

    call get_covgraph(sy,nl%nndist,nl%nnStruct,nl%nrnnstruct,gsp2%bml_type,3.0_dp,gcov_bml,gsp2%mdim)
    ! call initGraphPartitioning(gpat, "syprt", npat, sy%nats)
    ! call equalPartition(gpat, 400, sy%nats)

    call bml_print_matrix("gcov",gcov_bml,1,20,1,20)

    ! call bml_write_matrix(gcov_bml,"gcov.mtx")

    nparts_cov = 2

    allocate(xadj_cov(sy%nats+1))
    allocate(adjncy_cov(sy%nats*sy%nats))
    allocate(part_cov(sy%nats))
    allocate(core_count_cov(nparts_cov))
    allocate(CH_count_cov(nparts_cov))
    allocate(Halo_count_cov(sy%nats, sy%nats))

    call bml_adjacency(gcov_bml, xadj_cov, adjncy_cov, 1)

    call metisPartition(gpat, sy%nats, sy%estr%norbs, xadj_cov, adjncy_cov, nparts_cov, part_cov, core_count_cov,&
      CH_count_cov, Halo_count_cov, sumCubes, maxCH, smooth_maxCH, pnorm)                

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! WARNING: This refinements are not updating the gp core_pos structure !!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call simAnnealing(gpat, xadj_cov, adjncy_cov, part_cov, core_count_cov, CH_count_cov, &
    !       Halo_count_cov, sumCubes, maxCH, smooth_maxCH, pnorm, niter, seed)

    ! call KernLin2(gpat, xadj_cov, adjncy_cov, part_cov, core_count_cov, CH_count_cov, &
    !     Halo_count_cov, sumCubes, maxCH, smooth_maxCH, pnorm)

    call bml_copy_new(gcov_bml,aux_bml)
    do i=1,gpat%TotalParts
      call bml_matrix2submatrix_index(gcov_bml,&
        gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
        gpat%sgraph(i)%core_halo_index, &
        vsize,.true.)
      gpat%sgraph(i)%lsize = vsize(1)
      gpat%sgraph(i)%llsize = vsize(2)
    enddo
    call bml_deallocate(aux_bml)

    call printGraphPartitioning(gpat)

    allocate(syprt(gpat%nparts))    
    allocate(sy%resindex(sy%nats)) 
    sy%resindex(sy%nats)=-100

    !> For every partition get the partial CH systems.
    write(*,*)"Getting CH subsystems ..."
    do i=1,gpat%TotalParts
      call get_subsystem(sy,gpat%sgraph(i)%lsize,gpat%sgraph(i)%core_halo_index,syprt(i))      
      write(filename,*)i
      auxchar = adjustl(trim(filename))
      filename = "part_"//auxchar
      call write_system(syprt(i),filename,"pdb")      
      write(*,*)"SBY",i,gpat%nnodesInPart(i),gpat%sgraph(i)%llsize  
      ! do j=1,gpat%nnodesInPart(i)
      do j=1,gpat%sgraph(i)%llsize
        ! sy%resindex(gpat%sgraph(i)%nodeInPart(j)+1) = i
        ! sy%resindex(gpat%sgraph(i)%core_halo_index(gpat%sgraph(i)%core_pos(j))+1) = i
        sy%resindex(gpat%sgraph(i)%core_halo_index(j)+1) = i
        ! sy%resindex(gpat%sgraph(i)%core_pos(j)+1) = i
      enddo  
    enddo
    ! stop
    call write_system(sy,"system_parts","pdb")      

  end subroutine gpmd_Part


  !> Initialize the system parts.
  !!
  subroutine gpmd_InitPart

    write(*,*)""
    write(*,*)"         #######################"
    write(*,*)"           Computing partition =",ipt
    write(*,*)"         #######################"
    write(*,*)""

    !> Get the mapping of the Hamiltonian index with the atom index 
    !  hindex(1,i)=starting Hindex for atom i.
    !  hindex(2,i)=final Hindex for atom i.
    !  file://~/progress/build/doc/html/ham__latte__mod_8F90_source.html
    ! deallocate(hindex)
    if(allocated(hindex))deallocate(hindex)
    allocate(hindex(2,syprt(ipt)%nats))    
    allocate(syprt(ipt)%estr%hindex(2,syprt(ipt)%nats))    

    call get_hindex(syprt(ipt)%spindex,tb%norbi,hindex,norb)            
    syprt(ipt)%estr%hindex = hindex

    if(bml_get_N(ham_bml).GT.0)then
      call bml_deallocate(ham_bml)
      call bml_deallocate(over_bml)
    endif      
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ham_bml)    
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,over_bml)    

    call get_hsmat(ham_bml,over_bml,syprt(ipt)%coordinate,syprt(ipt)%lattice_vector,syprt(ipt)%spindex,&
      tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)

    nnodes = size(hindex, dim=2)

    if(lt%verbose.GE.5)then 
      if (printRank() .eq. 1) then
        write(*,*)
        call bml_print_matrix("H0",ham_bml,0,6,0,6)
        call bml_print_matrix("S",over_bml,0,6,0,6)
      endif
    endif

    !> Get occupation based on last shell population. 
    !  WARNING: This could change depending on the TB method being used.
    nel = sum(element_numel(syprt(ipt)%atomic_number(:)),&
      size(syprt(ipt)%atomic_number,dim=1))
    bndfil = nel/(2.0_dp*norb)

    if (printRank() .eq. 1) then
      write(*,*) "Filling fraction = ", bndfil
      write(*,*) "Number of electrons = ", nel
    endif

    call bml_deallocate(rhoat_bml)
    call build_atomic_density(rhoat_bml,tb%numel,hindex,syprt(ipt)%spindex,norb,&
      lt%bml_type)

    if(lt%verbose.GE.2) call bml_print_matrix("rhoat_bml",rhoat_bml,0,6,0,6)

    if(lt%zmat.EQ."ZSP")then 
      call init_ZSPmat(igenz,zk1_bml,zk2_bml,zk3_bml&
        ,zk4_bml,zk5_bml,zk6_bml,norb,lt%bml_type)
    endif

  end subroutine gpmd_InitPart


  !> First Charge computation.
  !! \brief Here we compute the non-scf charges.
  !!  
  subroutine gpmd_FirstCharges()
    implicit none

    !> Initialize the density matrix (rho_bml) and inverse overlap factor (zmat_bml).
    ! call init_pzmat(rho_bml,zmat_bml,lt%bml_type,lt%mdim,norb)
    if(bml_get_N(rho_bml).GT.0)then
      call bml_deallocate(rho_bml)
      call bml_deallocate(zmat_bml)
    endif  

    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml)                 
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zmat_bml)      

    !> Get the Inverse square root overlap matrix.
    call timer_start(dyn_timer,"Build Z")
    ! call buildzdiag(over_bml,zmat_bml,lt%threshold,lt%mdim,lt%bml_type)
    call gpmd_buildz()
    call timer_stop(dyn_timer,1)

    !> Initialize the orthogonal versions of ham and rho.
    ! call init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)
    if(bml_get_N(orthop_bml).GT.0)then
      call bml_deallocate(orthoh_bml)
      call bml_deallocate(orthop_bml)
    endif  
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthoh_bml)                 
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthop_bml)      

    call bml_print_matrix("H",ham_bml,0,6,0,6)

    if(lt%verbose.GE.5)then 
      if (printRank() .eq. 1) then
        write(*,*)
        call bml_print_matrix("Z",zmat_bml,0,6,0,6)
      endif
    endif

    !> Orthogonalize ham.
    call timer_start(ortho_timer)
    call orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
      lt%threshold,lt%bml_type,lt%verbose)
    call timer_stop(ortho_timer)

    call gpmd_rhosolver()

    !> Deorthogonalize rho.       
    call timer_start(deortho_timer)
    call deorthogonalize(orthop_bml,zmat_bml,rho_bml,&
      lt%threshold,lt%bml_type,lt%verbose)
    call timer_stop(deortho_timer)

    if(lt%verbose.GE.2)then
      if (printRank() .eq. 1) then
        call bml_print_matrix("rho0",rho_bml,0,6,0,6)       
      endif
    endif

    call bml_write_matrix(rho_bml,"rho.mtx")
    !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing 
    !  the charges.
    call prg_get_charges(rho_bml, over_bml, hindex, syprt(ipt)%net_charge, tb%numel, syprt(ipt)%spindex, lt%mdim, lt%threshold)

    if (printRank() .eq. 1) then
      write(*,*)""
      write(*,*)"Total charge =", sum(syprt(ipt)%net_charge(:),size(syprt(ipt)%net_charge,dim=1))
    endif

    if(lt%verbose.GE.2)then 
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"System charges:"
        do j=1,syprt(ipt)%nats
          write(*,*)syprt(ipt)%symbol(j),syprt(ipt)%net_charge(j)
        enddo
      endif
    endif

    syprt(ipt)%estr%norbs = norb
    call bml_copy_new(ham_bml,syprt(ipt)%estr%ham) !H0
    call bml_copy_new(rho_bml,syprt(ipt)%estr%rho)  !Rho0
    call bml_copy_new(over_bml,syprt(ipt)%estr%over) !S
    call bml_copy_new(zmat_bml,syprt(ipt)%estr%zmat) !Z

    ! write(filename,*)ipt
    ! auxchar = adjustl(trim(filename))
    ! filename = "dos_part_"//auxchar
    ! call get_eigenvalues(ham_bml,eigenvals,10)
    ! call write_tdos(eigenvals, 0.1_dp, 1000, -40.0_dp, 20.0_dp, filename)
    ! deallocate(eigenvals)

  end subroutine gpmd_FirstCharges


  !> Collect full system charges.
  !!
  subroutine gpmd_Gather()
    implicit none

    if(.not.allocated(sy%net_charge))allocate(sy%net_charge(sy%nats))
    if(.not.allocated(charges_old))allocate(charges_old(sy%nats))

    do i=1,gpat%TotalParts
      ! do j=1,gpat%sgraph(i)%llsize          
      do j=1,gpat%sgraph(i)%llsize          
        ! ii = gpat%sgraph(i)%core_pos(j)
        ! jj = gpat%sgraph(i)%core_halo_index(ii)
        jj = gpat%sgraph(i)%core_halo_index(j)+1
        sy%net_charge(jj) = syprt(i)%net_charge(j)
      enddo
      write(filename,*)i
      auxchar = adjustl(trim(filename))
      filename = "charge_"//auxchar
      call write_system(syprt(i),filename,"pdb")      
    enddo

    call write_system(sy,"charged_system","pdb")      

    write(*,*)""; write(*,*)"System charges:"
    do j=1,sy%nats
      write(*,*)j,sy%symbol(j),sy%net_charge(j)
    enddo

    charges_old = sy%net_charge

  end subroutine gpmd_Gather    


  !>  SCF loop.
  !!
  subroutine gpmd_SCFLoop(Nr_SCF,nguess)
    implicit none
    integer :: Nr_SCF
    real(dp):: nguess(:)

    sy%net_charge = nguess
    converged = .false.

    !> Beginning of the SCF loop.
    do iscf=1,Nr_SCF    

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This is done for the whole system 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"SCF iter", iscf
      endif

      !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"In real Coul ..."
      endif

      ! call get_ewald_real(sy%spindex,sy%splist,sy%coordinate&
      !   ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
      !   sy%volr,lt%coul_acc,coul_forces_r,coul_pot_r);
      call timer_start(dyn_timer,"Real coul")
      call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate&
        ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
        sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
        nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call timer_stop(dyn_timer,1)

      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"In recip Coul ..."    
      endif

      call timer_start(dyn_timer,"Recip coul")
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
        ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
        sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);  
      call timer_stop(dyn_timer,1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over the parts
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do ipt=1,gpat%TotalParts

        norb = syprt(ipt)%estr%norbs

        if(.not.allocated(syprt(ipt)%estr%coul_pot_k))then 
          allocate(syprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
          allocate(syprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))
        endif  

        !> Get Coulombic potential for the part and charges
        do j=1,gpat%sgraph(ipt)%lsize          
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          syprt(ipt)%estr%coul_pot_k(j) = coul_pot_k(jj)
          syprt(ipt)%estr%coul_pot_r(j) = coul_pot_r(jj)              
          syprt(ipt)%net_charge(j) = sy%net_charge(jj)              
        enddo

        if(bml_get_N(ham_bml).GT.0)then
          call bml_deallocate(ham_bml)
        endif      
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ham_bml)    

        !> Save actual hamiltonian as the non-scf Hamiltonian (H0)
        call bml_copy_new(syprt(ipt)%estr%ham,ham_bml)

        !> Get the scf hamiltonian. The outputs is ham_bml.
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"In get_hscf ..."
        endif
        call get_hscf(syprt(ipt)%estr%ham,syprt(ipt)%estr%over,ham_bml,syprt(ipt)%spindex,syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
          syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)

        !> Initialize the orthogonal versions of ham and rho.
        ! call init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)
        if(bml_get_N(orthop_bml).GT.0)then
          call bml_deallocate(orthoh_bml)
          call bml_deallocate(orthop_bml)
        endif  
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthoh_bml)                 
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthop_bml)      

        !> Orthogonalize ham.
        call timer_start(ortho_timer)
        call orthogonalize(ham_bml,syprt(ipt)%estr%zmat,orthoh_bml,&
          lt%threshold,lt%bml_type,lt%verbose)
        call timer_stop(ortho_timer)

        if(converged) call bml_copy(ham_bml,syprt(ipt)%estr%ham)

        call bml_deallocate(ham_bml)

        !> Now sove for the density matrix.
        call gpmd_rhosolver()
        call bml_deallocate(orthoh_bml)

        if(bml_get_N(rho_bml).GT.0)then
          call bml_deallocate(rho_bml)
        endif  
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml)                 

        !> Deorthogonalize orthop_bml to get the density matrix rho_bml.
        call timer_start(deortho_timer)
        call deorthogonalize(orthop_bml,syprt(ipt)%estr%zmat,rho_bml,&
          lt%threshold,lt%bml_type,lt%verbose)
        call timer_stop(deortho_timer)

        call bml_deallocate(orthop_bml)

        !> Get the system charges from rho
        call prg_get_charges(rho_bml,syprt(ipt)%estr%over,syprt(ipt)%estr%hindex,syprt(ipt)%net_charge,tb%numel,&
          syprt(ipt)%spindex,lt%mdim,lt%threshold)

        if(converged) call bml_copy(rho_bml,syprt(ipt)%estr%rho)

        call bml_deallocate(rho_bml)

        write(*,*)"Total charge of the part =", sum(syprt(ipt)%net_charge(:)),size(syprt(ipt)%net_charge,dim=1)

      enddo  

      do i=1,gpat%TotalParts
        do j=1,gpat%sgraph(i)%llsize          
          jj = gpat%sgraph(i)%core_halo_index(j)+1
          sy%net_charge(jj) = syprt(i)%net_charge(j)
        enddo
      enddo

      call qmixer(sy%net_charge,charges_old,dqin,&
        dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,lt%verbose)

      ! call linearmixer(sy%net_charge,charges_old,scferror,lt%mixcoeff,lt%verbose)

      charges_old = sy%net_charge

      if(lt%verbose.GE.1)then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"System charges:"
          do j=1,sy%nats
            write(*,*)sy%symbol(j),sy%net_charge(j)
          enddo
        endif
      endif

      if(converged)then ! To do a last extra step.
        exit
      else  
        if(scferror.lt.lt%scftol.and.iscf.gt.5) then 
          if (printRank() .eq. 1) then                
            write(*,*)""; write(*,*)"SCF converged within",iscf,"steps ..."
            write(*,*)"SCF error =",scferror
          endif
          converged = .true.
        endif 
      endif  

    enddo

    call write_system(sy,"charged_system","pdb")

    !> End of SCF loop.
    nguess = sy%net_charge

    if(lt%verbose.GE.1)then
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"System charges (After SCF):"
        do j=1,sy%nats
          write(*,*)sy%symbol(j),sy%net_charge(j)
        enddo
      endif
    endif  

  end subroutine gpmd_SCFLoop


  !> Solver for computing the density matrix. 
  !!
  subroutine gpmd_RhoSolver()
    implicit none

    call timer_start(dyn_timer,"Solver")

    if(lt%method.EQ."GSP2")then 
      call timer_start(graphsp2_timer)
      call subgraphSP2Loop(orthoh_bml, g_bml, orthop_bml, gp, lt%threshold)
      call timer_stop(graphsp2_timer)
      ! call sp2_alg1_seq(orthoh_bml,orthop_bml,lt%threshold, gp%pp, gp%maxIter, gp%vv)
    elseif(lt%method.EQ."SP2")then 
      call sp2_alg2(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
        ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
    elseif(lt%method.EQ."Diag")then  
      ! call build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
      ! call build_density_T(orthoh_bml,orthop_bml,lt%threshold,bndfil, 0.1_dp, Ef)
      call build_density_T_Fermi(orthoh_bml,orthop_bml,lt%threshold, 0.1_dp, Ef)      
      write(*,*)"ipt =",ipt,"Ef =",Ef
    else
      stop"No valid Method in LATTE parameters"
    endif

    call timer_stop(dyn_timer,1)

#ifdef DO_MPI_BLOCK
    call allGatherParallel(orthop_bml)
#endif

    if(lt%verbose.GE.2)then
      if (printRank() .eq. 1) then
        call bml_print_matrix("gsp2 orthop_bml",orthop_bml,0,6,0,6)
      endif
    endif  

  end subroutine gpmd_RhoSolver


  !> Graph partitioning subroutine.
  !!
  subroutine gpmd_GraphPart
    implicit none

    !> Symmetrize and Threshold the Matrix
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,copy_g_bml)
    call bml_threshold(g_bml, gsp2%gthreshold)
    call bml_transpose(g_bml, copy_g_bml)
    call bml_add(0.5_dp,g_bml,0.5_dp,copy_g_bml,0.0_dp)
    call bml_threshold(g_bml, gsp2%gthreshold)
    call bml_deallocate(copy_g_bml)
#ifdef DO_MPI_BLOCK
    call allGatherParallel(g_bml)
#endif

    !> Create graph partitioning - Use Block or METIS or METIS+SA or METIS+KL
    call timer_start(part_timer)

    !> Block partitioning
    if (gsp2%partition_type == "Block") then
      !> Partition by orbital or atom
      if (gsp2%graph_element == "Orbital") then
        call equalPartition(gp, gsp2%nodesPerPart, norb)
      else
        call equalGroupPartition(gp, hindex, nnodes, gsp2%nodesPerPart, norb)
      endif

      !> METIS, METIS+SA, or METIS+KL partitioning
    else

      !> Partition by orbital or atom
      if (gsp2%graph_element == "Orbital") then
        if(.not.allocated(xadj))allocate(xadj(norb+1))
        if(.not.allocated(adjncy))allocate(adjncy(norb*norb))
        call bml_adjacency(g_bml, xadj, adjncy, 1)
        nnodes = norb

      else
        write(*,*) "METIS with Atom graph elements is not available"
        stop
        !            allocate(xadj(nnodes+1))
        !            allocate(adjncy(nnodes*nnodes))
        !            call bml_adjacency_group(g_bml, hnode, nnodes, xadj, adjncy, 1)
      endif

#ifdef DO_GRAPHLIB
      nparts = gsp2%partition_count

      if (first_part) then
        allocate(part(nnodes))
        allocate(core_count(nparts))
      endif

      allocate(CH_count(nparts))
      allocate(Halo_count(nparts, nnodes))

      !> For METIS, if no refinement of first time, do full partitioning
      if (gsp2%partition_refinement == 'None' .or. first_part) then

        !> Which METIS partitioning
        select case(gsp2%partition_type)
        case("METIS")
          call metisPartition(gp, nnodes, norb, xadj, adjncy, nparts, part, core_count,&
            CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)                
        case("METIS+SA")
          call metisPartition(gp, nnodes, norb, xadj, adjncy, nparts, part, core_count,&
            CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm) 
          call simAnnealing(gp, xadj, adjncy, part, core_count, CH_count, &
            Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
        case("METIS+KL")
          call metisPartition(gp, nnodes, norb, xadj, adjncy, nparts, part, core_count,&
            CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm) 
          call KernLin2(gp, xadj, adjncy, part, core_count, CH_count, &
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
            call simAnnealing(gp, xadj, adjncy, part, core_count, CH_count, &
              Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
          case("KL")
            call KernLin2(gp, xadj, adjncy, part, core_count, CH_count, &
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

    call timer_stop(part_timer)

  end subroutine gpmd_GraphPart


  !> Construction of the congruence transformation.
  !!
  subroutine gpmd_BuildZ()
    implicit none

    igenz = igenz + 1

    write(*,*)zsp%nfirst  

    if(lt%zmat.eq."ZSP")then !Congruence transformation.

      call buildzsparse(over_bml,zmat_bml,igenz,lt%mdim,&
        lt%bml_type, zk1_bml,zk2_bml,zk3_bml&
        ,zk4_bml,zk5_bml,zk6_bml,zsp%nfirst,zsp%nrefi,zsp%nreff,&
        zsp%numthresi,zsp%numthresf,zsp%integration,zsp%verbose)

    else               

      !Build Z matrix using diagonalization (usual method).
      call buildzdiag(over_bml,zmat_bml,lt%threshold,lt%mdim,lt%bml_type)

    endif

  end subroutine gpmd_BuildZ


  !> Finalize gpmd 
  !!
  subroutine gpmd_Finalize()
    implicit none
    deallocate(gbnd)

    call bml_deallocate(ham_bml)
    call bml_deallocate(ham0_bml)
    call bml_deallocate(orthoh_bml)
    call bml_deallocate(g_bml)
    call bml_deallocate(rho_bml)
    call bml_deallocate(over_bml)
    call bml_deallocate(zmat_bml)

    ! Progress is done
    call progress_shutdown()

  end subroutine gpmd_Finalize

end program gpmd
