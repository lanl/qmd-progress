!> High-level program to perform SFC cycles in Extended Huckel Hamiltonian.
!! This program takes coordinates in xyz or pdb format and extracts information 
!! about the system. 
program gpscf_dist

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
  use densitymatrix_mod
  use neighborlist_latte_mod
  use ppot_latte_mod
  use hsderivative_latte_mod
  use slaterkosterforce_latte_mod
  use pulaycomponent_mod
  use nonorthocoulombforces_latte_mod
  use xlbo_mod
  use md_latte_mod
  use sp2parser_mod
  ! Graph partitioning mudules
  use parallel_mod
  use timer_mod
  use graph_sp2parser_mod
  use sp2_mod
  use graph_mod
  use partition_mod
  use subgraphLoop_mod
  use homolumo_mod
  use normalize_mod

  implicit none     

  integer, parameter     ::  dp = kind(1.0d0)
  integer                ::  seed = 1
  integer                ::  i, j, k, nel, norb, nnodes
  integer                ::  nparts, niter=500, npat, ipt
  integer                ::  ii, jj, iscf
  integer                ::  igenz
  integer, allocatable   ::  hindex(:,:), hnode(:), vsize(:)
  integer, allocatable   ::  xadj(:), adjncy(:), CH_count(:)
  integer, allocatable   ::  part(:), core_count(:), Halo_count(:,:)

  character(2), allocatable ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable ::  intKind(:)
  character(50)          ::  inputfile, filename

  real(dp)               ::  bndfil, scferror, coulcut
  real(dp)               ::  sumCubes, maxCH, Ef, smooth_maxCH, pnorm=6
  real(dp), allocatable  ::  charges_old(:), coul_forces_k(:,:), coul_forces_r(:,:), coul_pot_k(:)
  real(dp), allocatable  ::  n_5(:), onsitesH(:,:), onsitesS(:,:), rhoat(:)
  real(dp), allocatable  ::  coul_pot_r(:), dqin(:,:), dqout(:,:), eigenvals(:)
  real(dp), allocatable  ::  auxcharge(:), origin(:)

  type(bml_matrix_t)     ::  ham0_bml, ham_bml, orthoh_bml, orthop_bml
  type(bml_matrix_t)     ::  over_bml, rho_bml, zmat_bml, g_bml, eigenvects
  type(bml_matrix_t)     ::  gch_bml, rhoat_bml
  type(bml_matrix_t)     ::  zk1_bml, zk2_bml, zk3_bml
  type(bml_matrix_t)     ::  zk4_bml, zk5_bml, zk6_bml

  type(latte_type)       ::  lt
  type(sp2data_type)     ::  sp2
  type(gsp2data_type)    ::  graph2
  type(system_type)      ::  sy
  type(xlbo_type)        ::  xl
  type(tbparams_type)    ::  tb
  type(genZSPinp)        :: zsp

  type(neighlist_type)   ::  nl
  type(system_type), allocatable    ::  syprt(:)
  type(intpairs_type), allocatable  ::  intPairsH(:,:), intPairsS(:,:)
  type(ppot_type), allocatable      ::  ppot(:,:)

  type(graph_partitioning_t) :: gp

  integer :: icount, myRank
  integer :: pp(100)
  real(dp) :: vv(100)
  real(dp) :: ehomo, elumo, egap, traceMult
  real(dp), allocatable :: gbnd(:)

  logical :: first_part = .true.
  logical :: firstLocalPart

  !!!!!!!!!!!!!!!!!!!!!!!!
  !> Main program driver
  !!!!!!!!!!!!!!!!!!!!!!!!

  igenz = 0
  Ef = -2.5_dp

  !> Init
  call gpmd_init()
  
  !> First SP2 and charge calcualtion
  call gpmd_first()

  !> Graph partitioning
  call gpmd_part()

  !> SCF loop
  call gpmd_scf(lt%maxscf,sy%net_charge,.true.)

  !> Finalize
  call gpmd_final()
  
contains

  !> initialize gpmd
  !!
  subroutine gpmd_init()

    implicit none

    !> Start progress
    call progress_init()

    !> Get MPI rank
    myRank = getMyRank() + 1

    if (printRank() .eq. 1) then
      write(*,*) "gpdist start ..."
    endif

    call getarg(1, inputfile)
    if (printRank() .eq. 1) then
      write(*,*)""
      write(*,*)"Reading ",inputfile,"..."
      write(*,*)""
    endif

    !> Parsing input file. This file contains all the variables needed to 
    !  run the scf including the sp2 (solver) variables. lt is "latte_type"
    !  structure 
    !  containing all the variables. 
    !  file://~/progress/build/doc/html/structlatteparser__latte__mod_1_1latte__type.html
    call parse_latte(lt,inputfile)

    !> Parsing SP2 input paramenters. This will read the variables in the input
    !file.
    !  sp2 is the "sp2data_type".     
    call parse_sp2(sp2,inputfile)

    !> Parsing GSP2 input paramenters. This will read the variables in the input
    !file. 
    !  gsp2 is the "gsp2data_type". 
    call parse_gsp2(graph2,inputfile)

    !> Parsing Extended Lagrangian input paramenters. This will read the
    !variables in the input file. 
    !  xl is the "xlbo_type". 
    call parse_xlbo(xl,inputfile)

    !> Parsing Z sparse propagation. 
    !  note: no need to pass a structure. 
    ! call genZSPsolver%parse(zsp,inputfile)
    call parse_zsp(zsp,inputfile)

    !> Parsing system coordinates. This reads the coords.pdb file to get the
    !position of every 
    !  atom in the system. sy is the "system_type" structure containing all the
    !  variables.
    !  file://~/progress/build/doc/html/classsystem__latte__mod.html
    call parse_system(sy,lt%coordsfile)

    !> Center sytem inside the box and fold it by the lattice_vectors.
    call translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin)
    call write_system(sy,adjustl(trim(lt%jobname))//"_centered","pdb")

    !> Get the Coulombic cut off.
    call get_coulcut(lt%coul_acc,lt%timeratio,sy%nats,sy%lattice_vector,coulcut)

    !> Building the neighbor list.
    call build_nlist(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)

    !> LATTE Hamiltonian parameter 
    call load_latteTBparams(tb,sy%splist,lt%parampath)

    !> Allocate bounds vector.
    allocate(gbnd(2))

  end subroutine gpmd_init

  subroutine gpmd_first

    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>  First Charge computation 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Get the mapping of the Hamiltonian index with the atom index 
    !  hindex(1,i)=starting Hindex for atom i.
    !  hindex(2,i)=final Hindex for atom i.
    !  file://~/progress/build/doc/html/ham__latte__mod_8F90_source.html
    call get_hindex(sy%spindex,tb%norbi,hindex,norb)

    !> Bond integrals parameters for LATTE Hamiltonian.
    call load_bintTBparamsH(sy%splist,tb%onsite_energ,&
      typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,lt%parampath)
    call write_bintTBparamsH(typeA,typeB,&
      intKind,intPairsH,intPairsS,"mybondints.nonorth")

    !> Initialize hamiltonian matrix (ham_bml) and overlap matrix (over_bml)
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ham_bml, &
      lt%bml_dmode)
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,over_bml, &
      lt%bml_dmode)

    !> Create first full hamiltonian and overlap matrices
    call get_hsmat(ham_bml,over_bml,sy%coordinate,sy%lattice_vector,sy%spindex,&
      tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)

    nnodes = size(hindex, dim=2)
    allocate(hnode(nnodes))
    do i = 1, nnodes
      hnode(i) = hindex(1, i)
    enddo

    !> Get the reciprocal vectors
    call get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)

    if (printRank() .eq. 1) then
      write(*,*) "Number of orbitals = ", norb
      write(*,*)
    endif

    !> Get occupation based on last shell population. 
    !  WARNING: This could change depending on the TB method being used.
    nel = sum(element_numel(sy%atomic_number(:)),&
      size(sy%atomic_number,dim=1))
    bndfil = nel/(2.0_dp*norb)
    if (printRank() .eq. 1) then
      write(*,*) "bndfil = ", bndfil
      write(*,*) "nel = ", nel
    endif

    !> Load Pair potentials for latte TB.
    call load_PairPotTBparams(lt%parampath,sy%splist,ppot)

    call build_atomic_density(rhoat_bml,tb%numel,hindex,sy%spindex,norb,&
      lt%bml_type)

    call init_ZSPmat(igenz,zk1_bml,zk2_bml,zk3_bml&
      ,zk4_bml,zk5_bml,zk6_bml,norb,lt%bml_type)

    !> Initialize the inverse overlap factor (zmat_bml).
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zmat_bml, &
      lt%bml_dmode)

    !> Get the Inverse square root overlap matrix.
    ! over is distributed, but a full copy is on each node
    ! zmat will be sequential, same on each node
    call gpmd_buildz(over_bml, zmat_bml, lt%bml_type)

    !> Initialize the orthogonal version of ham.
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthoh_bml, &
      lt%bml_dmode)

    !> Orthogonalize ham.
    call timer_start(ortho_timer)
    call orthogonalize(ham_bml,zmat_bml,orthoh_bml,lt%threshold,lt%bml_type,lt%verbose)
    call timer_stop(ortho_timer)

    !> Calculate gershgorin bounds
    call bml_gershgorin(orthoh_bml, gbnd)
    if (printRank() .eq. 1) then
      write(*,*) "Gershgorin: mineval = ", gbnd(1), " maxeval = ", gbnd(2)
      write(*,*)
    endif

    !> SP2 algorithm.
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthop_bml, &
      lt%bml_dmode)
    call timer_start(sp2_timer)
    call sp2_alg2_genseq(orthoh_bml,orthop_bml,lt%threshold,bndfil,sp2%minsp2iter,&
      sp2%maxsp2iter,sp2%sp2conv,sp2%sp2tol, pp, icount, vv)
    call timer_stop(sp2_timer)

    if (printRank() .eq. 1) then
      write(*,*)
      write(*,*) "SP2 count = ", icount
      do i = 1, icount
        write(*,*)"rank = ", myRank, "i = ", i, " pp = ", pp(i), " vv = ", vv(i)
      enddo
      write(*,*)
    endif

    !> Calculate Homo-Lumo gap
    call homolumogap(vv, icount, pp, gbnd(1), gbnd(2), ehomo, elumo, egap)
    if (printRank() .eq. 1) then
      write(*,*) "rank = ", getMyRank(), "Homo-lumo: ehomo = ", ehomo, &
        " elumo = ", elumo, " egap = ", egap
      write(*,*)
    endif

#ifdef DO_MPI
    if (getNRanks() .gt. 1 .and. &
      bml_get_distribution_mode(orthop_bml) == BML_DMODE_DISTRIBUTED) then
      call allGatherParallel(orthop_bml)
    endif
#endif

    !> Make an atom-based version of the density matrix adjacency matrix
    ! for the first graph
    call bml_group_matrix(orthop_bml, g_bml, hnode, nnodes, graph2%gthreshold)

    !> Deorthogonalize rho.       
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml, &
      lt%bml_dmode)
    call timer_start(deortho_timer)
    call deorthogonalize(orthop_bml,zmat_bml,rho_bml,lt%threshold,lt%bml_type,lt%verbose)
    call timer_stop(deortho_timer)

#ifdef DO_MPI
    if (getNRanks() .gt. 1 .and. &
      bml_get_distribution_mode(rho_bml) == BML_DMODE_DISTRIBUTED) then
      call allGatherParallel(rho_bml)
    endif
#endif

    !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing 
    !  the charges.
    call prg_get_charges(rho_bml, over_bml, hindex, sy%net_charge, tb%numel, sy%spindex, lt%mdim, lt%threshold)
    charges_old = sy%net_charge

    if (printRank() .eq. 1) then
      write(*,*)"Rank = ", getMyRank(), " Total charges =", sum(sy%net_charge(:), &
        size(sy%net_charge,dim=1))
      write(*,*)
    endif

    call bml_deallocate(ham_bml)
    call bml_deallocate(rho_bml)
    call bml_deallocate(over_bml)
    call bml_deallocate(zmat_bml)
    call bml_deallocate(orthoh_bml)
    call bml_deallocate(orthop_bml)

  end subroutine gpmd_first

  subroutine gpmd_part
 
    implicit none

    !> Create graph partitioning - Use Block or METIS or METIS+SA or METIS+KL
    call timer_start(part_timer)

    !> Block partitioning
    if (graph2%partition_type == "Block") then
      !> Partition by orbital or atom
      if (graph2%graph_element == "Orbital") then
        call equalPartition(gp, graph2%nodesPerPart, norb)
      else
        call equalGroupPartition(gp, hindex, nnodes, graph2%nodesPerPart, &
          norb)
      endif
 
     !> METIS, METIS+SA, or METIS+KL partitioning
    else
!      !> Partition by orbital or atom
      if (graph2%graph_element == "Orbital") then
        allocate(xadj(norb+1))
        allocate(adjncy(norb*norb))
        call bml_adjacency(g_bml, xadj, adjncy, 1)
        nnodes = norb
      else
        allocate(xadj(nnodes+1))
        allocate(adjncy(nnodes*nnodes))
        call bml_adjacency(g_bml, xadj, adjncy, 1)
      endif

      write(*,*) "nnodes =", nnodes
      write(*,*) "norb =", norb

#ifdef DO_GRAPHLIB
      nparts = graph2%partition_count

      if (first_part) then
        allocate(part(nnodes))
        allocate(core_count(nparts))
      endif

      allocate(CH_count(nparts))
      allocate(Halo_count(nparts, nnodes))

      !> For METIS, if no refinement or first time, do full partitioning
      if (graph2%partition_refinement == 'None' .or. first_part) then

        !> Which METIS partitioning
        select case(graph2%partition_type)
          case("METIS")
            call metisPartition(gp, nnodes, nnodes, xadj, adjncy, nparts, part, &
              core_count, CH_count, Halo_count, sumCubes, maxCH, &
              smooth_maxCH, pnorm)           
          case("METIS+SA")
            call metisPartition(gp, nnodes, nnodes, xadj, adjncy, nparts, part, &
              core_count, CH_count, Halo_count, sumCubes, maxCH, &
              smooth_maxCH, pnorm)
            call simAnnealing(gp, xadj, adjncy, part, core_count, CH_count, &
              Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
          case("METIS+KL")
            call metisPartition(gp, nnodes, nnodes, xadj, adjncy, nparts, part, core_count,&
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
          select case(graph2%partition_refinement)
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

    !> Determine halo indeces for all partitions for graph
    call getPartitionHalosFromGraph(gp, g_bml, graph2%double_jump)

    !> Calculate part ordering and reordering
    call partOrdering(gp)

    call timer_stop(part_timer)

  end subroutine gpmd_part

  !> Initialize the system parts.
  !!
  subroutine gpmd_initPart(ipt)

    implicit none

    integer :: ipt

    !> For each partition get the partial systems.
    call get_subsystem(sy,gp%sgraph(ipt)%lsize,gp%sgraph(ipt)%core_halo_index, &
      syprt(ipt))

    !> Get the mapping of the Hamiltonian index with the atom index 
    !  hindex(1,i)=starting Hindex for atom i.
    !  hindex(2,i)=final Hindex for atom i.
    !  file://~/progress/build/doc/html/ham__latte__mod_8F90_source.html
    if(allocated(hindex))deallocate(hindex)
    allocate(hindex(2,syprt(ipt)%nats))
    if (allocated(syprt(ipt)%estr%hindex)) &
      deallocate(syprt(ipt)%estr%hindex)
    allocate(syprt(ipt)%estr%hindex(2,syprt(ipt)%nats))

    call get_hindex(syprt(ipt)%spindex,tb%norbi,hindex,norb)

    syprt(ipt)%estr%hindex = hindex
    syprt(ipt)%estr%norbs = norb

    if(bml_get_N(syprt(ipt)%estr%ham0).GT.0)then
      call bml_deallocate(syprt(ipt)%estr%ham0)
      call bml_deallocate(syprt(ipt)%estr%over)
    endif
    call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb, &
      syprt(ipt)%estr%ham0)
    call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb, &
      syprt(ipt)%estr%over)

    call get_hsmat(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over, &
      syprt(ipt)%coordinate,syprt(ipt)%lattice_vector,syprt(ipt)%spindex,&
      tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)

    nnodes = size(hindex, dim=2)

    !> Get occupation based on last shell population. 
    !  WARNING: This could change depending on the TB method being used.
    nel = sum(element_numel(syprt(ipt)%atomic_number(:)),&
      size(syprt(ipt)%atomic_number,dim=1))
    bndfil = nel/(2.0_dp*norb)

    if (printRank() .eq. 1) then
      write(*,*)
      write(*,*) "ipt = ", ipt
      write(*,*) "norb = ", norb
      write(*,*) "nnodes = ", nnodes
      write(*,*) "Filling fraction = ", bndfil
      write(*,*) "Number of electrons = ", nel
      write(*,*)
    endif

    call bml_deallocate(rhoat_bml)
    call build_atomic_density(rhoat_bml,tb%numel,hindex,syprt(ipt)%spindex,norb,&
      BML_MATRIX_DENSE)

!    if(lt%verbose.GE.2) call bml_print_matrix("rhoat_bml",rhoat_bml,0,6,0,6)

    !> Initialize the inverse overlap factor (zmat_bml).
    if(bml_get_N(syprt(ipt)%estr%zmat).GT.0) then 
      call bml_deallocate(syprt(ipt)%estr%zmat)
    endif
    call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb, &
      syprt(ipt)%estr%zmat)

    !> Get the Inverse square root overlap matrix.
    call timer_start(buildz_timer)
    call gpmd_buildz(syprt(ipt)%estr%over, syprt(ipt)%estr%zmat, &
      BML_MATRIX_DENSE)
    call timer_stop(buildz_timer)

  end subroutine gpmd_initPart

  !> SCF loop
  subroutine gpmd_scf(Nr_SCF,nguess,mix)

    implicit none

    integer :: Nr_SCF
    real(dp):: nguess(:)
    logical :: mix

    charges_old = nguess

    ! Make charge changes in aux
    if(.not.allocated(auxcharge))allocate(auxcharge(sy%nats))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>  SCF loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> For every partition get the partial systems.
    if (printRank() .eq. 1) write(*,*)"Getting subsystems ..."
    if(allocated(syprt))deallocate(syprt)
    allocate(syprt(gp%totalParts))
    if(allocated(sy%resindex))deallocate(sy%resindex)
    allocate(sy%resindex(sy%nats))
    sy%resindex(sy%nats)=-100

    !> Beginning of the SCF loop.
    do iscf=1, Nr_SCF   

      if (printRank() .eq. 1) then
        write(*,*)"SCF iter", iscf
      endif

      !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
      if (printRank() .eq. 1) then
        write(*,*)"In real Coul ..."
      endif
      call timer_start(realcoul_timer)
      call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate,&
        nguess,tb%hubbardu,sy%lattice_vector,&
        sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
        nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call timer_stop(realcoul_timer)

      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      if (printRank() .eq. 1) then
        write(*,*)"In recip Coul ..."    
      endif
      call timer_start(recipcoul_timer)
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
        ,nguess,tb%hubbardu,sy%lattice_vector,&
        sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
      call timer_stop(recipcoul_timer)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over the local parts
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      firstLocalPart = .true.
      auxcharge = 0.0_dp
      do ipt = gp%localPartMin(myRank), gp%localPartMax(myRank)

        call gpmd_initPart(ipt)

        norb = syprt(ipt)%estr%norbs

        if(.not.allocated(syprt(ipt)%estr%coul_pot_k))then
          allocate(syprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
          allocate(syprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))
        endif

        syprt(ipt)%estr%coul_pot_k = 0.0_dp
        syprt(ipt)%estr%coul_pot_r = 0.0_dp
        syprt(ipt)%net_charge = 0.0_dp

        !> Get Coulombic potential for the part and charges
        do j=1,gp%sgraph(ipt)%lsize
          jj = gp%sgraph(ipt)%core_halo_index(j)+1
          syprt(ipt)%estr%coul_pot_k(j) = coul_pot_k(jj)
          syprt(ipt)%estr%coul_pot_r(j) = coul_pot_r(jj)
          syprt(ipt)%net_charge(j) = nguess(jj)
        enddo

        if(bml_get_N(syprt(ipt)%estr%ham).GT.0)then
          call bml_deallocate(syprt(ipt)%estr%ham)
        endif
        call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb, &
          syprt(ipt)%estr%ham)

        !> Save actual hamiltonian as the non-scf Hamiltonian (H0)
        call bml_copy_new(syprt(ipt)%estr%ham0, syprt(ipt)%estr%ham)

        !> Get the scf hamiltonian.
        if (printRank() .eq. 1) then
          write(*,*)"In get_hscf ..."
        endif
        call get_hscf(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over, &
          syprt(ipt)%estr%ham, syprt(ipt)%spindex,syprt(ipt)%estr%hindex, &
          tb%hubbardu, syprt(ipt)%net_charge,syprt(ipt)%estr%coul_pot_r, &
          syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)

        !> Orthogonalize the Hamiltonian
        if (printRank() .eq. 1) then
          write(*,*)"in orthogonalize H ..."
        endif
        !> Initialize the orthogonal versions of ham and rho.
        if(bml_get_N(syprt(ipt)%estr%oham).GT.0)then
          call bml_deallocate(syprt(ipt)%estr%oham)
        endif
        call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb, &
          syprt(ipt)%estr%oham)         

        !> Orthogonalize ham.
        call timer_start(ortho_timer)
        call orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat, &
          syprt(ipt)%estr%oham, lt%threshold,BML_MATRIX_DENSE,lt%verbose)
        call timer_stop(ortho_timer)

        !> Calculate local part gershgorin bounds
        call bml_gershgorin_partial(syprt(ipt)%estr%oham, gbnd, gp%sgraph(ipt)%llsize)
        if (firstLocalPart) then
          gp%mineval = gbnd(1)
          gp%maxeval = gbnd(2) 
          firstLocalPart = .false.
        else
          if (gbnd(1) .lt. gp%mineval) then
            gp%mineval = gbnd(1)
          endif
          if (gbnd(2) .gt. gp%maxeval) then
            gp%maxeval = gbnd(2)
          endif
        endif
       
      enddo

      !if (printRank() .eq. 1) then
        write(*,*) "Rank", myRank-1, "Local Gershgorin Bounds:", &
          " mineval = ", gp%mineval, " maxeval = ", gp%maxeval
        write(*,*)
      !endif

      !> Do a reduction of gershgorin bounds across all ranks.
      call gershgorinReduction(gp)

      if (printRank() .eq. 1) then
        write(*,*) "Gershgorin Bounds:", &
          " mineval = ", gp%mineval, " maxeval = ", gp%maxeval
        write(*,*) 
      endif

      !! Calculate SP2 sequence
      call sp2sequence(gp%pp, gp%maxIter, gp%mineval, gp%maxeval, ehomo, elumo, &
        graph2%errlimit)
      if (printRank() .eq. 1) then
        write(*,*)
        write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter
        write(*,*)
      endif

    gp%vv = 0.0_dp
    do ipt = gp%localPartMin(myRank), gp%localPartMax(myRank)

      if(bml_get_N(syprt(ipt)%estr%orho).GT.0)then
        call bml_deallocate(syprt(ipt)%estr%orho)
      endif
      call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb,&
        syprt(ipt)%estr%orho)

      !> Solve for the density matrix
      call gpmd_rhoSolver(ipt, syprt(ipt)%estr%oham, &
        syprt(ipt)%estr%orho, gp, gp%sgraph(ipt)%llsize)
      call bml_deallocate(syprt(ipt)%estr%oham)

      if(bml_get_N(rho_bml).GT.0)then
        call bml_deallocate(syprt(ipt)%estr%rho)
      endif
      call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norb,norb, &
        syprt(ipt)%estr%rho)         

      !> Deorthogonalize orthop_bml to get the density matrix rho_bml.
      call timer_start(deortho_timer)
      call deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat, &
        syprt(ipt)%estr%rho, lt%threshold,BML_MATRIX_DENSE,lt%verbose)
      call timer_stop(deortho_timer)

      !> Get the system charges from rho
      call prg_get_charges(syprt(ipt)%estr%rho,syprt(ipt)%estr%over, &
        syprt(ipt)%estr%hindex, syprt(ipt)%net_charge,tb%numel, &
        syprt(ipt)%spindex,lt%mdim,lt%threshold)

      write(*,*)"Total charge of the part ", ipt, "  = ", &
        sum(syprt(ipt)%net_charge(:)),size(syprt(ipt)%net_charge,dim=1)

      do j=1,gp%sgraph(ipt)%llsize
          jj = gp%sgraph(ipt)%core_halo_index(j)+1
          auxcharge(jj) = syprt(ipt)%net_charge(j)
      enddo

    enddo

    !> Sum trace norm from all parts
    call fnormGraph(gp)

    !> Calculate Homo-Lumo gap
    call homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, egap)

    if (printRank() .eq. 1) then
      write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
        " egap = ", egap
      write(*,*)
    endif

    !! Calculate Trace[HP]
!    traceMult = bml_traceMult(orthoh_bml, orthop_bml)
!#ifdef DO_MPI_BLOCK
!    if (getNRanks() > 1) then
!      call sumRealReduce(traceMult)
!    endif
!#endif
!    if (printRank() .eq. 1) then
!      write(*,*) "Trace[HP] for SP2 = ", traceMult
!      write(*,*) "Band energy per atom = ", traceMult/gsp2%natoms
!      write(*,*)
!    endif

    !> Assemble density matrix parts on a local node
    !! Threshold for next graph.
    do ipt = gp%localPartMin(myRank), gp%localPartMax(myRank)

      call get_partial_atomgraph(syprt(ipt)%estr%orho, &
        syprt(ipt)%estr%hindex,gch_bml,graph2%gthreshold)
      call bml_submatrix2matrix(gch_bml, g_bml, &
        gp%sgraph(ipt)%core_halo_index, gp%sgraph(ipt)%lsize, &
        gp%sgraph(ipt)%llsize, graph2%gthreshold)

    enddo

    !> Collect matrix together from distributed sets of parts
    call collectMatrixFromParts(gp, g_bml)

    !> Do new partitioning
    call destroyGraphPartitioning(gp)
    call gpmd_part()

    !> Gather charges from all nodes
    call sumRealReduceN(auxcharge, gp%totalNodes)
    write(*,*)"Total system charge = ", sum(auxcharge(:))

    nguess = auxcharge

    !> Pulay mixing
    if (mix) then
      call qmixer(nguess,charges_old,dqin,&
        dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,lt%verbose)
    endif

    charges_old = nguess
    
    if(scferror.lt.lt%scftol.and.iscf.gt.5) then
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"SCF converged within",iscf,"steps ..."
        write(*,*)"SCF error =",scferror
      endif
      exit
    endif

  enddo

  call write_system(sy,"charged_system","pdb")

  !> End of SCF loop.

  deallocate(auxcharge)

  if (printRank() .eq. 1) then
    write(*,*)""
    write(*,*)"Total charge of the system (After SCF) =", &
      sum(nguess(:),size(nguess,dim=1))
    write(*,*)""
  endif

  if(lt%verbose.GE.1)then
    if (printRank() .eq. 1) then
      write(*,*)""; write(*,*)"System charges (After SCF):"
      do j=1,sy%nats
        write(*,*)sy%symbol(j),nguess(j)
      enddo
    endif
  endif

  end subroutine gpmd_scf


  !> Solver for computing the density matrix. 
  !!
  subroutine gpmd_rhoSolver(ipt, orthoh_bml, orthop_bml, gp, core_size)

    implicit none

    type (bml_matrix_t), intent(in) :: orthoh_bml
    type (bml_matrix_t), intent(inout) :: orthop_bml
    type (graph_partitioning_t), intent(inout) :: gp
    integer, intent(in) :: ipt, core_size

    real(dp) :: thresh0

    thresh0 = 0.0_dp

    ! Subgraph SP2
    if(lt%method .EQ. "GSP2")then

      call timer_start(graphsp2_timer)
      call sp2_submatrix(orthoh_bml, orthop_bml, thresh0, gp%pp, &
        gp%maxIter, gp%vv, gp%mineval, gp%maxeval, &
        core_size)
      call timer_stop(graphsp2_timer)

    ! Diagonalization
    elseif(lt%method .EQ. "Diag")then

      call build_density_T_Fermi(orthoh_bml,orthop_bml,lt%threshold, 0.1_dp, Ef)
      write(*,*)"ipt =",ipt,"Ef =",Ef

    else
      stop"No valid Method in LATTE parameters"
    endif

    if(lt%verbose.GE.2)then
      if (printRank() .eq. 1) then
        call bml_print_matrix("orthop_bml",orthop_bml,0,6,0,6)
      endif
    endif

  end subroutine gpmd_rhoSolver

  subroutine gpmd_buildz(over_bml, zmat_bml, bml_type)

    implicit none

    type (bml_matrix_t), intent(inout) :: over_bml
    type (bml_matrix_t), intent(inout) :: zmat_bml
    character(len=*), intent(in) :: bml_type
    
    igenz = igenz + 1

    if(lt%zmat.eq."ZSP")then !Congruence transformation.

      call buildzsparse(over_bml,zmat_bml,igenz,lt%mdim,&
        lt%bml_type, zk1_bml,zk2_bml,zk3_bml&
        ,zk4_bml,zk5_bml,zk6_bml,zsp%nfirst,zsp%nrefi,zsp%nreff,&
        zsp%numthresi,zsp%numthresf,zsp%integration,zsp%verbose)

    else

      !Build Z matrix using diagonalization (usual method).
      call buildzdiag(over_bml,zmat_bml,lt%threshold,lt%mdim,bml_type)

    endif

  end subroutine gpmd_buildz

  subroutine gpmd_final

    implicit none

    if (allocated(part) .eqv. .true.) deallocate(part)
    if (allocated(core_count) .eqv. .true.) deallocate(core_count)

    deallocate(gbnd)
    deallocate(hnode)

    call bml_deallocate(ham_bml)
    call bml_deallocate(ham0_bml)
    call bml_deallocate(orthoh_bml)
    call bml_deallocate(g_bml)
    call bml_deallocate(rho_bml)
    call bml_deallocate(over_bml)
    call bml_deallocate(zmat_bml)

    !> Progress is done
    call progress_shutdown()

  end subroutine gpmd_final

end program gpscf_dist
