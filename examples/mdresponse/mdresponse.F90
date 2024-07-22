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
program mdresponse

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
  use prg_sp2parser_mod
  use prg_pulaycomponent_mod
  use nonorthocoulombforces_latte_mod
  use prg_xlbo_mod
  ! Graph partitioning modules
  use prg_parallel_mod
  use prg_timer_mod
  use prg_graphsp2parser_mod
  use prg_sp2_mod
  use prg_graph_mod
  use prg_subgraphLoop_mod
  use prg_homolumo_mod
  use md_latte_mod
  use prg_partition_mod

  implicit none

  integer, parameter                ::  dp = kind(1.0d0)
  integer                           ::  seed = 1
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable         ::  intKind(:)
  character(50)                     ::  inputfile, filename
  character(2)                      ::  auxchar
  integer                           ::  MD_step, Nr_SCF_It, i, icount
  integer                           ::  j, nel, norb, pp(100), nnodes, count
  integer                           ::  nparts, niter=500, npat
  integer, allocatable              ::  hindex(:,:), hnode(:)
  integer, allocatable              ::  xadj(:), adjncy(:), CH_count(:)
  integer, allocatable              ::  part(:), core_count(:), Halo_count(:,:)
  real(dp)                          ::  C0, C1, C2, C3
  real(dp)                          ::  C4, C5, ECoul, EKIN
  real(dp)                          ::  EPOT, ERep, Energy, Etot
  real(dp)                          ::  F2V, KE2T, MVV2KE, M_prg_init
  real(dp)                          ::  TRRHOH, Temp, Time, alpha
  real(dp)                          ::  bndfil, cc, coulcut, dt
  real(dp)                          ::  dx, egap, ehomo, elumo
  real(dp)                          ::  kappa, scferror, trace_mult, vv(100)
  real(dp)                          ::  sumCubes, maxCH, Ef, smooth_maxCH, pnorm=6
  real(dp)                          ::  dvdw, d, gc(3), lasterror
  real(dp), allocatable             ::  FPUL(:,:), FSCOUL(:,:), FTOT(:,:), PairForces(:,:), trace(:)
  real(dp), allocatable             ::  SKForce(:,:), VX(:), VY(:), VZ(:)
  real(dp), allocatable             ::  charges_old(:), coul_forces(:,:), coul_forces_k(:,:), coul_forces_r(:,:)
  real(dp), allocatable             ::  coul_pot_k(:), coul_pot_r(:), dqin(:,:), dqout(:,:)
  real(dp), allocatable             ::  eigenvals(:), gbnd(:), n(:), n_0(:)
  real(dp), allocatable             ::  n_1(:), n_2(:), n_3(:), n_4(:)
  real(dp), allocatable             ::  n_5(:), onsitesH(:,:), onsitesS(:,:), rhoat(:)
  real(dp), allocatable             ::  origin(:), row(:)
  type(bml_matrix_t)                ::  aux_bml, dH0x_bml, dH0y_bml, dH0z_bml
  type(bml_matrix_t)                ::  dSx_bml, dSy_bml, dSz_bml, eigenvects
  type(bml_matrix_t)                ::  g_bml, ham0_bml, ham_bml, orthoh_bml
  type(bml_matrix_t)                ::  orthop_bml, over_bml, rho_bml, rhoat_bml
  type(bml_matrix_t)                ::  rhoh_bml, zmat_bml
  type(bml_matrix_t)                ::  copy_g_bml, gcov_bml
  type(graph_partitioning_t)        ::  gp
  type(graph_partitioning_t)        ::  gpat
  type(gsp2data_type)               ::  gsp2
  type(md_type)                     ::  md
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

  type(bml_matrix_t)                :: ZK1_bml, ZK2_bml, ZK3_bml
  type(bml_matrix_t)                :: ZK4_bml, ZK5_bml, ZK6_bml
  integer                           :: igenz !Couter to keep track of the times zmat is computed.
  type(genZSPinp)                   :: zsp

  integer, allocatable              :: xadj_cov(:), adjncy_cov(:), CH_count_cov(:)
  integer, allocatable              :: part_cov(:), core_count_cov(:), Halo_count_cov(:,:)
  integer                           :: vsize(2)
  integer                           :: nparts_cov


!!!!!!!!!!!!!!!!!!!!!!!!
  !> Main program driver
!!!!!!!!!!!!!!!!!!!!!!!!

  call gpmd_Init()

  call gpmd_FirstCharges()

  call gpmd_SCFLoop(lt%maxscf,sy%net_charge)

  call prg_write_system(sy,adjustl(trim(lt%jobname))//"_SCFcharges","pdb")

  call gpmd_EnergAndForces(sy%net_charge)

  call gpmd_PrepareMD()

  call gpmd_MDloop()

  call gpmd_Finalize()

contains

  !> Initialize gpmd.
  !!
  subroutine gpmd_Init()
    implicit none

    !> Start progress
    call prg_progress_init()

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
    call prg_parse_sp2(sp2,inputfile)

    !> Parsing GSP2 input paramenters. This will read the variables in the input file.
    !  gsp2 is the "gsp2data_type".
    call prg_parse_gsp2(gsp2,inputfile)

    !> Parsing MD input paramenters. This will read the variables in the input file.
    !  md is the "md_type".
    call parse_md(md,inputfile)

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

    !> Center sytem inside the box and fold it by with lattice_vectors.
    call prg_translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin)
    call prg_write_system(sy,adjustl(trim(lt%jobname))//"_centered","pdb")

    !> Building the neighbor list.
    call get_coulcut(lt%coul_acc,lt%timeratio,sy%nats,sy%lattice_vector,coulcut)
    call build_nlist(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)

    !> Allocate bounds vector.
    allocate(gbnd(2))

    !> Get Huckel hamiltonian. Computes the Extended Huckel Hamiltonian from the
    !  atom coordinates. The main inputs are the huckelTBparams and the system coordinate (sy%coordinate)
    !  The main outputs are Hamiltonian (ham_bml) and Overlap (over_bml) matrices.
    !   call get_hshuckel(ham_bml,over_bml,sy%coordinate,sy%spindex,sy%spatnum,&
    !     "../../huckelTBparams",lt%bml_type,lt%mdim,lt%threshold&
    !     ,tb%nsp,tb%splist,tb%basis,tb%numel,tb%onsite_energ,&
    !     tb%norbi,tb%hubbardu)

    !> LATTE Hamiltonian
    call load_latteTBparams(tb,sy%splist,lt%parampath)

    !> Get the mapping of the Hamiltonian index with the atom index
    !  hindex(1,i)=starting Hindex for atom i.
    !  hindex(2,i)=final Hindex for atom i.
    !  file://~/progress/build/doc/html/ham__latte__mod_8F90_source.html
    call get_hindex(sy%spindex,tb%norbi,hindex,norb)

    call load_bintTBparamsH(sy%splist,tb%onsite_energ,&
         typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,lt%parampath)
    call write_bintTBparamsH(typeA,typeB,&
         intKind,intPairsH,intPairsS,adjustl(trim(lt%jobname))//"_bondints.nonorth")
    ! call prg_init_hsmat(ham_bml,over_bml,lt%bml_type,lt%mdim,norb)
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,ham_bml)
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,over_bml)

    call get_hsmat(ham_bml,over_bml,sy%coordinate,sy%lattice_vector,sy%spindex,&
         tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)

    nnodes = size(hindex, dim=2)
    allocate(hnode(nnodes))
    do i = 1, nnodes
      hnode(i) = hindex(1, i)
    enddo

    !> Get the reciprocal vectors
    call prg_get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)

    if(lt%verbose.GE.5)then
      if (printRank() .eq. 1) then
        write(*,*)
        call bml_print_matrix("H0",ham_bml,0,6,0,6)
        call bml_print_matrix("S",over_bml,0,6,0,6)
      endif
    endif

    !> Get occupation based on last shell population.
    !  WARNING: This could change depending on the TB method being used.
    nel = sum(element_numel(sy%atomic_number(:)),&
         size(sy%atomic_number,dim=1))
    bndfil = nel/(2.0_dp*norb)

    if (printRank() .eq. 1) then
      write(*,*) "Filling fraction = ", bndfil
      write(*,*) "Number of electrons = ", nel
    endif

    !> Load Pair potentials for latte TB.
    call load_PairPotTBparams(lt%parampath,sy%splist,ppot)

    call prg_build_atomic_density(rhoat_bml,tb%numel,hindex,sy%spindex,norb,&
         lt%bml_type)

    if(lt%verbose.GE.2) call bml_print_matrix("rhoat_bml",rhoat_bml,0,6,0,6)

    call prg_init_ZSPmat(igenz,zk1_bml,zk2_bml,zk3_bml&
         ,zk4_bml,zk5_bml,zk6_bml,norb,lt%bml_type)

  end subroutine gpmd_Init


  !> First Charge computation.
  !! \brief Here we compute the non-scf charges.
  !!
  subroutine gpmd_FirstCharges()
    implicit none

    !> Initialize the density matrix (rho_bml) and inverse overlap factor (zmat_bml).
    call prg_init_pzmat(rho_bml,zmat_bml,lt%bml_type,lt%mdim,norb)

    !> Get the Inverse square root overlap matrix.
    call prg_timer_start(dyn_timer,"Build Z")
    call gpmd_buildz()
    call prg_timer_stop(dyn_timer,1)

    !> Initialize the orthogonal versions of ham and rho.
    call prg_init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)

    !> Orthogonalize ham.
    call prg_timer_start(ortho_timer)
    call prg_orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
         lt%threshold,lt%bml_type,lt%verbose)
    call prg_timer_stop(ortho_timer)

#ifdef DO_MPI_BLOCK
    call prg_allGatherParallel(orthoh_bml)
#endif

    !> Calculate gershgorin bounds
    call bml_gershgorin(orthoh_bml, gbnd)
    if (printRank() .eq. 1) then
      write(*,*)
      write(*,*) "Gershgorin: mineval = ", gbnd(1), " maxeval = ", gbnd(2)
    endif

    !> SP2 algorithm to get the density matrix
    pp=0
    vv=0.0_dp

    if(lt%method.NE."Diag")then
      call prg_timer_start(sp2_timer)
      call prg_sp2_alg2_genseq(orthoh_bml,orthop_bml,lt%threshold,bndfil,sp2%minsp2iter,&
           sp2%maxsp2iter,sp2%sp2conv,sp2%sp2tol, pp, icount, vv,lt%verbose)
      call prg_timer_stop(sp2_timer)

#ifdef DO_MPI_BLOCK
      call prg_allGatherParallel(orthop_bml)
#endif

      if (printRank() .eq. 1) then
        write(*,*)
        write(*,*) "SP2 vv:"
        do i = 1, icount
          write(*,*) "i = ", i, " vv = ", vv(i)
        enddo
        write(*,*)
      endif

      !> Calculate Homo-Lumo gap
      call prg_homolumogap(vv, icount, pp, gbnd(1), gbnd(2), ehomo, elumo, egap,lt%verbose)
      if (printRank() .eq. 1) then
        write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
             " egap = ", egap
        write(*,*)
      endif

    else
      call gpmd_RhoSolver
    endif

    !> Save a copy of the density matrix adjacency matrix for the first graph
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,g_bml)
    call bml_copy(orthop_bml, g_bml)

    !> Deprg_orthogonalize rho.
    call prg_timer_start(deortho_timer)
    call prg_deorthogonalize(orthop_bml,zmat_bml,rho_bml,&
         lt%threshold,lt%bml_type,lt%verbose)
    call prg_timer_stop(deortho_timer)
#ifdef DO_MPI_BLOCK
    call prg_allGatherParallel(rho_bml)
#endif

    if(lt%verbose.GE.2)then
      if (printRank() .eq. 1) then
        call bml_print_matrix("rho0",rho_bml,0,6,0,6)
      endif
    endif

    !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing
    !  the charges.
    call prg_get_charges(rho_bml, over_bml, hindex, sy%net_charge, tb%numel, sy%spindex, lt%mdim, lt%threshold)
    charges_old = sy%net_charge

    if (printRank() .eq. 1) then
      write(*,*)""
      write(*,*)"Total charge =", sum(sy%net_charge(:),size(sy%net_charge,dim=1))
    endif

    if(lt%verbose.GE.2)then
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"System charges:"
        do j=1,sy%nats
          write(*,*)sy%symbol(j),sy%net_charge(j)
        enddo
      endif
    endif

    !call prg_get_eigenvalues(ham_bml,eigenvals,10)
    !call prg_write_tdos(eigenvals, 0.1_dp, 1000, -40.0_dp, 20.0_dp, "tdos.dos")
    !stop

  end subroutine gpmd_FirstCharges


  !>  SCF loop
  !!
  subroutine gpmd_SCFLoop(Nr_SCF,nguess)
    implicit none
    integer :: Nr_SCF
    real(dp):: nguess(:)

    sy%net_charge = nguess

    !> Save actual hamiltonian as the non-scf Hamiltonian (H0)
    call bml_copy_new(ham_bml,ham0_bml)

    charges_old = sy%net_charge

    !> Beginning of the SCF loop.
    do i=1,Nr_SCF

      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"SCF iter", i
      endif

      !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"In real Coul ..."
      endif

      ! call get_ewald_real(sy%spindex,sy%splist,sy%coordinate&
      !   ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
      !   sy%volr,lt%coul_acc,coul_forces_r,coul_pot_r);
      call prg_timer_start(dyn_timer,"Real coul")
      call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate&
           ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
           nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call prg_timer_stop(dyn_timer,1)

      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"In recip Coul ..."
      endif

      call prg_timer_start(dyn_timer,"Recip coul")
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
      call prg_timer_stop(dyn_timer,1)

      !> Get the scf hamiltonian. The outputs is ham_bml.
      if (printRank() .eq. 1) then
        write(*,*)""; write(*,*)"In prg_get_hscf ..."
      endif
      call prg_get_hscf(ham0_bml,over_bml,ham_bml,sy%spindex,hindex,tb%hubbardu,sy%net_charge,&
           coul_pot_r,coul_pot_k,lt%mdim,lt%threshold)

      !> Initialize the orthogonal versions of ham and rho.
      call prg_init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)

      !> Orthogonalize the Hamiltonian
      !    if (printRank() .eq. 1) then
      write(*,*)""; write(*,*)"In prg_orthogonalize H ..."
      !    endif
      call prg_timer_start(ortho_timer)
      call prg_orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
           lt%threshold,lt%bml_type,lt%verbose)
      call prg_timer_stop(ortho_timer,1)
#ifdef DO_MPI_BLOCK
      call prg_allGatherParallel(orthoh_bml)
#endif

      if(lt%verbose.GE.2)then
        if (printRank() .eq. 1) then
          call bml_print_matrix("orthoh_bml",orthoh_bml,0,6,0,6)
          call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
          call bml_print_matrix("ham0_bml",ham0_bml,0,6,0,6)
          call bml_print_matrix("zmat_bml",zmat_bml,0,6,0,6)
          call bml_print_matrix("orthop_bml",orthop_bml,0,6,0,6)
        endif
      endif

      if(lt%method.EQ."GSP2")call gpmd_graphpart()


      !> Calculate gershgorin bounds
      call bml_gershgorin(orthoh_bml, gbnd)
      gp%mineval = gbnd(1)
      gp%maxeval = gbnd(2)
      if (printRank() .eq. 1) then
        write(*,*)""
        write(*,*) "Gershgorin: mineval = ", gbnd(1), " maxeval = ", gbnd(2)
      endif

      write(*,*)ehomo,elumo,gp%mineval, gp%maxeval,gsp2%errlimit

      ! stop

      if(lt%method.eq."GSP2")then
        !! Calculate SP2 sequence
        call prg_sp2sequence(gp%pp, gp%maxIter, gp%mineval, gp%maxeval, ehomo, elumo, &
             gsp2%errlimit)
        if (printRank() .eq. 1) then
          write(*,*)
          write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter
          write(*,*)
        endif
      endif

      !> Now sove for the desity matrix.
      call gpmd_rhosolver()


      if(lt%method.eq."GSP2")then
        !> Calculate Homo-Lumo gap
        call prg_homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, egap)

        if (printRank() .eq. 1) then
          write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
               " egap = ", egap
          write(*,*)
        endif
      endif


      !! Calculate Trace[HP]
      trace_mult = bml_trace_mult(orthoh_bml, orthop_bml)
#ifdef DO_MPI_BLOCK
      if (getNRanks() > 1) then
        call prg_sumRealReduce(trace_mult)
      endif
#endif
      if (printRank() .eq. 1) then
        write(*,*) "Trace[HP] for SP2 = ", trace_mult
        write(*,*) "Band energy per atom = ", trace_mult/gsp2%natoms
        write(*,*)
      endif

      !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
      call prg_timer_start(deortho_timer)
      call prg_deorthogonalize(orthop_bml,zmat_bml,rho_bml,&
           lt%threshold,lt%bml_type,lt%verbose)
      call prg_timer_stop(deortho_timer)
#ifdef DO_MPI_BLOCK
      call prg_allGatherParallel(rho_bml)
#endif

      !> Copy density matrix for the graph for the next iteration
      call bml_copy(orthop_bml, g_bml)

      call prg_destroyGraphPartitioning(gp)

      call bml_deallocate(orthop_bml)

      !> Get the system charges.
      call prg_get_charges(rho_bml,over_bml,hindex,sy%net_charge,tb%numel,&
           sy%spindex,lt%mdim,lt%threshold)

      !if (printRank() .eq. 1) then
      write(*,*)"Total charge", sum(sy%net_charge(:)),size(sy%net_charge,dim=1)
      !endif

      call prg_qmixer(sy%net_charge,charges_old,dqin,&
           dqout,scferror,i,lt%pulaycoeff,lt%mpulay,lt%verbose)

      ! call prg_linearmixer(sy%net_charge,charges_old,scferror,lt%mixcoeff,lt%verbose)

      charges_old = sy%net_charge

      if(lt%verbose.GE.1)then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"System charges:"
          do j=1,sy%nats
            write(*,*)sy%symbol(j),sy%net_charge(j)
          enddo
        endif
      endif

      if(scferror.lt.lt%scftol.and.i.gt.5) then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"SCF converged within",i,"steps ..."
          write(*,*)"SCF error =",scferror
        endif
        exit
      endif

    enddo
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

  !> Get the energies and forces.
  !!
  subroutine gpmd_EnergAndForces(charges)
    Implicit none
    real(dp) :: charges(:)

    !> Get Repulsive energy and forces
    call prg_timer_start(dyn_timer,"Get pair pot")
    call get_PairPot_contrib(sy%coordinate,sy%lattice_vector,sy%spindex,ppot,PairForces,ERep)
    call prg_timer_stop(dyn_timer,1)
    write(*,*)"Energy Repulsive = ", ERep

    !> Get Coulombic energy
    ECoul = 0.0;
    do i = 1,sy%nats
      ECoul = ECoul + charges(i)*(tb%hubbardu(sy%spindex(i))*charges(i) + coul_pot_r(i) + coul_pot_k(i) );
    enddo
    write(*,*)"Energy Coulomb = ", ECoul

    !> Get Electronic energy
    call bml_copy_new(rho_bml,aux_bml)
    call bml_add_deprecated(1.0_dp,aux_bml,-1.0_dp,rhoat_bml,lt%threshold)
    TRRHOH  = bml_trace_mult(aux_bml, ham_bml)
    write(*,*)"Energy Band = ", TRRHOH
    call bml_deallocate(aux_bml)

    Etot = TRRHOH - 0.5_dp*ECoul  + ERep
    write(*,*)"Energy Electronic =", Etot

    EPOT = Etot;

    coul_forces =  coul_forces_r + coul_forces_k

    dx = 0.0001_dp;

    call prg_timer_start(dyn_timer,"get_dH and get_dS")
    call get_dH(dx,sy%coordinate,hindex,sy%spindex,intPairsH,onsitesH,sy%symbol,&
         sy%lattice_vector, norb, tb%norbi, lt%bml_type, &
         lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)

    if(lt%verbose.GE.10)then
      call bml_print_matrix("dH0x_bml",dH0x_bml,0,10,0,10)
      call bml_print_matrix("dH0y_bml",dH0y_bml,0,10,0,10)
      call bml_print_matrix("dH0z_bml",dH0z_bml,0,10,0,10)
    endif

    ! stop
    call prg_timer_start(dyn_timer,"SK forces")
    call get_skforce(sy%nats,rho_bml,dH0x_bml,dH0y_bml,&
         dH0z_bml,hindex,SKForce,lt%threshold)
    call prg_timer_stop(dyn_timer,1)

    ! Deallocate these matrices to avoid memory problems.
    call bml_deallocate(dH0x_bml)
    call bml_deallocate(dH0y_bml)
    call bml_deallocate(dH0z_bml)

    call get_dS(dx,sy%coordinate,hindex,sy%spindex,intPairsS,onsitesS,sy%symbol,&
         sy%lattice_vector, norb, tb%norbi, lt%bml_type, &
         lt%threshold, dSx_bml,dSy_bml,dSz_bml)
    call prg_timer_stop(dyn_timer,1)

    if(lt%verbose.GE.10)then
      call bml_print_matrix("dSx_bml",dSx_bml,0,10,0,10)
      call bml_print_matrix("dSy_bml",dSy_bml,0,10,0,10)
      call bml_print_matrix("dSz_bml",dSz_bml,0,10,0,10)
    endif

    call prg_timer_start(dyn_timer,"Pulay forces")
    call prg_get_pulayforce(sy%nats,zmat_bml,ham_bml,rho_bml,&
         dSx_bml,dSy_bml,dSz_bml,hindex,FPUL,lt%threshold)
    call prg_timer_stop(dyn_timer,1)

    call prg_timer_start(dyn_timer,"Non orth Coul F")
    call get_nonortho_coul_forces(sy%nats, norb, dSx_bml,dSy_bml,dSz_bml,&
         hindex,sy%spindex,rho_bml,charges,coul_pot_r,coul_pot_k,tb%hubbardu,FSCOUL,lt%threshold)
    call prg_timer_stop(dyn_timer,1)

    ! Deallocate these matrices to avoid memory problems.
    call bml_deallocate(dSx_bml)
    call bml_deallocate(dSy_bml)
    call bml_deallocate(dSz_bml)

    if(.not.allocated(FTOT))allocate(FTOT(3,sy%nats))

    FTOT = SKForce + PairForces + FPUL + coul_forces + FSCOUL

    if(lt%verbose.GE.10)then
      write(*,*)""; write(*,*)"SK forces:"
      do i = 1,sy%nats
        write(*,*)i,SKForce(1,i),SKForce(2,i),SKForce(3,i)
      enddo

      write(*,*)""; write(*,*)"Pulay forces:"
      do i = 1,sy%nats
        write(*,*)i,FPUL(1,i),FPUL(2,i),FPUL(3,i)
      enddo

      write(*,*)""; write(*,*)"Coulomb forces:"
      do i = 1,sy%nats
        write(*,*)"Coul Force",i,coul_forces(1,i),coul_forces(2,i),coul_forces(3,i)
      enddo

      write(*,*)""; write(*,*)"Repulsive forces:"
      do i = 1,sy%nats
        write(*,*)i,PairForces(1,i),PairForces(2,i),PairForces(3,i)
      enddo

      write(*,*)""; write(*,*)"Nonorthogonal Coulomb forces:"
      do i = 1,sy%nats
        write(*,*)i,FSCOUL(1,i),FSCOUL(2,i),FSCOUL(3,i)
      enddo

      write(*,*)""; write(*,*) "Total forces"
      do i = 1,sy%nats
        write(*,*)i,FTOT(1,i),FTOT(2,i),FTOT(3,i)
      enddo
    endif

  end subroutine gpmd_EnergAndForces


  !>  Main MD loop
  !!  This routine performs the MD loops up to "ls%mdsteps"
  !!
  subroutine gpmd_MDloop()
    implicit none

    do MD_step = 1,lt%mdsteps

      write(*,*)""
      write(*,*)"         #######################"
      write(*,*)"           MDStep =",MD_step
      write(*,*)"         #######################"
      write(*,*)""

      !> Get Kinetic energy
      EKIN = 0.0_dp
      do i=1,sy%nats
        EKIN = EKIN + sy%mass(i)*(VX(i)**2+VY(i)**2+VZ(i)**2)
      enddo
      EKIN = 0.5_dp*MVV2KE*EKIN

      !! Statistical temperature in Kelvin
      Temp = (2.0_dp/3.0_dp)*KE2T*EKIN/real(sy%nats,dp);
      !! Total Energy in eV
      Energy = EKIN + EPOT;
      !! Time in fs
      Time = MD_step*lt%timestep;

      write(*,*)"Time [fs] = ",Time
      write(*,*)"Energy Kinetic [eV] = ",EKIN
      write(*,*)"Energy Potential [eV] = ",EPOT
      write(*,*)"Energy Total [eV] = ",Energy
      write(*,*)"Temperature [K] = ",Temp

      !> First 1/2 of Leapfrog step
      call prg_timer_start(dyn_timer,"Half Verlet")
      call halfVerlet(sy%mass,FTOT,lt%timestep,VX,VY,VZ)
      call prg_timer_stop(dyn_timer,1)


      if(lt%verbose.GE.5)then
        do i = 1,sy%nats
          write(*,*)"Velocities",i,VX(i),VY(i),VZ(i)
        enddo
      endif

      !> Update positions
      call prg_timer_start(dyn_timer,"Update positions")
      call updatecoords(origin,sy%lattice_vector,lt%timestep,VX,VY,VZ,sy%coordinate)
      call prg_timer_stop(dyn_timer,1)

      !> Update neighbor list
      call prg_timer_start(dyn_timer,"Build Nlist")
      call build_nlist(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)
      call prg_timer_stop(dyn_timer,1)

      !> Get new H0 and S
      call prg_timer_start(dyn_timer,"Build H0 and S")
      call get_hsmat(ham_bml,over_bml,sy%coordinate,sy%lattice_vector,sy%spindex,&
           tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)
      call prg_timer_stop(dyn_timer,1)


      !> Get the Inverse square root overlap matrix.
      call prg_timer_start(dyn_timer,"Build Z")
      call gpmd_buildz()
      call prg_timer_stop(dyn_timer,1)

      if(lt%verbose.GE.2)then
        call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
        call bml_print_matrix("over_bml",over_bml,0,6,0,6)
        call bml_print_matrix("zmat_bml",zmat_bml,0,6,0,6)
      endif

      if(lt%verbose.GE.1)then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"before xlbo:"
          do j=1,sy%nats
            write(*,*)sy%symbol(j),sy%net_charge(j)
          enddo
        endif
      endif

      call prg_xlbo_nint(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,MD_step,xl)

      if(lt%verbose.GE.1)then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"after xlbo:"
          do j=1,sy%nats
            write(*,*)sy%symbol(j),n(j)
          enddo
        endif
      endif

      Nr_SCF_It = xl%maxscfiter;

      !> Use SCF the first M_prg_init MD steps
      if(MD_step < xl%minit) Nr_SCF_It = xl%maxscfInitIter

      !> SCF loop
      call gpmd_SCFLoop(Nr_SCF_It,n)

      lasterror = scferror

      !> Recompute the neighbors list
      call prg_timer_start(dyn_timer,"Real coul")
      call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate&
           ,n,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
           nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call prg_timer_stop(dyn_timer,1)

      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      if (printRank() .eq. 1) then
        write(*,*)"In recip Coul ..."
      endif

      call prg_timer_start(dyn_timer,"Recip coul")
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,n,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
      call prg_timer_stop(dyn_timer,1)

      coul_forces = coul_forces_r + coul_forces_k

      call prg_timer_start(dyn_timer,"Get Hscf")
      call prg_get_hscf(ham0_bml,over_bml,ham_bml,sy%spindex,hindex,tb%hubbardu,n,&
           coul_pot_r,coul_pot_k,lt%mdim,lt%threshold)
      call prg_timer_stop(dyn_timer,1)

      !> Initialize the orthogonal versions of ham and rho.
      call prg_init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)

      !> Orthogonalize the Hamiltonian
      if (printRank() .eq. 1) then
        write(*,*)"In prg_orthogonalize H ..."
      endif


      call prg_timer_start(ortho_timer)
      call prg_orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
           lt%threshold,lt%bml_type,lt%verbose)
      call prg_timer_stop(ortho_timer,1)

#ifdef DO_MPI_BLOCK
      call prg_allGatherParallel(orthoh_bml)
#endif

      if(lt%verbose.GT.2)then
        if (printRank() .eq. 1) then
          call bml_print_matrix("orthoh_bml",orthoh_bml,0,6,0,6)
          call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
          call bml_print_matrix("zmat_bml",zmat_bml,0,6,0,6)
        endif
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !       if(MD_step == 10)then
      !         allocate(trace(norb))
      !         call bml_threshold(orthoh_bml,5.0d-3)
      !         if(bml_get_N(aux_bml) > 1)call bml_deallocate(aux_bml)
      !         call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,aux_bml)
      !
      !         call bml_multiply(orthoh_bml, orthoh_bml, aux_bml, 1.0_dp, 0.0_dp, 0.0_dp)
      !
      !         allocate(row(norb))
      !           count = 0
      !         do i=1,norb
      !           row = 0.0_dp
      !           call bml_get_row(aux_bml,i,row)
      !           do j=1,norb
      !             if(abs(row(j))>0) count = count + 1
      !           enddo
      !         enddo
      !
      !         write(*,*)count/norb
      !         stop
      !       endif

      if(lt%method.EQ."GSP2")call gpmd_graphpart()


      !> Calculate gershgorin bounds
      call bml_gershgorin(orthoh_bml, gbnd)
      gp%mineval = gbnd(1)
      gp%maxeval = gbnd(2)
      if (printRank() .eq. 1) then
        write(*,*)
        write(*,*) "Gershgorin: mineval = ", gbnd(1), " maxeval = ", gbnd(2)
      endif

      if(lt%method.eq."GSP2")then
        !! Calculate SP2 sequence
        call prg_sp2sequence(gp%pp, gp%maxIter, gp%mineval, gp%maxeval, ehomo, elumo, &
             gsp2%errlimit)
        if (printRank() .eq. 1) then
          write(*,*)
          write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter
          write(*,*)
        endif
      endif

      !> Now we construct the density matrix.
      call gpmd_rhosolver()

      if(lt%method.eq."GSP2")then
        !> Calculate Homo-Lumo gap
        call prg_homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, egap)

        if (printRank() .eq. 1) then
          write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
               " egap = ", egap
          write(*,*)
        endif
      endif

      !! Calculate Trace[HP]
      trace_mult = bml_trace_mult(orthoh_bml, orthop_bml)
#ifdef DO_MPI_BLOCK
      if (getNRanks() > 1) then
        call prg_sumRealReduce(trace_mult)
      endif
#endif
      if (printRank() .eq. 1) then
        write(*,*) "Trace[HP] for SP2 = ", trace_mult
        write(*,*) "Band energy per atom = ", trace_mult/gsp2%natoms
        write(*,*)
      endif

      !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
      call prg_timer_start(deortho_timer)
      call prg_deorthogonalize(orthop_bml,zmat_bml,rho_bml,&
           lt%threshold,lt%bml_type,lt%verbose)
      call prg_timer_stop(deortho_timer,1)
#ifdef DO_MPI_BLOCK
      call prg_allGatherParallel(rho_bml)
#endif

      !> Copy density matrix for the graph for the next iteration
      call bml_copy(orthop_bml, g_bml)

      call prg_destroyGraphPartitioning(gp)

      call bml_deallocate(orthop_bml)

      !> Get the system charges.
      call prg_get_charges(rho_bml,over_bml,hindex,sy%net_charge,tb%numel,&
           sy%spindex,lt%mdim,lt%threshold)

      !if (printRank() .eq. 1) then
      write(*,*)"Total charge", sum(sy%net_charge(:)),size(sy%net_charge,dim=1)
      !endif

      if(lt%verbose.GE.1)then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"System charges forint:"
          do j=1,sy%nats
            write(*,*)sy%symbol(j),sy%net_charge(j)
          enddo
        endif
      endif

      if(lt%verbose.GE.1)then
        if (printRank() .eq. 1) then
          write(*,*)""; write(*,*)"System ncharges:"
          do j=1,sy%nats
            write(*,*)sy%symbol(j),n(j)
          enddo
        endif
      endif

      call gpmd_EnergAndForces(n)

      !> Adjust forces for the linearized XLBOMD functional
      call prg_xlbo_fcoulupdate(Coul_Forces,sy%net_charge,n)

      !> Total XLBOMD force
      FTOT = SKForce + PairForces + FPUL + Coul_Forces + FSCOUL;

      !> Integrate second 1/2 of leapfrog step
      call halfVerlet(sy%mass,FTOT,lt%timestep,VX,VY,VZ)

      ! call prg_write_trajectory(sy,MD_step,md%writeeach,lt%timestep,adjustl(trim(lt%jobname))//"_traj","pdb")
      call prg_write_trajectory(sy,MD_step,md%writeeach,lt%timestep,adjustl(trim(lt%jobname))//"_traj","pdb")

    enddo
    ! End of MD loop.

  end subroutine gpmd_MDloop


  !> Finalize the gpmd
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
    call prg_progress_shutdown()

  end subroutine gpmd_Finalize


  !>  Preparing for MD
  !!
  subroutine gpmd_prepareMD()
    Implicit none

    !> Initialize velocities
    allocate(VX(sy%nats))
    allocate(VY(sy%nats))
    allocate(VZ(sy%nats))

    VX = 0.0_dp; VY = 0.0_dp; VZ = 0.0_dp;    ! Initialize velocities

    !> Kinetic energy in eV (MVV2KE: unit conversion)
    MVV2KE = 166.0538782_dp/1.602176487_dp

    KE2T = 1.0_dp/0.000086173435_dp

  end subroutine gpmd_prepareMD


  !> Solver for computing the density matrix
  !!
  subroutine gpmd_RhoSolver()
    implicit none

    call prg_timer_start(dyn_timer,"Solver")

    if(lt%method.EQ."GSP2")then
      call prg_timer_start(graphsp2_timer)
      call prg_subgraphSP2Loop(orthoh_bml, g_bml, orthop_bml, gp, lt%threshold)
      call prg_timer_stop(graphsp2_timer)
      ! call prg_sp2_alg1_seq(orthoh_bml,orthop_bml,lt%threshold, gp%pp, gp%maxIter, gp%vv)
    elseif(lt%method.EQ."SP2")then
      if(sp2%flavor.eq."Basic")then
        call prg_sp2_basic(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
             ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
      elseif(sp2%flavor.eq."Alg1")then
        call prg_sp2_alg1(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
             ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
      elseif(sp2%flavor.eq."Alg2")then
        call prg_sp2_alg2(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
             ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
      else
        stop "No valid SP2 flavor"
      endif
    elseif(lt%method.EQ."Diag")then
      call prg_build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
      ! call prg_build_density_T(orthoh_bml,orthop_bml,lt%threshold,bndfil, 0.1_dp, Ef)
      write(*,*)"Fermi Level =",Ef
    else
      stop "No valid Method in LATTE parameters"
    endif

    call prg_timer_stop(dyn_timer,1)

#ifdef DO_MPI_BLOCK
    call prg_allGatherParallel(orthop_bml)
#endif

    if(lt%verbose.GE.2)then
      if (printRank() .eq. 1) then
        call bml_print_matrix("gsp2 orthop_bml",orthop_bml,0,6,0,6)
      endif
    endif

  end subroutine gpmd_RhoSolver


  !> Driver for the grpah-partitioning.
  !!
  subroutine gpmd_graphpart
    implicit none

    !> Symmetrize and Threshold the Matrix
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,copy_g_bml)
    call bml_threshold(g_bml, gsp2%gthreshold)
    call bml_transpose_new(g_bml, copy_g_bml)
    call bml_add_deprecated(0.5_dp,g_bml,0.5_dp,copy_g_bml,0.0_dp)
    call bml_threshold(g_bml, gsp2%gthreshold)
    call bml_deallocate(copy_g_bml)
#ifdef DO_MPI_BLOCK
    call prg_allGatherParallel(g_bml)
#endif

    !> Create graph partitioning - Use Block or METIS or METIS+SA or METIS+KL
    call prg_timer_start(part_timer)

    !> Block partitioning
    if (gsp2%partition_type == "Block") then
      !> Partition by orbital or atom
      if (gsp2%graph_element == "Orbital") then
        call prg_equalPartition(gp, gsp2%nodesPerPart, norb)
      else
        call prg_equalGroupPartition(gp, hindex, nnodes, gsp2%nodesPerPart, norb)
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
          call prg_metisPartition(gp, nnodes, norb, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
        case("METIS+SA")
          call prg_metisPartition(gp, nnodes, norb, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
          call prg_simAnnealing(gp, xadj, adjncy, part, core_count, CH_count, &
               Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
        case("METIS+KL")
          call prg_metisPartition(gp, nnodes, norb, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
          call prg_KernLin2(gp, xadj, adjncy, part, core_count, CH_count, &
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
            call prg_simAnnealing(gp, xadj, adjncy, part, core_count, CH_count, &
                 Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
          case("KL")
            call prg_KernLin2(gp, xadj, adjncy, part, core_count, CH_count, &
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

    call prg_timer_stop(part_timer)

  end subroutine gpmd_graphpart


  !> Driver for the construction of Z
  !!
  subroutine gpmd_buildz()
    implicit none

    igenz = igenz + 1

    write(*,*)zsp%nfirst

    if(lt%zmat.eq."ZSP")then !Congruence transformation.

      call prg_buildzsparse(over_bml,zmat_bml,igenz,lt%mdim,&
           lt%bml_type, zk1_bml,zk2_bml,zk3_bml&
           ,zk4_bml,zk5_bml,zk6_bml,zsp%nfirst,zsp%nrefi,zsp%nreff,&
           zsp%numthresi,zsp%numthresf,zsp%integration,zsp%verbose)

    else

      !Build Z matrix using diagonalization (usual method).
      call prg_buildzdiag(over_bml,zmat_bml,lt%threshold,lt%mdim,lt%bml_type)

      ! stop
    endif

  end subroutine gpmd_buildz


end program mdresponse
