!> High-level program to perform SFC cycles in Extended Huckel Hamiltonian.
!!
program gploop

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
  ! Graph partitioning modules
  use prg_parallel_mod
  use prg_timer_mod
  use prg_graphsp2parser_mod
  use prg_sp2_mod
  use prg_graph_mod
  use prg_partition_mod
  use prg_subgraphLoop_mod
  use prg_homolumo_mod

  implicit none

  integer, parameter     ::  dp = kind(1.0d0)
  integer                ::  seed = 1
  integer                ::  i, j, k, nel, norb, nnodes
  integer                ::  nparts, niter=500
  integer, allocatable   ::  hindex(:,:), hnode(:)
  integer, allocatable   ::  xadj(:), adjncy(:), CH_count(:)
  integer, allocatable   ::  part(:), core_count(:), Halo_count(:,:)

  real(dp)               ::  bndfil, scferror
  real(dp)               ::  sumCubes, maxCH, smooth_maxCH, pnorm=6
  real(dp), allocatable  ::  charges_old(:), coul_forces_k(:,:), coul_forces_r(:,:), coul_pot_k(:)
  real(dp), allocatable  ::  coul_pot_r(:), dqin(:,:), dqout(:,:), eigenvals(:)

  type(bml_matrix_t)     ::  ham0_bml, ham_bml, orthoh_bml, orthop_bml
  type(bml_matrix_t)     ::  over_bml, rho_bml, zmat_bml, g_bml
  type(bml_matrix_t)     ::  copy_g_bml
  type(latte_type)       ::  lt
  type(gsp2data_type)    ::  gsp2
  type(system_type)      ::  sy
  type(tbparams_type)    ::  tb
  type(graph_partitioning_t) :: gp

  integer :: icount
  integer :: pp(100)
  real(dp) :: vv(100)
  real(dp) :: ehomo, elumo, egap, trace_mult
  real(dp), allocatable :: gbnd(:)

  logical :: first_part = .true.

  !> Start progress
  call prg_progress_init()

  if (printRank() .eq. 1) then
    write(*,*) "gploop start ..."
  endif

  !> Parsing input file. This file contains all the variables needed to
  !  run the scf including the sp2 (solver) variables. lt is "latte_type" structure
  !  containing all the variables.
  !  file://~/progress/build/doc/html/structlatteparser__latte__mod_1_1latte__type.html
  call parse_latte(lt,"input.in")

  !> Parsing system coordinates. This reads the coords.pdb file to get the position of every
  !  atom in the system. sy is the "system_type" structure containing all the variables.
  !  file://~/progress/build/doc/html/classsystem__latte__mod.html
  call prg_parse_system(sy,"coords","pdb")

  !> Allocate bounds vactor.
  allocate(gbnd(2))

  !> Get Huckel hamiltonian. Computes the Extended Huckel Hamiltonian from the
  !  atom coordinates. The main inputs are the huckelTBparams and the system coordinate (sy%coordinate)
  !  The main outputs are Hamiltonian (ham_bml) and Overlap (over_bml) matrices.
  call get_hshuckel(ham_bml,over_bml,sy%coordinate,sy%spindex,sy%spatnum,&
       "../../huckelTBparams",lt%bml_type,lt%mdim,lt%threshold&
       ,tb%nsp,tb%splist,tb%basis,tb%numel,tb%onsite_energ,&
       tb%norbi,tb%hubbardu)

  !> Get the mapping of the Hamiltonian index with the atom index
  !  hindex(1,i)=starting Hindex for atom i.
  !  hindex(2,i)=final Hindex for atom i.
  !  file://~/progress/build/doc/html/ham__latte__mod_8F90_source.html
  call get_hindex(sy%spindex,tb%norbi,hindex,norb)

  nnodes = size(hindex, dim=2)
  allocate(hnode(nnodes))
  do i = 1, nnodes
    hnode(i) = hindex(1, i)
  enddo

  if (printRank() .eq. 1) then
    write(*,*) "Number of orbitals = ", norb
    write(*,*)
    call bml_print_matrix("ham0_bml",ham_bml,0,6,0,6)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>  First Charge computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the density matrix (rho_bml) and inverse overlap factor (zmat_bml).
  call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml)
  call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,zmat_bml)

  !> Get the Inverse square root overlap matrix.
  call prg_buildzdiag(over_bml,zmat_bml,lt%threshold,lt%mdim,lt%bml_type)

  !> Initialize the orthogonal versions of ham and rho.
  call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthoh_bml)
  call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,orthop_bml)

  !> Orthogonalize ham.
  call prg_timer_start(ortho_timer)
  call prg_orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
       lt%threshold,lt%bml_type,lt%verbose)
  call prg_timer_stop(ortho_timer)

#ifdef DO_MPI_BLOCK
  call prg_allGatherParallel(orthoh_bml)
#endif

  !> The SP2 algorithm is used to get the first orthogonal Density matrix (orthop).
  call prg_parse_gsp2(gsp2,"input.in")

  !> Calculate gershgorin bounds
  call bml_gershgorin(orthoh_bml, gbnd)
  if (printRank() .eq. 1) then
    write(*,*) "Gershgorin: mineval = ", gbnd(1), " maxeval = ", gbnd(2)
    write(*,*)
  endif

  !> SP2 algorithm.
  call prg_timer_start(sp2_timer)
  call prg_sp2_alg2_genseq(orthoh_bml,orthop_bml,lt%threshold,bndfil,gsp2%minsp2iter,&
       gsp2%maxsp2iter,gsp2%sp2conv,gsp2%sp2tol, pp, icount, vv)
  call prg_timer_stop(sp2_timer)
#ifdef DO_MPI_BLOCK
  call prg_allGatherParallel(orthop_bml)
#endif

  if (printRank() .eq. 1) then
    write(*,*)
    write(*,*) "SP2 count = ", icount
    do i = 1, icount
      write(*,*) "i = ", i, " pp = ", pp(i), " vv = ", vv(i)
    enddo
    write(*,*)
  endif

  !> Calculate Homo-Lumo gap
  call prg_homolumogap(vv, icount, pp, gbnd(1), gbnd(2), ehomo, elumo, egap)
  if (printRank() .eq. 1) then
    write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
         " egap = ", egap
    write(*,*)
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

  if (printRank() .eq. 1) then
    call bml_print_matrix("rho_bml",rho_bml,0,6,0,6)
  endif

  !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing
  !  the charges.
  call prg_get_charges(rho_bml, over_bml, hindex, sy%net_charge, tb%numel, sy%spindex, lt%mdim, lt%threshold)
  charges_old = sy%net_charge

  if (printRank() .eq. 1) then
    write(*,*)"Total charges =", sum(sy%net_charge(:),size(sy%net_charge,dim=1))
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>  SCF loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Save actual hamiltonian as the non-scf Hamiltonian (H0)
  call bml_copy_new(ham_bml,ham0_bml)

  !> Get the reciprocal vector (this is needed to compute the Coulombic interactions)
  call prg_get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)

  !> Beginning of the SCF loop.
  do i=1,lt%maxscf

    if (printRank() .eq. 1) then
      write(*,*)"SCF iter", i
    endif

    !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
    if (printRank() .eq. 1) then
      write(*,*)"In real Coul ..."
    endif
    call get_ewald_real(sy%spindex,sy%splist,sy%coordinate&
         ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
         sy%volr,lt%coul_acc,coul_forces_r,coul_pot_r);

    !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
    if (printRank() .eq. 1) then
      write(*,*)"In recip Coul ..."
    endif
    call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
         ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
         sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);

    !> Get the scf hamiltonian. The outputs is ham_bml.
    if (printRank() .eq. 1) then
      write(*,*)"in prg_get_hscf ..."
    endif
    call prg_get_hscf(ham0_bml,over_bml,ham_bml,sy%spindex,hindex,tb%hubbardu,sy%net_charge,&
         coul_pot_r,coul_pot_k,lt%mdim,lt%threshold)

    !> Orthogonalize the Hamiltonian
    if (printRank() .eq. 1) then
      write(*,*)"in prg_orthogonalize H ..."
    endif

    call prg_timer_start(ortho_timer)
    call prg_orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
         lt%threshold,lt%bml_type,lt%verbose)
    call prg_timer_stop(ortho_timer)
#ifdef DO_MPI_BLOCK
    call prg_allGatherParallel(orthoh_bml)
#endif

    if (printRank() .eq. 1) then
      call bml_print_matrix("orthoh_bml",orthoh_bml,0,6,0,6)
      call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
      call bml_print_matrix("zmat_bml",zmat_bml,0,6,0,6)
    endif

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
        allocate(xadj(norb+1))
        allocate(adjncy(norb*norb))
        call bml_adjacency(g_bml, xadj, adjncy, 1)
        nnodes = norb
      else
        !write(*,*) "METIS with Atom graph elements is not available"
        !stop
        allocate(xadj(nnodes+1))
        allocate(adjncy(nnodes*nnodes))
        call bml_adjacency_group(g_bml, hnode, nnodes, xadj, adjncy, 1)
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

    !> Calculate gershgorin bounds
    call bml_gershgorin(orthoh_bml, gbnd)
    gp%mineval = gbnd(1)
    gp%maxeval = gbnd(2)
    if (printRank() .eq. 1) then
      write(*,*) "Gershgorin: mineval = ", gbnd(1), " maxeval = ", gbnd(2)
      write(*,*)
    endif

    !! Calculate SP2 sequence
    call prg_sp2sequence(gp%pp, gp%maxIter, gp%mineval, gp%maxeval, ehomo, elumo, &
         gsp2%errlimit)
    if (printRank() .eq. 1) then
      write(*,*)
      write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter
      write(*,*)
    endif

    !> Now use the graph-based SP2 algorithm to get the orthogonal Density
    ! matrix.
    call prg_timer_start(graphsp2_timer)

    call prg_subgraphSP2Loop(orthoh_bml, g_bml, orthop_bml, gp, lt%threshold)
    !    call prg_sp2_alg1_seq(orthoh_bml,orthop_bml,lt%threshold, gp%pp, gp%maxIter, gp%vv)

    call prg_timer_stop(graphsp2_timer)
#ifdef DO_MPI_BLOCK
    call prg_allGatherParallel(orthop_bml)
#endif

    if (printRank() .eq. 1) then
      call bml_print_matrix("gsp2 orthop_bml",orthop_bml,0,6,0,6)
    endif

    !> Calculate Homo-Lumo gap
    call prg_homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, egap)

    if (printRank() .eq. 1) then
      write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
           " egap = ", egap
      write(*,*)
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

    !> Get the system charges.
    call prg_get_charges(rho_bml,over_bml,hindex,sy%net_charge,tb%numel,&
         sy%spindex,lt%mdim,lt%threshold)

    !if (printRank() .eq. 1) then
    write(*,*)"Total charge", sum(sy%net_charge(:)),size(sy%net_charge,dim=1)
    !endif

    call prg_qmixer(sy%net_charge,charges_old,dqin,&
         dqout,scferror,i,lt%pulaycoeff,lt%mpulay,lt%verbose)

    charges_old = sy%net_charge

    if (printRank() .eq. 1) then
      write(*,*)"System charges:"
      do j=1,4
        write(*,*)sy%symbol(j),sy%net_charge(j)
      enddo
    endif
    if(scferror.lt.lt%scftol.and.i.gt.5) then
      if (printRank() .eq. 1) then
        write(*,*)"SCF converged within",i,"steps ..."
        write(*,*)"SCF error =",scferror
      endif
      exit
    endif

  enddo
  !> End of SCF loop.

  allocate(eigenvals(norb))
  call prg_get_eigenvalues(ham_bml,eigenvals,lt%verbose)
  call prg_write_tdos(eigenvals, 0.01_dp, 1000, -25.0_dp, 20.0_dp, "dos.dos")

  deallocate(gbnd)
  deallocate(hnode)

  if (allocated(part) .eqv. .true.) deallocate(part)
  if (allocated(core_count) .eqv. .true.) deallocate(core_count)

  call bml_deallocate(ham_bml)
  call bml_deallocate(ham0_bml)
  call bml_deallocate(orthoh_bml)
  call bml_deallocate(g_bml)
  call bml_deallocate(rho_bml)
  call bml_deallocate(over_bml)
  call bml_deallocate(zmat_bml)

  !> Progress is done
  call prg_progress_shutdown()

end program gploop
