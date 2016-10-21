!> High-level program to perform SFC cycles in Extended Huckel Hamiltonian.
!! This program takes coordinates in xyz or pdb format and extracts information 
!! about the system. 
program main

  !BML lib.
  use bml 
  !LATTE lib modes.
  use system_mod
  use ptable_mod
  use latteparser_latte_mod
  use huckel_latte_mod
  use tbparams_latte_mod 
  use ham_latte_mod
  use coulomb_latte_mod
  use charges_mod
  use sp2_mod
  use sp2parser_mod
  use initmatrices_mod
  use genz_mod
  use nonortho_mod
  use pulaymixer_mod
  use density_mod
  use timer_mod
  
  implicit none     

  integer, parameter     ::  dp = kind(1.0d0)
  integer                ::  i, nel, norb, j
  integer, allocatable   ::  hindex(:,:)
  character(3), allocatable         ::  intKind(:)  
  real(dp)               ::  bndfil, scferror
  real(dp), allocatable  ::  charges_old(:), coul_forces_k(:,:), coul_forces_r(:,:), coul_pot_k(:)
  real(dp), allocatable  ::  coul_pot_r(:), dqin(:,:), dqout(:,:)
  type(bml_matrix_t)     ::  ham0_bml, ham_bml, orthoh_bml, orthop_bml
  type(bml_matrix_t)     ::  over_bml, rho_bml, zmat_bml
  type(latte_type)       ::  lt
  type(sp2data_type)     ::  sp2
  type(system_type)      ::  sy
  type(tbparams_type)    ::  tb
  type(intpairs_type), allocatable  ::  intPairsH(:,:), intPairsS(:,:)  
  real(dp), allocatable             ::  onsitesH(:,:), onsitesS(:,:)
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)  
  
  !> Initialize timer
!   call timer_init()
  
  !> Parsing input file. This file contains all the variables needed to 
  !  run the scf including the sp2 (solver) variables. lt is "latte_type" structure 
  !  containing all the variables. 
  !  file://~/progress/build/doc/html/structlatteparser__latte__mod_1_1latte__type.html
  call parse_latte(lt,"input.in") 

  !> SP2 algorithm to get the orthogonal Density matrix (orthop).   
  call parse_sp2(sp2,"input.in") 
  
  !> Parsing system coordinates. This reads the coords.pdb file to get the position of every 
  !  atom in the system. sy is the "system_type" structure containing all the variables.
  !  file://~/progress/build/doc/html/classsystem__latte__mod.html
  call parse_system(sy,"coords","pdb") 

  !> Get Huckel hamiltonian. Computes the Extended Huckel Hamiltonian from the 
  !  atom coordinates. The main inputs are the huckelTBparams and the system coordinate (sy%coordinate)
  !  The main outputs are Hamiltonian (ham_bml) and Overlap (over_bml) matrices.
  call get_hshuckel(ham_bml,over_bml,sy%coordinate,sy%spindex,sy%spatnum,&
    "../../../huckelTBparams",lt%bml_type,lt%mdim,lt%threshold&
    ,tb%nsp,tb%splist,tb%basis,tb%numel,tb%onsite_energ,&
    tb%norbi,tb%hubbardu)    

  !> LATTE Hamiltonian
!   call load_latteTBparams(tb,sy%splist,"/home/christian/progress/latteTBparams") 
!   call get_hindex(sy%spindex,tb%norbi,hindex,norb)        
!   call load_bintTBparamsH(sy%splist,tb%onsite_energ,&
!     typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,"/home/christian/progress/latteTBparams") 
!   call write_bintTBparamsH(typeA,typeB,&
!       intKind,intPairsH,intPairsS,"mybondints.nonorth")  
!   call init_hsmat(ham_bml,over_bml,lt%bml_type,lt%mdim,norb)
!   call get_hsmat(ham_bml,over_bml,sy%coordinate,sy%lattice_vector,sy%spindex,&
!     tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)
    
    
  !> Get the mapping of the Hamiltonian index with the atom index 
  !  hindex(1,i)=starting Hindex for atom i.
  !  hindex(2,i)=final Hindex for atom i.
  !  file://~/progress/build/doc/html/ham__latte__mod_8F90_source.html
   call get_hindex(sy%spindex,tb%norbi,hindex,norb)        

  call bml_print_matrix("ham0_bml",ham_bml,0,6,0,6)

  !> Get occupation based on last shell population. 
  !  WARNING: This could change depending on the TB method being used.
  nel = sum(element_numel(sy%atomic_number(:)),&
    size(sy%atomic_number,dim=1))
  bndfil = nel/(2.0_dp*norb)
  write(*,*)"nel,norb,bndfil",nel,norb,bndfil

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>  First Charge computation 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the density matrix (rho_bml) and inverse overlap factor (zmat_bml).
  call init_pzmat(rho_bml,zmat_bml,lt%bml_type,lt%mdim,norb)

  !> Get the Inverse square root overlap matrix.
  call buildzdiag(over_bml,zmat_bml,lt%threshold,lt%mdim,lt%bml_type)

  !> Initialize the orthogonal versions of ham and rho.
  call init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)

  !> Orthogonalize ham.
  call orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
    lt%threshold,lt%bml_type,lt%verbose)

  ! SP2 algorithm.
    call sp2_alg2(orthoh_bml,orthop_bml,lt%threshold,bndfil,sp2%minsp2iter,& 
      sp2%maxsp2iter,sp2%sp2conv,sp2%sp2tol)
!   call build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
    
  !> Deorthogonalize rho.       
  call deorthogonalize(orthop_bml,zmat_bml,rho_bml,&
    lt%threshold,lt%bml_type,lt%verbose)

  call bml_print_matrix("rho_bml",rho_bml,0,6,0,6)       

  !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing 
  !  the charges.
  call prg_get_charges(rho_bml, over_bml, hindex, sy%net_charge, tb%numel, sy%spindex, lt%mdim, lt%threshold)
  charges_old = sy%net_charge

  write(*,*)"Total charges =", sum(sy%net_charge(:),size(sy%net_charge,dim=1))
  
  write(*,*)"System charges:"
  do i=1,sy%nats
    write(*,*)sy%symbol(i),sy%net_charge(i)
  enddo

!   write(*,*)nel,"1404",norb
! stop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>  SCF loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Save actual hamiltonian as the non-scf Hamiltonian (H0)
  call bml_copy_new(ham_bml,ham0_bml)

  !> Get the reciprocal vector (this is needed to compute the Coulombic interactions)
  call get_recip_vects(sy%lattice_vector,sy%recip_vector,sy%volr,sy%volk)
!   write(*,*)sy%lattice_vector,sy%recip_vector
!    stop

  !> Initialize timer
  call timer_init()
  
  !> Beginning of the SCF loop.
  do i=1,lt%maxscf    

    write(*,*)"SCF iter", i

    !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
    write(*,*)"In real Coul ..."
    call get_ewald_real(sy%spindex,sy%splist,sy%coordinate&
      ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
      sy%volr,lt%coul_acc,coul_forces_r,coul_pot_r);

    !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
    write(*,*)"In recip Coul ..."    
!     if(.not.allocated(coul_pot_k))allocate(coul_pot_k(sy%nats))
!     coul_pot_k = 0.0_dp
    call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
      ,sy%net_charge,tb%hubbardu,sy%lattice_vector,&
      sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);  

    !> Get the scf hamiltonian. The outputs is ham_bml.
    write(*,*)"in get_hscf ..."
    call get_hscf(ham0_bml,over_bml,ham_bml,sy%spindex,hindex,tb%hubbardu,sy%net_charge,&
      coul_pot_r,coul_pot_k,lt%mdim,lt%threshold)

    !> Initialize the orthogonal versions of ham and rho.
    call init_ortho(orthoh_bml,orthop_bml,lt%bml_type,lt%mdim,norb)

    !> Orthogonalize the Hamiltonian
    write(*,*)"in orthogonalize H ..."
    call orthogonalize(ham_bml,zmat_bml,orthoh_bml,&
      lt%threshold,lt%bml_type,lt%verbose)

    call bml_print_matrix("orthoh_bml",orthoh_bml,0,6,0,6)
    call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
    call bml_print_matrix("zmat_bml",zmat_bml,0,6,0,6)    

    !> Solvers    
    !> SP2 algorithm to get the orthogonal Density matrix (orthop).
    write(*,*)sp2_timer

    call timer_start(dyn_timer,"mytag")
    call sp2_alg2(orthoh_bml,orthop_bml,lt%threshold,bndfil,sp2%minsp2iter,& 
      sp2%maxsp2iter,sp2%sp2conv,sp2%sp2tol)
    call timer_stop(dyn_timer,1)      
    
!     call build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
    
    call bml_print_matrix("sp2 orthop_bml",orthop_bml,0,6,0,6)

    !> Deorthogonalize orthop_bml to get the density matrix rho_bml.
    call deorthogonalize(orthop_bml,zmat_bml,rho_bml,&
      lt%threshold,lt%bml_type,lt%verbose)

    call bml_deallocate(orthop_bml)

    !> Get the system charges.
    call prg_get_charges(rho_bml,over_bml,hindex,sy%net_charge,tb%numel,&
      sy%spindex,lt%mdim,lt%threshold)

    write(*,*)"Total charge", sum(sy%net_charge(:)),size(sy%net_charge,dim=1)

    call qmixer(sy%net_charge,charges_old,dqin,&
      dqout,scferror,i,lt%pulaycoeff,lt%mpulay,lt%verbose)
      
!     call linearmixer(sy%net_charge,charges_old,&
!       scferror,lt%pulaycoeff,lt%verbose)
      
    write(*,*)"System charges:"
    do j=1,4
      write(*,*)sy%symbol(j),sy%net_charge(j)
    enddo
      
    if(scferror.lt.lt%scftol) then 
      write(*,*)"SCF converged within",i,"steps ..."
      write(*,*)"SCF error =",scferror
      exit
    endif 

  enddo
  !> End of SCF loop.

  end 

