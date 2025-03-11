  !> Finalize gpmd
  subroutine gpmdcov_Finalize()

    use gpmdcov_vars
   
    implicit none 
    integer :: io_max, io, m, isiostat
    integer :: k, l
    logical :: isopened
    
    if(bml_allocated(ham_bml))call bml_deallocate(ham_bml)
    if(bml_allocated(ham0_bml))call bml_deallocate(ham0_bml)
    if(bml_allocated(g_bml))call bml_deallocate(g_bml)
    if(bml_allocated(rho_bml))call bml_deallocate(rho_bml)
    if(bml_allocated(over_bml))call bml_deallocate(over_bml)
    if(bml_allocated(zmat_bml))call bml_deallocate(zmat_bml)
   
    if(allocated(ppot))then
      do k=1,sy%nsp
        do l=1,sy%nsp
          deallocate(ppot(k,l)%potparams)
        enddo
      enddo
      deallocate(ppot)
    endif

    if(allocated(gpat%sgraph))call  prg_destroyGraphPartitioning(gpat)
    if(allocated(sy%coordinate))call prg_destroy_system(sy)
    if(allocated(sy%symbol)) deallocate(sy%symbol)
    if(allocated(sy%atomic_number)) deallocate(sy%atomic_number)
    if(allocated(sy%coordinate)) deallocate(sy%coordinate)
    !if(allocated(sy%velocity)) deallocate(sy%velocity)
    !if(allocated(sy%force)) deallocate(sy%force)
    if(allocated(sy%net_charge)) deallocate(sy%net_charge)
    if(allocated(sy%mass)) deallocate(sy%mass)
    if(allocated(sy%lattice_vector)) deallocate(sy%lattice_vector)
    !if(allocated(sy%recip_vector)) deallocate(sy%recip_vector)
    if(allocated(sy%spindex)) deallocate(sy%spindex)
    if(allocated(sy%splist)) deallocate(sy%splist)
    if(allocated(sy%spatnum)) deallocate(sy%spatnum)
    if(allocated(sy%spmass)) deallocate(sy%spmass)
    !if(allocated(sy%userdef)) deallocate(sy%userdef)
    !if(allocated(sy%resindex)) deallocate(sy%resindex)
    !if(allocated(sy%resname)) deallocate(sy%resname)
    !if(allocated(sy%atomname)) deallocate(sy%atomname)

    if(allocated(nl%nntype))call gpmdcov_destroy_nlist(nl,lt%verbose)

    if(allocated(syprt))then
      do ipt=1,gpat%TotalParts
        call prg_destroy_estr(syprt(ipt)%estr)
      enddo
      do ipt=1,gpat%TotalParts
        call prg_destroy_subsystems(syprt(ipt),lt%verbose)
      enddo
      deallocate(syprt)
    endif

    !call prg_progress_shutdown()
    
  if(allocated(TypeA))deallocate(TypeA)
  if(allocated(TypeB))deallocate(TypeB)
  if(allocated(intKind))deallocate(intKind)
  if(allocated(hindex))deallocate(hindex)
  if(allocated(hnode))deallocate(hnode)
  if(allocated( TypeA ))deallocate( TypeA )
  if(allocated( TypeB ))deallocate( TypeB )
  if(allocated( intKind ))deallocate( intKind )
  if(allocated( CH_count))deallocate( CH_count )
  if(allocated( CH_count_cov))deallocate( CH_count_cov )
  if(allocated( Halo_count ))deallocate( Halo_count )
  if(allocated( Halo_count_cov ))deallocate( Halo_count_cov )
  if(allocated( PartsInRankI ))deallocate( PartsInRankI )
  if(allocated( adjncy ))deallocate( adjncy )
  if(allocated( adjncy_cov ))deallocate( adjncy_cov )
  if(allocated( core_count ))deallocate( core_count )
  if(allocated( core_count_cov ))deallocate( core_count_cov)
  if(allocated( displ ))deallocate( displ)
  if(allocated( hindex ))deallocate( hindex)
  if(allocated( hnode ))deallocate( hnode )
  if(allocated( norbsInEachCH ))deallocate( norbsInEachCH )
  if(allocated( norbsInEachCHAtRank ))deallocate( norbsInEachCHAtRank )
  if(allocated( npartsVect ))deallocate( npartsVect )
  if(allocated( norbsInEachRank ))deallocate( norbsInEachRank )
  if(allocated( part ))deallocate( part )
  if(allocated( part_cov ))deallocate( part_cov )
  if(allocated( partsInEachRank ))deallocate( partsInEachRank )
  if(allocated( partsSizes))deallocate( partsSizes)
  if(allocated( reshuffle ))deallocate( reshuffle )
  if(allocated( thisPartSize))deallocate( thisPartSize )
  if(allocated( vectorint ))deallocate( vectorint )
  if(allocated( xadj ))deallocate( xadj )
  if(allocated( xadj_cov ))deallocate( xadj_cov )
  if(allocated( GFPUL ))deallocate( GFPUL )
  if(allocated( GFSCOUL ))deallocate( GFSCOUL )
  if(allocated( Ker ))deallocate( Ker )
  if(allocated( PairForces ))deallocate( PairForces )
  if(allocated( DispForces ))deallocate( DispForces )
  if(allocated( SKForce ))deallocate( SKForce )
  if(allocated( VX ))deallocate( VX)
  if(allocated( VY ))deallocate( VY)
  if(allocated( VZ ))deallocate( VZ)
  if(allocated( acceprat ))deallocate( acceprat )
  if(allocated( auxcharge ))deallocate( auxcharge )
  if(allocated( auxcharge1 ))deallocate( auxcharge1 )
  if(allocated( charges_old ))deallocate( charges_old )
  if(allocated( collectedforce ))deallocate( collectedforce )
  if(allocated( coul_forces ))deallocate( coul_forces )
  if(allocated( coul_forces_k ))deallocate( coul_forces_k )
  if(allocated( coul_forces_r ))deallocate( coul_forces_r )
  if(allocated( coul_pot_k ))deallocate( coul_pot_k )
  if(allocated( coul_pot_r ))deallocate( coul_pot_r )
  if(allocated( dqin ))deallocate( dqin )
  if(allocated( dqout ))deallocate( dqout  )
  if(allocated( dvals ))deallocate( dvals  )
  if(allocated( dvalsAll ))deallocate( dvalsAll  )
  if(allocated( dvalsInRank ))deallocate( dvalsInRank  )
  if(allocated( eigenvals ))deallocate( eigenvals  )
  if(allocated( eigenvalues ))deallocate( eigenvalues  )
  if(allocated( evals ))deallocate( evals  )
  if(allocated( evalsAll ))deallocate( evalsAll  )
  if(allocated( evalsInRank ))deallocate( evalsInRank  )
  if(allocated( fvals ))deallocate( fvals  )
  if(allocated( fvalsAll ))deallocate( fvalsAll  )
  if(allocated( fvalsInRank ))deallocate( fvalsInRank  )
  if(allocated( g_dense ))deallocate( g_dense  )
  if(allocated( gbnd ))deallocate( gbnd )
  if(allocated( n ))deallocate( n )
  if(allocated( n_0 ))deallocate( n_0 )
  if(allocated( n_1 ))deallocate( n_1 )
  if(allocated( n_2 ))deallocate( n_2 )
  if(allocated( n_2 ))deallocate( n_2 )
  if(allocated( n_3 ))deallocate( n_3 )
  if(allocated( n_4 ))deallocate( n_4 )
  if(allocated( n_5 ))deallocate( n_5 )
  if(allocated( onsitesH ))deallocate( onsitesH )
  if(allocated( onsitesS ))deallocate( onsitesS )
  if(allocated( origin ))deallocate( origin )
  if(allocated( rhoat ))deallocate( rhoat )
  if(allocated( row ))deallocate( row )
  if(allocated( row1 ))deallocate( row1 )
  if(allocated( smdForce ))deallocate( smdForce )
  if(allocated(tb%splist))deallocate(tb%splist)
  if(allocated(tb%numel))deallocate(tb%numel)
  if(allocated(tb%basis))deallocate(tb%basis)
  if(allocated(tb%onsite_energ))deallocate(tb%onsite_energ)
  if(allocated(tb%mass))deallocate(tb%mass)
  if(allocated(tb%HubbardU))deallocate(tb%HubbardU)
  if(allocated(tb%w))deallocate(tb%w)
  if(allocated(tb%norbi))deallocate(tb%norbi)
  if(allocated(intPairsH))deallocate(intPairsH)
  if(allocated(intPairsS))deallocate(intPairsS)
  
  if(bml_allocated(  aux_bml  ))call bml_deallocate(   aux_bml      ) 
  if(bml_allocated(  dH0x_bml ))call bml_deallocate(   dH0x_bml     ) 
  if(bml_allocated(  dH0y_bml ))call bml_deallocate(   dH0y_bml     ) 
  if(bml_allocated(  dH0z_bml ))call bml_deallocate(   dH0z_bml     ) 
  if(bml_allocated(  dSx_bml  ))call bml_deallocate(   dSx_bml      ) 
  if(bml_allocated(  dSy_bml  ))call bml_deallocate(   dSy_bml      ) 
  if(bml_allocated(  dSz_bml  ))call bml_deallocate(   dSz_bml      ) 
  if(bml_allocated(  g_bml    ))call bml_deallocate(   g_bml        ) 
  if(bml_allocated(  ham0_bml ))call bml_deallocate(   ham0_bml     ) 
  if(bml_allocated(  ham_bml  ))call bml_deallocate(   ham_bml      ) 
  if(bml_allocated(  over_bml ))call bml_deallocate(   over_bml     ) 
  if(bml_allocated(  rho_bml  ))call bml_deallocate(   rho_bml      ) 
  if(bml_allocated(  rhoat_bml))call bml_deallocate(   rhoat_bml   )  
  if(bml_allocated(  rhoh_bml ))call bml_deallocate(   rhoh_bml     ) 
  if(bml_allocated(  zmat_bml ))call bml_deallocate(   zmat_bml     ) 
  if(bml_allocated(  gch_bml  ))call bml_deallocate(   gch_bml      ) 
  if(bml_allocated(  copy_g_bml))call bml_deallocate( copy_g_bml   )  
  if(bml_allocated(  gcov_bml ))call bml_deallocate(   gcov_bml     ) 
  if(bml_allocated(  aux1_bml ))call bml_deallocate(   aux1_bml     ) 
  if(bml_allocated(  ZK1_bml  ))call bml_deallocate(   ZK1_bml      ) 
  if(bml_allocated(  ZK2_bml  ))call bml_deallocate(   ZK2_bml      ) 
  if(bml_allocated(  ZK3_bml  ))call bml_deallocate(   ZK3_bml      ) 
  if(bml_allocated(  ZK4_bml  ))call bml_deallocate(   ZK4_bml      ) 
  if(bml_allocated(  ZK5_bml  ))call bml_deallocate(   ZK5_bml      ) 
  if(bml_allocated(  ZK6_bml  ))call bml_deallocate(   ZK6_bml      ) 


    io_max = 100
    do io = io_max,1,-1
      inquire (unit=io, opened=isopened, iostat=isiostat)
      if(isopened)then
        if((io .ne. 5) .and. (io .ne. 6))then
       write(*,*)"io",lib_init,io,isopened
       !call flush(io)
       !close(io)
     endif
   endif
    end do
 


  end subroutine gpmdcov_Finalize

