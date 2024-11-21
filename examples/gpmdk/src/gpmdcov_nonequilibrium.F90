module gpmdcov_nonequilibrium_mod

  use bml
  use latteparser_latte_mod
  use gpmdcov_vars
  use prg_quantumdynamics_mod

  public :: gpmdcov_get_currents, gpmdcov_apply_voltage

  real(dp),parameter :: hbar = 0.65821188926_dp

contains

  !> Computing currents beween atoms 
  !! \brief This will compute the commutator between the Hamiltonian and the 
  !! Desity matrix. It will also multiply by a factor -1/hbar to get 
  !! units of currents. The routine will the compute the effective 
  !! current that would circulate between each bond given the new 
  !! Hamiltonian at time t and the previus Density matrix at time t - 1
  !! \param h_bml The hamiltonian of the core+halo subsystem
  !! \param p_bml The density matrix for the core+halo ar previous time
  !! \param z_bml The overlap matriz at time t
  !! \param subsyhindex The hindex of the subsystem. subsyhindex(1,i) = 
  !! Initial orbital index for atom i
  !! nats_core Number of atoms of the core subsystem
  !! norb_core Number of orbitals for the core subsystem 
  !! core_halo_index Indices that maps subsystem to full system. atom "i" in the 
  !! subsystem will have "core_halo_index(i) + 1" index in the full system 
  !! symbols Symbols list for the full system
  !! currthr Threshold to plot the current
  !! 
  subroutine gpmdcov_get_currents(h_bml,p_bml,z_bml,subsyhindex,nats_core,&
      &norb_core,core_halo_index,symbols,currthr)

    implicit none

    character(2), allocatable        ::  symbols(:)
    integer                          ::  ati, atj, ii, jj
    integer                          ::  nats_core, norbH, norbP, norb_core
    integer, allocatable             ::  core_halo_index(:)
    integer, allocatable, intent(in)  ::  subsyhindex(:,:)
    real(dp)                         ::  currthr, factor
    real(dp), allocatable            ::  aux_dense(:,:), com_dense(:,:), currents(:,:)
    real(dp), allocatable            ::  hCore_dense(:,:), h_dense(:,:), pCore_dense(:,:), p_dense(:,:)
    real(dp), allocatable            ::  sinv_dense(:,:), zCore_dense(:,:), z_dense(:,:)
    type(bml_matrix_t)               ::  h_bml, p_bml, zCore_bml, z_bml


    !Since p is from step MDstep - 1 it could have a 
    !different number of orbitals compared to h
    norbH = bml_get_N(h_bml) 
    norbP = bml_get_N(p_bml) 

    allocate(h_dense(norbH,norbH))
    allocate(p_dense(norbP,norbP))
    allocate(z_dense(norbH,norbH))
    call bml_export_to_dense(h_bml,h_dense)
    call bml_export_to_dense(p_bml,p_dense)
    call bml_export_to_dense(z_bml,z_dense)
    
    allocate(hCore_dense(norb_core,norb_core))
    hCore_dense(1:norb_core,1:norb_core) = h_dense(1:norb_core,1:norb_core)
    deallocate(h_dense)

    allocate(pCore_dense(norb_core,norb_core))
    pCore_dense(1:norb_core,1:norb_core) = p_dense(1:norb_core,1:norb_core)
    deallocate(p_dense)

    allocate(zCore_dense(norb_core,norb_core))
    zCore_dense(1:norb_core,1:norb_core) = z_dense(1:norb_core,1:norb_core)
    deallocate(z_dense)

    allocate(com_dense(norb_core,norb_core))
    allocate(aux_dense(norb_core,norb_core))

    factor =  -1.0_dp/hbar

    sinv_dense = matmul(zCore_dense,zCore_dense)    
    com_dense = matmul(hCore_dense,pCore_dense)
    com_dense = matmul(sinv_dense,com_dense)
    aux_dense = matmul(pCore_dense,hCore_dense)
    com_dense = com_dense - matmul(aux_dense,sinv_dense)
    allocate(currents(nats_core,nats_core))

    currents = 0.0_dp
    do ati = 1,nats_core-1
      do atj = ati+1,nats_core
        currents(ati,atj)=0.0_dp
          do ii = subsyhindex(1,ati),subsyhindex(2,ati)
            do jj = subsyhindex(1,atj),subsyhindex(2,atj)
            currents(ati,atj) = currents(ati,atj) + com_dense(ii,jj)
          enddo
        enddo
      enddo
    enddo

    do ii=1,nats_core
      do jj=ii+1,nats_core
        ati = core_halo_index(ii)+1
        atj = core_halo_index(jj)+1
        if(abs(currents(ii,jj)) .gt. currthr)then
          write(*,*)"Currents:",MDstep,ati,atj,symbols(ati),symbols(atj),&
          &Currents(ii,jj)
        endif
      enddo
    enddo

    deallocate(aux_dense, com_dense, currents)
    deallocate(hCore_dense, pCore_dense)
    deallocate(sinv_dense, zCore_dense)

  end subroutine gpmdcov_get_currents


  !> Applying a voltage (a shift in the chemical potential)
  !! \brief This will apply a shift in the chemical potentials of preselected 
  !! atomic sites.
  !! \param fullSystemNats Number of atoms for the full system 
  !! \param subSystemNats Number of atoms fro the subsystem 
  !! \param subsyhindex hindex for the subsystem 
  !! \param core_halo_index atomic indices for the core+halo subsystem
  !! \param h_bml Non-scc hamiltonian for the subsystem
  !! \param s_bml Overlap matrix for the subsystem 
  !! 
  subroutine gpmdcov_apply_voltage(fullSystemNats,subSystemNats,subsyhindex,core_halo_index,h_bml,s_bml)

    implicit none

    integer, allocatable, intent(in)  ::  subsyhindex(:,:)
    character(2), allocatable        ::  symbols(:)
    character(20) :: bml_type
    integer                          ::  ati, atj, ii, jj, nvpoints
    integer                          ::  nats_core, norbH, norbP, norbs, ich, chNats
    integer, allocatable, intent(in)             ::  core_halo_index(:)
    integer, intent(in) :: fullSystemNats, subSystemNats 
    real(dp) :: realtmp
    real(dp), allocatable :: h_dense(:,:),s_dense(:,:),hamv(:),vspsv_dense(:,:)
    type(bml_matrix_t)               ::  h_bml, s_bml

    if(.not. vinit)then 
        open(1,file=gpmdt%voltagef,status='OLD')
        read(1,*)nvpoints
        allocate(voltagev(fullSystemNats))
        do ii = 1,nvpoints
           read(1,*)ati,realtmp
           voltagev(ati) = realtmp
        enddo
        vinit = .true.
        write(*,*)"Applying a shift in the chemical potentials using information in ",gpmdt%voltagef
    endif

    norbs = bml_get_N(h_bml)

    allocate(h_dense(norbs,norbs))
    allocate(s_dense(norbs,norbs))
    allocate(vspsv_dense(norbs,norbs))
    allocate(hamv(norbs))
    hamv = 0.0_dp
    call bml_export_to_dense(h_bml,h_dense)
    call bml_export_to_dense(s_bml,s_dense)
 
    do ich = 1,subSystemNats 
        ati = core_halo_index(ich)+1 !Get the atom from the ch index 
        do ii = subsyhindex(1,ich),subsyhindex(2,ich)
                hamv(ii) = voltagev(ati)  
        enddo
    enddo

    do ii = 1,norbs
        vspsv_dense(ii,:) = 0.5_dp*(hamv(ii)*s_dense(ii,:) + s_dense(i,:)*hamv(:))
    enddo

    deallocate(s_dense,hamv)

    h_dense = h_dense + vspsv_dense 

    deallocate(vspsv_dense)  

    bml_type = bml_get_type(h_bml)

    call bml_import_from_dense(bml_type, h_dense, h_bml, 0.0_dp, norbs)

    deallocate(h_dense)

  end subroutine gpmdcov_apply_voltage


end module gpmdcov_nonequilibrium_mod
