!> Module to compute the density matrix response and related quantities.
!! \ingroup PROGRESS
!! \todo Add the response scf
!! \todo Change name response_SP2 to dm_prt_response
!! \todo Change name response_rs to rs_prt_response
!! More information about the theory can be found at \cite Niklasson2005 and Niklasson2015
module prg_response_mod

  use bml
  use prg_kernelparser_mod
  use prg_densitymatrix_mod
  use prg_openfiles_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950_dp

  type, public :: RespData_type
     character(20) :: RespMode
     character(20) :: TypeOfPert
     character(20) :: BmlType
     integer :: Mdim
     real(dp) :: NumThresh
     logical :: ComputeDipole
     logical :: GetResponse
     real(dp) :: FieldIntensity
     real(dp) :: field(3)
  end type RespData_type

  public :: prg_pert_from_file, prg_compute_dipole, prg_pert_constant_field
  public :: prg_compute_response_RS, prg_parse_response, prg_compute_response_SP2
  public :: prg_compute_response_FD, prg_compute_polarizability, prg_write_dipole_tcl
  public :: prg_project_response, prg_pert_sin_pot, prg_pert_cos_pot

contains

  !> The parser for the calculation of the DM response.
  !! \param RespData Response data type.
  !! \param filename Name of the file to parse.
  !!
  subroutine prg_parse_response(RespData, filename)
    implicit none
    type(RespData_type) :: RespData
    integer, parameter  :: nkey_char = 3, nkey_int = 4, nkey_re = 5, nkey_log = 2
    character(len=*)    :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'TypeOfPert=', 'BMLType=','RespMode=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'None', 'Dense','RayleighSchrodinger']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: 'Mdim=' &
         , 'p2=', 'p3=', 'p4=']
    integer :: valvector_int(nkey_int) = (/0 &
         , 0, 0, 1 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Fieldx=','Fieldy=','Fieldz=' &
         ,'FieldIntensity=','NumThresh=' ]
    real(dp) :: valvector_re(nkey_re) = (/0.0,0.0,0.0&
         ,0.001,0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'ComputeDipole=', 'GetResponse=']
    logical :: valvector_log(nkey_log) = (/ .false., .false./)
    !End of the library

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'RESPONSE{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    RespData%TypeofPert = valvector_char(1)

    if(valvector_char(2) == "Dense")then
       RespData%BmlType = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
       RespData%BmlType = BML_MATRIX_ELLPACK
    endif

    RespData%RespMode = valvector_char(3)

    !Reals
    RespData%Field(1) = valvector_re(1)
    RespData%Field(2) = valvector_re(2)
    RespData%Field(3) = valvector_re(3)
    RespData%FieldIntensity = valvector_re(4)
    RespData%NumThresh = valvector_re(5)

    !Logicals
    RespData%ComputeDipole = valvector_log(1)
    RespData%GetResponse = valvector_log(2)

    !Integers
    RespData%Mdim = valvector_int(1)

  end subroutine prg_parse_response

  !> To compute the dipole moment of the system.
  !! The units of the dipole moment are determined by the units of
  !! the coordinates and charges that are given.
  !! \param charges Charges on each atomic position.
  !! \param coordinate Coordinates of the atoms.
  !! \param nats Number of atoms.
  !! \param dipoleMoment Dipole moment vector.
  !! \param factor Unit conversion factor (use 1.0 is no conversion is required).
  !! \param verbose To give different verbosity levels.
  !! If coordinates are in \f$ \AA \f$ and charges are in fractions of electron, then
  !! transformation ea2debye form LATTE lib can be used to change units to Debye.
  !!
  subroutine prg_compute_dipole(charges,coordinate,dipoleMoment,factor,verbose)

    integer                  ::  cont, i, nats, verbose
    real(dp), intent(in)     ::  charges(:), coordinate(:,:), factor
    real(dp), intent(inout)  ::  dipoleMoment(3)

    if(verbose.gt.1)write(*,*)"Computing dipole moment ..."

    cont=0

    nats = size(charges,dim=1)

    dipoleMoment = 0.0_dp

    do i=1,nats
       dipoleMoment(1)=dipoleMoment(1)+coordinate(1,i)*charges(i)
       dipoleMoment(2)=dipoleMoment(2)+coordinate(2,i)*charges(i)
       dipoleMoment(3)=dipoleMoment(3)+coordinate(3,i)*charges(i)
    enddo

    dipoleMoment=factor*dipoleMoment

    write(*,*)"dipoleMomentx=",dipoleMoment(1)
    write(*,*)"dipoleMomenty=",dipoleMoment(2)
    write(*,*)"dipoleMomentz=",dipoleMoment(3)

  end subroutine prg_compute_dipole

  !> To visualize a dipole moment using VMD.
  !! This will prg_generate a .tcl script that could be run using VMD
  !! To visualize with VMD:
  !!    $ vmd -e dipole.tcl
  !!
  !! \param dipoleMoment Dipole moment vector.
  !! \param file PDB/XYZ file to load for visualization.
  !! \param factor Arbitrary scale for visualization.
  !! \param verbose To give different verbosity levels.
  !!
  subroutine prg_write_dipole_tcl(dipoleMoment,file,factor,verbose)
    implicit none
    integer                  ::  verbose, io_unit
    real(dp), intent(in)     ::  factor
    real(dp), intent(in)     ::  dipoleMoment(3)
    character(*), intent(in) ::  file

    call prg_open_file(io_unit, "dipole.tcl")

    write(io_unit,*)"mol new ",file
    write(io_unit,*)"set molid 0"

    write(io_unit,*)"proc vmd_draw_arrow {mol start end} {"
    write(io_unit,*)" set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]"
    write(io_unit,*)" graphics $mol cylinder $start $middle radius 0.15"
    write(io_unit,*)" graphics $mol cone $middle $end radius 0.25"
    write(io_unit,*)"}"

    write(io_unit,*)"draw arrow  {0 0 0} {",&
         factor*dipoleMoment(1), factor*dipoleMoment(2), factor*dipoleMoment(3),"}"

    close(io_unit)

  end subroutine prg_write_dipole_tcl

  !> To compute the polarizability of the system.
  !! The units of the directional polarizability are determined by the units of
  !! the perturbation and Hamiltonian. This equation can be found in \cite Weber2005 equation 4a.
  !! Note that in equation 4a of the reference there is a 2 that account for the double occupancy
  !! which is not present in this case cause the density matrix construction is done by taking the occupancy
  !! into account.
  !! \param charges Charges on each atomic position.
  !! \param coordinate Coordinates of the atoms.
  !! \param nats Number of atoms.
  !! \param dipoleMoment Dipole moment vector.
  !! \param factor Unit conversion factor (use 1.0 is no conversion is required).
  !! \param verbose To give different verbosity levels.
  !! If coordinates are in \f$ \AA \f$ and charges are in fractions of electron, then
  !! transformation ea2debye form LATTE lib can be used to change units to Debye.
  !!
  subroutine prg_compute_polarizability(rsp_bml,prt_bml,polarizability,factor,verbose)
    implicit none
    integer                  ::  cont, i, nats, verbose
    real(dp), intent(in)     ::  factor
    real(dp), intent(inout)  ::  polarizability
    real(dp)                 ::  dipolez, dipolex, dipoley
    type(bml_matrix_t), intent(in)     ::  prt_bml, rsp_bml
    type(bml_matrix_t) :: aux_bml

    if(verbose.gt.1)write(*,*)"Computing polarizability ..."

    call bml_copy_new(prt_bml,aux_bml)
    call bml_multiply(rsp_bml,prt_bml,aux_bml,1.0_dp,0.0_dp,0.0_dp)

    polarizability = -factor*bml_trace(aux_bml)

    write(*,*)"Polarizability=",polarizability

    call bml_deallocate(aux_bml)

  end subroutine prg_compute_polarizability

  !> Read perturbation from file.
  !! \todo Add read perturbation from file
  !!
  subroutine prg_pert_from_file(prt_bml,norb)
    implicit none
    integer :: Norb, verbose
    type(bml_matrix_t), intent(inout) :: prt_bml
    if(verbose.gt.1) write(*,*)'Reading perturbation V from file ...'

  end subroutine prg_pert_from_file

  !> Computes the first order response density matrix using Rayleigh Schrodinger Perturbation
  !! theory
  !! The transformation hereby performed are:
  !! - \f$ V = C^{\dagger} H^{(1)} C \f$
  !! - \f$ \tilde{V}_{ij} = \frac{V_{ij}}{\epsilon_j - \epsilon_i} \f$, with
  !! \f$ \tilde{V}_{ii} = 0 \, \forall i \f$.
  !! - \f$ C^{(1)} = C \tilde{V} \f$
  !! - And finally: \f$ \rho^{(1)} = Cf(C^{(1)})^{\dagger} + C^{(1)}fC^{\dagger}\f$
  !! \param ham_bml Hamiltonian in bml format (\f$ H^{(0)} \f$).
  !! \param prt_bml Perturbation in bml format (\f$ H^{(1)} \f$).
  !! \param rsp_bml First order response to the perturbation (\f$ \rho^{(1)} \f$).
  !! \param bndfil Filing factor.
  !! \param threshold Threshold value for matrix elements.
  !! \param verbose Different levels of verbosity.
  !! \warning This works only for the prg_orthogonalized form of ham_bml.
  !! \warning The response must be in the prg_orthogonalized form.
  !!
  subroutine prg_compute_response_RS(ham_bml,prt_bml,rsp_bml,lambda&
       ,bndfil,threshold,verbose)

    character(20)                      ::  bml_type
    integer                            ::  i, j, l, norb
    integer                            ::  verbose
    real(dp)                           ::  lambda, nocc
    real(dp), allocatable              ::  aux(:,:), evals(:), row(:)
    real(dp), intent(in)               ::  bndfil, threshold
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, occupation_bml, umat_bml
    type(bml_matrix_t)                 ::  aux2_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml, prt_bml
    type(bml_matrix_t), intent(inout)  ::  rsp_bml

    norb = bml_get_N(ham_bml)

    if(verbose.ge.1) write(*,*)'In prg_compute_response_RS...  Norbs = ', norb

    allocate(evals(norb))

    bml_type = bml_get_type(ham_bml)

    !Ensure dense type to diagonalize
    allocate(aux(norb,norb))
    call bml_export_to_dense(ham_bml,aux)
    call bml_zero_matrix(bml_matrix_dense,bml_element_real,dp,norb,norb,aux_bml)
    call bml_import_from_dense(bml_matrix_dense,aux,aux_bml,threshold,norb)
    deallocate(aux)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,umat_bml)

    !Eigenvectors and eigevalues of H
    call bml_diagonalize(aux_bml,evals,umat_bml)
    call bml_deallocate(aux_bml)

    !Ensure the umat type is in bml_type.
    allocate(aux(norb,norb))
    call bml_export_to_dense(umat_bml,aux)
    call bml_deallocate(umat_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,umat_bml)
    call bml_import_from_dense(bml_type,aux,umat_bml,threshold,norb)
    deallocate(aux)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)

    !Obtaining matrix $\f V = C^{\dagger} H^{(1)} C^{\dagger} $\f
    call bml_transpose(umat_bml,aux_bml) !C^{\dagger}
    !First product C^{\dagger} H^{(1)}
    call bml_multiply(aux_bml,prt_bml,aux1_bml,1.0_dp,0.0_dp,threshold)
    !Seccond product to obtain V (aux_bml)
    call bml_multiply(aux1_bml,umat_bml,aux_bml,1.0_dp,0.0_dp,threshold)

    !\f$ \tilde{V}_{ij} = \frac{V_{ij}}{\epsilon_j - \epsilon_i} \f$
    allocate(row(norb))
    do i=1,norb
       call bml_get_row(aux_bml,i,row)
       do j=1,norb
          if(j.ne.i)then
             row(j) = row(j)/(evals(j)-evals(i))
          else
             row(j) = 0.0_dp
          endif
       enddo
       call bml_set_row(aux_bml,i,row,threshold)
    enddo
    deallocate(row)

    !\f$ C^{(1)} = C \tilde{V} \f$
    call bml_multiply(umat_bml,aux_bml,aux1_bml,1.0_dp,0.0_dp,threshold)

    !Get matrix f
    nocc = norb*bndfil
    write(*,*)"Trace CV",bml_trace(aux_bml), nocc

    do i=1,norb                 !Reusing eigenvalues to apply the theta function.
       if(i.le.nocc) then
          evals(i) = 2.0_dp
       else
          evals(i) = 0.0_dp
       endif
    enddo
    if(abs(nocc - int(nocc)).gt.0.01_dp)then
       evals(int(nocc)+1) = 1.0_dp
    endif

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_set_diagonal(occupation_bml, evals) !eps(i,i) = eps(i)
    deallocate(evals)

    !Cf(C^{(1)})^{t}$
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux2_bml)
    call bml_multiply(umat_bml,occupation_bml,aux_bml,1.0_dp,0.0_dp,threshold)
    call bml_transpose(aux1_bml,aux2_bml)
    call bml_multiply(aux_bml,aux2_bml,rsp_bml,1.0_dp,0.0_dp,threshold)

    !CfC^{(1)}$ + C^{(1)}fC
    call bml_transpose(umat_bml,aux2_bml)
    call bml_multiply(occupation_bml,aux2_bml,aux_bml,1.0_dp,0.0_dp,threshold)

    call bml_deallocate(occupation_bml)
    call bml_deallocate(umat_bml)

    call bml_multiply(aux1_bml,aux_bml,rsp_bml,1.0_dp,1.0_dp,threshold)
    call bml_scale(1.0_dp,rsp_bml)

    write(*,*)"Trace",bml_trace(rsp_bml)

    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_compute_response_RS

  !> Computes the first order response density matrix using finite differences.
  !! The transformation hereby performed are:
  !! - \f$ H^+ = H^{(0)} + \prg_delta H^{(1)} \f$
  !! - \f$ H^- = H^{(0)} - \prg_delta H^{(1)} \f$
  !! - \f$ \rho^+ = f(H^+) \f$
  !! - \f$ \rho^- = f(H^-) \f$
  !! - \f$ \rho^{(1)} =  (\rho^+ - \rho^-)/(2\prg_delta) \f$.
  !! Where f denotes the Fermi function (construction of the density matrix)
  !! \param ham_bml Hamiltonian in bml format (\f$ H^{(0)} \f$).
  !! \param prt_bml Perturbation in bml format (\f$ H^{(1)} \f$).
  !! \param rsp_bml First order response to the perturbation (\f$ \rho^{(1)} \f$).
  !! \param bndfil Filing factor.
  !! \param threshold Threshold value for matrix elements.
  !! \param verbose Different levels of verbosity.
  !! \warning This works only for the prg_orthogonalized form of ham_bml.
  !! \warning The response must be in the prg_orthogonalized form.
  !!
  subroutine prg_compute_response_FD(ham_bml,prt_bml,rsp_bml,prg_delta&
       ,bndfil,threshold,verbose)

    character(20)                      ::  bml_type
    integer                            ::  i, j, l, norb
    integer                            ::  verbose
    real(dp)                           ::  prg_delta, nocc
    real(dp), allocatable              ::  resp(:,:), evals(:), row(:)
    real(dp), intent(in)               ::  bndfil, threshold
    type(bml_matrix_t)                 ::  rhoplus_bml, rhominus_bml, occupation_bml, umat_bml
    type(bml_matrix_t)                 ::  aux_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml, prt_bml
    type(bml_matrix_t), intent(inout)  ::  rsp_bml

    norb = bml_get_N(ham_bml)

    if(verbose.ge.1) write(*,*)'In prg_compute_response_RS...  Norbs = ', norb

    allocate(resp(norb,norb))

    call bml_copy_new(ham_bml,aux_bml)

    bml_type = bml_get_type(ham_bml)

    !! - \f$ H^+ = H^{(1)} + \prg_delta I \f$
    call bml_add_deprecated(1.0_dp, aux_bml, prg_delta, prt_bml, threshold)

    call bml_copy_new(ham_bml,rhoplus_bml)

    !! - \f$ \rho^+ = f(H + H^+) \f$
    call prg_build_density_t0(aux_bml,rhoplus_bml,threshold,bndfil)

    !! - \f$ H^- = H^{(1)} - \prg_delta I \f$
    call bml_add_deprecated(1.0_dp, aux_bml, -2.0_dp*prg_delta, prt_bml, threshold)

    call bml_copy_new(ham_bml,rhominus_bml)

    !! - \f$ \rho^- = f(H + H^-) \f$
    call prg_build_density_t0(aux_bml,rhominus_bml,threshold,bndfil)

    !! - \f$ \rho^{(1)} =  (\rho^+ - \rho^-)/(2\prg_delta) \f$.
    call bml_add_deprecated(1.0_dp, rhoplus_bml, -1.0_dp, rhominus_bml, threshold)
    call bml_scale(1.0_dp/(2.0_dp*prg_delta),rhoplus_bml,rsp_bml)

  end subroutine prg_compute_response_FD

  !> Apply a constant field perturbation through the dipole moment operator
  !! (\f$ \hat{\mu} = e \hat{\textbf{r}} \f$). In the matrix representation, this is:
  !! \f$ H^{(1)} = \lambda \frac{1}{2}(\,S \, e \textbf{r} \cdot \textbf{E} +  \, e \textbf{r} \cdot \textbf{E}S) \f$.
  !! The symmetrization is done in order to preserve the Hermiticity of H.
  !! In this case the whole system will be affected by the field. In a
  !! latter version we will add the possibility of applying this
  !! field to a region of the system. In this implementation \f$ e= 1 \f$
  !! and units can be transformed by using the parameter \f$ \lambda \f$.
  !! \note If the Hamiltonian is already in the prg_orthogonalized form, then parameter over_bml
  !! can be omitted.
  !! \param field Direction of the applied field (\f$ \hat{\textbf{E}} \f$).
  !! \param intensity Intensity of the field (\f$ ||\textbf{E}|| \f$)..
  !! \param coordinate Coordinates of the system (\f$ \textbf{r}\f$).
  !! \param lambda Constant to premultiply the perturbation (\f$ \lambda\f$).
  !! \param prt_bml Perturbation in bml format (\f$ H^{(1)} \f$).
  !! \param threshold Threshold value for bml format matrices.
  !! \param spindex Species index. It gives the species index of a particular atom.
  !! \param norbi Number of orbitals for each atomic site.
  !! \param verbose Different levels of verbosity.
  !! \param over_bml It has to be present for a nonorthogonal representation (\f$ S \f$).
  !!
  subroutine prg_pert_constant_field(field,intensity,coordinate,lambda&
       ,prt_bml,threshold,spindex,norbi,verbose,over_bml)

    integer                                  ::  cont, i, ii, nats
    integer                                  ::  norb
    integer, intent(in)                      ::  spindex(:), norbi(:), verbose
    real(dp)                                 ::  fieldnorm, fieldx, fieldy, fieldz
    real(dp)                                 ::  intensity, lambda, threshold
    real(dp), allocatable                    ::  diag(:)
    real(dp), intent(in)                     ::  coordinate(:,:), field(3)
    type(bml_matrix_t)                       ::  aux_bml
    type(bml_matrix_t), intent(inout)        ::  prt_bml
    type(bml_matrix_t), optional, intent(in)  ::  over_bml

    write(*,*)'Applying a constant field perturbation ...'

    fieldx = field(1)
    fieldy = field(2)
    fieldz = field(3)

    nats=size(spindex,dim=1)
    norb=sum(norbi(spindex(:)))

    write(*,*)nats,norb

    fieldnorm=sqrt(fieldx**2+fieldy**2+fieldz**2) !Get the norm

    fieldx=intensity*fieldx/fieldnorm   !Normalize and add intensity
    fieldy=intensity*fieldy/fieldnorm
    fieldz=intensity*fieldz/fieldnorm

    cont=0;

    allocate(diag(norb))

    do i=1,nats
       do ii=1,norbi(spindex(i))
          cont=cont+1
          diag(cont)=lambda*(fieldx*coordinate(1,i) + &
               fieldy*coordinate(2,i) + fieldz*coordinate(3,i))
       enddo
    enddo

    if(bml_get_N(prt_bml) < 0) stop "bml_pert not allocated"

    call bml_set_diagonal(prt_bml,diag,threshold)

    deallocate(diag)

    if(present(over_bml))then  !(S*V + V*S)/2
       call bml_copy_new(prt_bml,aux_bml)
       call bml_multiply(over_bml,aux_bml,prt_bml,0.5_dp, 0.0_dp,threshold)
       call bml_multiply(aux_bml,over_bml,prt_bml,0.5_dp,1.0_dp,threshold)
       call bml_deallocate(aux_bml)
    endif

  end subroutine prg_pert_constant_field

  !> Apply a sinusoidal length dependent potential
  !! (\f$ \sin(\tilde{\textbf{r}}_x) \f$) where \f$ \textbf{r}_x \f$ is the x coordinate.
  !! The Hamiltonian gets modified as follows:
  !! \f$ H^{(1)} = \frac{1}{2}\lambda (S \sin(\tilde{\textbf{r}}_x) + \sin(\tilde{\textbf{r}}_x) S) \f$.
  !! \f$ \tilde{\textbf{r}}_x = 2\pi(\textbf{r}/l_x) - \pi \f$.
  !! The symmetrization is done in order to preserve the Hermiticity of H.
  !! Units can be transformed by using the parameter \f$ \lambda \f$.
  !! \note If the Hamiltonian is already in the prg_orthogonalized form, then parameter over_bml
  !! can be omitted.
  !! \param direction Direction of the potential gradient (x,y or z).
  !! \param lx Length of the box in x direction.
  !! \param coordinate Coordinates of the system (\f$ \textbf{r}\f$).
  !! \param lambda Constant to premultiply the perturbation (\f$ \lambda\f$).
  !! \param prt_bml Perturbation in bml format (\f$ H^{(1)} \f$).
  !! \param threshold Threshold value for bml format matrices.
  !! \param norbi Number of orbitals for each atomic site.
  !! \param verbose Different levels of verbosity.
  !! \param over_bml It has to be present for a nonorthogonal representation (\f$ S \f$).
  !!
  subroutine prg_pert_sin_pot(direction,lx,coordinate,lambda&
       ,prt_bml,threshold,spindex,norbi,verbose,over_bml)

    integer                                  ::  cont, i, ii, nats
    integer                                  ::  norb
    integer, intent(in)                      ::  spindex(:), norbi(:), verbose
    real(dp)                                 ::  fieldnorm, fieldx, fieldy, fieldz
    real(dp)                                 ::  lx
    real(dp)                                 ::  intensity, lambda, threshold
    real(dp), allocatable                    ::  diag(:)
    real(dp), intent(in)                     ::  coordinate(:,:)
    type(bml_matrix_t)                       ::  aux_bml
    type(bml_matrix_t), intent(inout)        ::  prt_bml
    type(bml_matrix_t), optional, intent(in) ::  over_bml
    character                                ::  direction

    write(*,*)'Applying a constant field perturbation ...'

    nats=size(spindex,dim=1)
    norb=sum(norbi(spindex(:)))

    write(*,*)nats,norb

    cont=0;

    allocate(diag(norb))

    do i=1,nats
       do ii=1,norbi(spindex(i))
          cont=cont+1
          diag(cont)=lambda*sin(2.0_dp*pi*((coordinate(1,i)-lx)/lx)-pi)
       enddo
    enddo

    if(bml_get_N(prt_bml) < 0) stop "bml_pert not allocated"

    call bml_set_diagonal(prt_bml,diag,threshold)

    deallocate(diag)

    if(present(over_bml))then  !(S*V + V*S)/2
       call bml_copy_new(prt_bml,aux_bml)
       call bml_multiply(over_bml,aux_bml,prt_bml,0.5_dp, 0.0_dp,threshold)
       call bml_multiply(aux_bml,over_bml,prt_bml,0.5_dp,1.0_dp,threshold)
       call bml_deallocate(aux_bml)
    endif

  end subroutine prg_pert_sin_pot

  !> Apply a cosine length dependent potential
  !! (\f$ \cos(\tilde{\textbf{r}}_x) \f$) where \f$ \textbf{r}_x \f$ is the x coordinate.
  !! The Hamiltonian gets modified as follows:
  !! \f$ H^{(1)} = \frac{1}{2}\lambda (S \sin(\tilde{\textbf{r}}_x) + \sin(\tilde{\textbf{r}}_x) S) \f$.
  !! \f$ \tilde{\textbf{r}}_x = 2\pi(\textbf{r}/l_x) - \pi \f$.
  !! The symmetrization is done in order to preserve the Hermiticity of H.
  !! Units can be transformed by using the parameter \f$ \lambda \f$.
  !! \note If the Hamiltonian is already in the prg_orthogonalized form, then parameter over_bml
  !! can be omitted.
  !! \param direction Direction of the potential gradient (x,y or z).
  !! \param lx Lenght of the box in x direction.
  !! \param coordinate Coordinates of the system (\f$ \textbf{r}\f$).
  !! \param lambda Constant to premultiply the perturbation (\f$ \lambda\f$).
  !! \param prt_bml Perturbation in bml format (\f$ H^{(1)} \f$).
  !! \param threshold Threshold value for bml format matrices.
  !! \param norbi Number of orbitals for each atomic site.
  !! \param verbose Different levels of verbosity.
  !! \param over_bml It has to be present for a nonorthogonal representation (\f$ S \f$).
  !!
  subroutine prg_pert_cos_pot(direction,lx,coordinate,lambda&
       ,prt_bml,threshold,spindex,norbi,verbose,over_bml)

    integer                                  ::  cont, i, ii, nats
    integer                                  ::  norb
    integer, intent(in)                      ::  spindex(:), norbi(:), verbose
    real(dp)                                 ::  fieldnorm, fieldx, fieldy, fieldz
    real(dp)                                 ::  lx
    real(dp)                                 ::  intensity, lambda, threshold
    real(dp), allocatable                    ::  diag(:)
    real(dp), intent(in)                     ::  coordinate(:,:)
    type(bml_matrix_t)                       ::  aux_bml
    type(bml_matrix_t), intent(inout)        ::  prt_bml
    type(bml_matrix_t), optional, intent(in) ::  over_bml
    character                                ::  direction

    write(*,*)'Applying a constant field perturbation ...'

    nats=size(spindex,dim=1)
    norb=sum(norbi(spindex(:)))

    write(*,*)nats,norb

    cont=0;

    allocate(diag(norb))

    do i=1,nats
       do ii=1,norbi(spindex(i))
          cont=cont+1
          diag(cont)=lambda*cos(2.0_dp*pi*((coordinate(1,i)-lx)/lx)-pi)
       enddo
    enddo

    if(bml_get_N(prt_bml) < 0) stop "bml_pert not allocated"

    call bml_set_diagonal(prt_bml,diag,threshold)

    deallocate(diag)

    if(present(over_bml))then  !(S*V + V*S)/2
       call bml_copy_new(prt_bml,aux_bml)
       call bml_multiply(over_bml,aux_bml,prt_bml,0.5_dp, 0.0_dp,threshold)
       call bml_multiply(aux_bml,over_bml,prt_bml,0.5_dp,1.0_dp,threshold)
       call bml_deallocate(aux_bml)
    endif

  end subroutine prg_pert_cos_pot

  !> Finds the first order response matrix from a Hamiltonian matrix
  ! and a perturbation V by purification. The method implemented here is the SP2 method
  ! [A. Niklasson, Phys. Rev. B, 66, 155115 (2002)].
  !! \param ham_bml Hamiltonian in bml format (\f$ H^{(0)} \f$).
  !! \param prt_bml Perturbation in bml format (\f$ H^{(1)} \f$).
  !! \param rsp_bml First order response to the perturbation (\f$ \rho^{(1)} \f$).
  !! \param bndfil Filing factor.
  !! \param threshold Threshold value for matrix elements.
  !! \param verbose Different levels of verbosity.
  !! \warning This works only for the prg_orthogonalized form of ham_bml.
  !! \warning The response must be in the prg_orthogonalized form.
  !!
  subroutine prg_compute_response_SP2(ham_bml,prt_bml,rsp_bml,rho_bml,lambda&
       ,bndfil,minsp2iter,maxsp2iter,sp2conv,idemtol,threshold,verbose)

    character(20)                      ::  bml_type
    integer                            ::  i, j, l, norb
    integer                            ::  verbose, mdim
    integer, intent(in)                ::  minsp2iter, maxsp2iter
    real(dp)                           ::  lambda, occ, maxminusmin
    real(dp)                           ::  alpha, beta, trx, trd
    real(dp), allocatable              ::  aux(:,:), evals(:), row(:), gbnd(:)
    real(dp), intent(in)               ::  bndfil, threshold,idemtol
    character(len=*), intent(in)       ::  sp2conv
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, occupation_bml
    type(bml_matrix_t)                 ::  x2_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml, prt_bml
    type(bml_matrix_t), intent(inout)  ::  rsp_bml, rho_bml
    real(dp)                           ::  idemperr, idemperr1, idemperr2
    integer                            ::  iter, breakloop
    real(dp)                           ::  trx2, trxOld, tr2xx2
    real(dp)                           ::  limdiff
    real(dp), allocatable              ::  trace(:)

    if(verbose.ge.1) write(*,*)'In prg_compute_response_SP2... '

    norb = bml_get_N(ham_bml)
    bml_type = bml_get_type(ham_bml)
    mdim = bml_get_M(ham_bml)

    idemperr = 0.0_dp
    idemperr1 = 0.0_dp
    idemperr2 = 0.0_dp

    if(bml_get_N(rsp_bml).le.0)then
       call bml_zero_matrix(bml_matrix_dense,bml_element_real,dp,norb,norb,rsp_bml)
    endif

    occ = bndfil*float(norb)

    call bml_copy(ham_bml, rho_bml)
    call bml_copy(prt_bml,rsp_bml)

    allocate(gbnd(2))
    call bml_gershgorin(ham_bml, gbnd)
    maxminusmin = gbnd(2) - gbnd(1)
    alpha = -1.00_dp/maxminusmin
    beta = gbnd(2)/maxminusmin
    call bml_scale_add_identity(rho_bml, alpha, beta, 0.00_dp)
    call bml_scale(alpha, rsp_bml) !Same as SP2 but with the perturbation.

    allocate(trace(2))
    trx = bml_trace(rho_bml)
    trd = bml_trace(rsp_bml)

    iter = 0
    breakloop = 0

    ! X2 <- X
    call bml_copy_new(rho_bml, x2_bml)
    call bml_copy_new(rsp_bml, aux_bml)

    do while (breakloop .eq. 0 .and. iter .lt. maxsp2iter)
       iter = iter + 1

       call bml_print_matrix("rsp_bml",rsp_bml,0,1,0,1)

       ! X2 <- X * X
       call bml_multiply_x2(rho_bml, x2_bml, threshold, trace)
       trd = bml_trace(rsp_bml)

       !Compute anticonmutator {X_n^0,Delta_n}
       call bml_multiply(rho_bml, rsp_bml, aux_bml, 1.0d0, 0.0d0)
       call bml_multiply(rsp_bml,rho_bml, aux_bml, 1.0d0, 1.0d0)

       trx2 = trace(2)
       write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2
       write(*,*) 'iter = ', iter, 'tr(resp) = ', trd

       tr2xx2 = 2.0_dp*trx - trx2
       trXOld = trx
       limDiff = abs(trx2 - occ) - abs(tr2xx2 - occ)

       if (limdiff .ge. idemtol) then

          ! X <- 2 * X - X2
          call bml_add_deprecated(2.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)
          call bml_add_deprecated(2.0_dp,rsp_bml,-1.0_dp, aux_bml, threshold)

          trx = 2.0_dp * trx - trx2

       elseif(limdiff .lt. -idemtol) then

          ! X <- X2
          call bml_copy(x2_bml, rho_bml)
          call bml_copy(aux_bml, rsp_bml)

          trx = trx2

       else

          trx = trxOld
          breakloop = 1

       end if

       idemperr2 = idemperr1
       idemperr1 = idemperr
       idemperr = abs(trx - trxOld)

       if (sp2conv .eq. "Rel" .and. iter .ge. minsp2iter .and. &
            (idemperr .ge. idemperr2 .or. idemperr .lt. idemtol)) then
          breakloop = 1
       end if

       if (iter .eq. maxsp2iter) then
          write(*,*) "SP2 purification is not converging: STOP!"
          stop
       end if

    enddo

    ! X <- 2 * X
    call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D
    call bml_scale(2.0_dp, rsp_bml)

    call bml_deallocate(x2_bml)
    call bml_deallocate(aux_bml)

    deallocate(trace)

  end subroutine prg_compute_response_SP2

  !> Project the response onto atomic positions.
  !! First order response to the perturbation (\f$ \rho^{(1)} \f$)
  !! projected onto the atomic position. Basically:
  !! \f$ rsp(i) = \sum_{\alpha \in i}\rho^{(1)}_{\alpha \alpha} \f$, where
  !! orbital \f$ \alpha \f$ belong to atom \f$ i \f$.
  !!
  !! \param rsp_bml First order response density matrix.
  !! \param spindex It gives the species index of a particular atom.
  !! \param norbi Number of orbitals of species i.
  !! \param coordinates Atomic coordinates.
  !! \param rspfunc Response function at atomic positions.
  !! \param verbose Different levels of verbosity.
  !!
  subroutine prg_project_response(rsp_bml,over_bml,spindex,norbi,coordinates,rspfunc,verbose)

    integer                              ::  cont, i, j, nats
    integer                              ::  norbs
    integer, intent(in)                  ::  norbi(:), spindex(:), verbose
    real(dp), allocatable                ::  diagonal(:)
    real(dp), allocatable, intent(inout)  ::  rspfunc(:)
    real(dp), intent(in)                 ::  coordinates(:,:)
    type(bml_matrix_t)                   ::  aux_bml
    type(bml_matrix_t), intent(in)       ::  over_bml
    type(bml_matrix_t), intent(inout)    ::  rsp_bml

    nats=size(spindex, dim=1)
    norbs=bml_get_N(rsp_bml)
    call bml_copy_new(rsp_bml,aux_bml)
    call bml_multiply(rsp_bml,over_bml,aux_bml,1.0_dp,0.0_dp,0.0_dp)

    if(.not. allocated(rspfunc))allocate(rspfunc(nats))
    allocate(diagonal(norbs))
    diagonal=0.0_dp
    call bml_get_diagonal(aux_bml,diagonal)
    ! call bml_get_diagonal(rsp_bml,diagonal)

    cont=0
    rspfunc = 0.0_dp
    do i=1,nats
       do j=1,norbi(spindex(i))
          cont=cont+1
          rspfunc(i) = rspfunc(i) + diagonal(cont)
       enddo
    enddo

    deallocate(diagonal)

    call bml_deallocate(aux_bml)

  end subroutine prg_project_response

end module prg_response_mod
