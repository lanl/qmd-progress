
!> \defgroup LATTE (LATTE related routines)
!> \defgroup PROGRESS (PROGRESS related routines)
!> \defgroup EXTERNAL (EXTERNAL related routines)
!> \defgroup PROGRAMS (High-level codes using PROGRESS/LATTE modules)


!> Glossary for the code 
!! 
module variable_names

  !> Vector containing a set of eigenvalues. 
  read(dp),  allocatable  ::  eigenvalues(:)                   

  !> Number of orbitals. Its usually the dimension of the Hamiltonian.
  integer                 ::  norb                             

  !> Number of occupied orbitals. 
  integer                 ::  nocc                   

  !> Number of atoms.
  integer                 ::  nats       

  !> Threshold value for sparse matrices.
  real(dp)                ::  threshold         

  !> Coefficient matrix in bml format.
  type(bml_matrix_t)      ::    eigenvectors_bml         

  !> Character containing the bml format type.
  character(20)           ::    bml_type                  

  !> System type containing coordinates, atom types, velocities, etc.     
  type(system_type)       ::    system                            

  !> Hamiltonian in bml format.    
  type(bml_matrix_t)      ::   ham_bml                 

  !> Density matrix in bml format.         
  type(bml_matrix_t)      ::  rho_bml                

  !> Inverse overlap factor in bml format.         
  type(bml_matrix_t)      ::  zmat_bml                

  !> Diagonal matrix containing the occupations.          
  type(bml_matrix_t)      ::   occupation_bml         

  !> Auxiliary bml matrix.         
  type(bml_matrix_t)      ::  aux_bml                

  !> A general bml matrix.          
  type(bml_matrix_t)      ::  mat_bml                

  !> Integer caring different levels of verbosity.          
  integer                 ::     verbose                

  !> File extension ("xyz" or "pdb").       
  character(3)            ::    extension 

  !> File extension ("xyz" or "pdb").       
  character(3)            ::  ext

  !> Variable to contain a file name.           
  character(20)           ::    filename      

  !> Input unit Number to write or read.        
  integer                 ::       io_unit                

  !Dummy variable to read.     
  character(20)           ::  dummy                  

end module variable_names

