!> Computes density of state (DOS) for molecules or 2yy
!! user-defined groups of atoms within a molecular system
module gpmdcov_dos_mod
use bml
contains
        !> Computes density of state (DOS) for a moleculer system 
        !!
        subroutine compute_dos(num_points, sigma, emin, emax, evals, mu, filename)
                use prg_openfiles_mod
                
                implicit none
                integer, parameter :: dp = kind(1.0d0)
                integer, intent(in) :: num_points
                real(dp), intent(in) :: sigma, emin, emax, mu
                real(dp), allocatable, intent(in) :: evals(:)
                character(*), intent(in) :: filename

                integer :: fileunit, i, j, nevals
                real(dp), allocatable :: dos(:)
                real(dp) :: ecut, de, eps
               
                !write(*,*) "gpmdcov_dos_mod: Input variables ",num_points,sigma,&
                !        &emin,emax,filename
                !>
                !! Open output file
                call prg_open_file(fileunit, filename)

                ecut = 2.0_dp
                
                allocate(dos(num_points))
                nevals = size(evals)
                de = (emax - emin) / real(num_points)
                !>
                !! Compute all contributions from all eigenvalues (evals)
                !! For every point, compute the Gaussian contributions
                !! of each eval
                dos = 0.0_dp
                do i = 1, num_points
                        eps = emin + de * (i - 1)
                        do j = 1, nevals
                                dos(i) = dos(i) + gaussian(eps, evals(j), sigma)
                        enddo 
                        write(fileunit,*) eps - mu, dos(i)
                enddo

                deallocate(dos)
                close(fileunit)

        end subroutine compute_dos

        !> Computes local density of state (DOS) for a user-defined group of atoms
        !! within a  moleculer system 
        !!
        subroutine compute_local_dos(num_points, atoms, hindex, sigma, emin, &
                        &emax, evects, evals, overlaps, mu, zmat_bml, symbols, lfilename)
                use prg_openfiles_mod

                implicit none
                integer, parameter :: dp = kind(1.0d0)
                integer, intent(in) :: num_points
                integer, allocatable :: atoms(:)
                integer, allocatable, intent(in) :: hindex(:,:)
                real(dp), intent(in) :: sigma, emin, emax, mu
                real(dp), allocatable, intent(in) :: evals(:)
                Type(bml_matrix_t), intent(in) :: evects, overlaps, zmat_bml
                character(*), intent(in) :: symbols(:)
                character(*), intent(in) :: lfilename

                integer :: fileunit, a, i, j, k, m, nevals, norb, nstates, at
                real(dp), allocatable :: dos(:), row(:), weights_dense(:,:), weights_2(:,:)
                real(dp), allocatable :: subsystem_weight(:), atom_weights(:)
                real(dp), allocatable :: evects_dense(:,:), overlaps_dense(:,:)
                !real(dp) :: ecut, de, eps, sum_squares, subsystem_weight
                real(dp) :: ecut, de, eps, sum_squares
                Type(bml_matrix_t) :: weights_matrix,z_transpose, zc_matrix
                
                !! Allocate row matrix and dense versions of eivenvectors(C), overlaps (S), and weigths matricies
                allocate(row(size(evals)))
                allocate(subsystem_weight(size(evals)))
                nstates = bml_get_N(evects)
                allocate(weights_dense(nstates,nstates))
                allocate(weights_2(nstates,nstates))
                allocate(evects_dense(nstates,nstates))
                allocate(overlaps_dense(nstates,nstates))
               
                !! Create new bml matricies for SC product and  weights based on the number of orbitals 
                call bml_zero_matrix("dense",bml_element_real,dp,nstates,nstates,weights_matrix)
                call bml_zero_matrix("dense",bml_element_real,dp,nstates,nstates,zc_matrix)

                call bml_print_matrix("evects",evects,0,6,0,6)
                call bml_print_matrix("overlap",overlaps,0,6,0,6)
                !call bml_get_row(evects,1,row)
                !write(*,*) "Row of evects ",row
                
                !! Transpose Z matrix (Z)
                call bml_transpose_new(zmat_bml, z_transpose)
                !! Multiply Z matrix (Z) and eigenvectors (C) to form zc_matrix
                call bml_multiply(z_transpose,evects,zc_matrix,1.0_dp,0.0_dp,0.0_dp)
                
                !! Multiply overlaps matrix (S) and ZC matrix to form weights_matrix
                call bml_multiply(zc_matrix,overlaps,weights_matrix,1.0_dp,0.0_dp,0.0_dp)
                
                !! Print new weights matrix to check
                !call bml_print_matrix("weights",weights_matrix,0,6,0,6)
                
                !! Create dense matrices from BML matricies
                call bml_export_to_dense(overlaps,overlaps_dense) 
                call bml_export_to_dense(evects,evects_dense) 
                call bml_export_to_dense(weights_matrix,weights_dense) 
                weights_2 = matmul(evects_dense, evects_dense)
                
                sum_squares = 0.0_dp
                do i = 1,nstates
                    sum_squares = sum_squares + weights_dense(2,i)**2
                    !sum_squares = sum_squares + evects_dense(2,i)**2
                enddo
                write(*,*) "Sum Squares ",sum_squares 
                
                !! Open output file
                call prg_open_file(fileunit, lfilename)

                !! Compute values of w_i
                subsystem_weight=0.0_dp

                !do i = 1,size(atoms)
                do k = 1,nstates
                  do a = 1,size(atoms)
                  !if (symbols(atoms(i)) .eq. "C") then
                    do j = hindex(1,atoms(a)),hindex(2,atoms(a))
                      subsystem_weight(k) = subsystem_weight(k) + weights_dense(k,j)**2
                      !subsystem_weight(k) = subsystem_weight(k) + evects_dense(k,j)**2
                      write(*,*) "Add to subsystem weight ",weights_dense(k,j)**2
                      !write(*,*) "Add to subsystem weight ",evects_dense(k,j)**2
                    enddo
                  !endif
                  enddo
                enddo

                write(*,*) "Subsystem weight ",subsystem_weight

                ecut = 2.0_dp
                
                allocate(dos(num_points))
                nevals = size(evals)
                write(*,*) "Nstates and nevals ",nstates,nevals 
                de = (emax - emin) / real(num_points,dp)
                !>
                !! Compute all contributions from all eigenvalues (evals)
                !! For every point, compute the Gaussian contributions
                !! of each eval
                dos = 0.0_dp
                do m = 1, num_points
                        eps = emin + de * (m - 1)
                        do j = 1, nevals
                               ! dos(i) = dos(i) + gaussian(eps, evals(j), sigma)
                                dos(m) = dos(m) + gaussian(eps, evals(j), sigma) * subsystem_weight(j)
                                
                        enddo
                        write(fileunit,*) eps - mu, dos(m)
                enddo

                !>
                ! Print out weights of different atoms
                ! Want to see what all atoms have the most weight right around Fermi level
                !allocate(atom_weights(sy%nats))
                !atom_weights = 0.0_dp

                !do k = 1, nstates
                !  do at = 1,sy%nats
                !    do i = hindex(at,1), hindex(at,2)
                !      atom_weights(at) = atom_weights(at) + weights_dense(k,j)
                !    enddo
                !  enddo
                !enddo

                !do at = 1,sy%nats
                !  write(*,*) "Atom weight ",at,atom_weights(at)
                !enddo

                !deallocate(atom_weights)
                deallocate(dos)
                close(fileunit)

        end subroutine compute_local_dos

        !> Function to compute the Gaussian contibution for an eigenvalue
        !!
        function gaussian(eps, eps_i, sigma) result(G)
                implicit none
                integer, parameter :: dp = kind(1.0d0)
                real(dp), intent(in) :: eps, eps_i, sigma
                real(dp) :: pi, delta, arg, myexp
                real(dp) :: G
                delta = 10E-5
                pi = 3.1415926535897932384626_dp
                arg = ((eps - eps_i)/sigma)**2
                if (arg .lt. 1E-5) then
                        myexp = 1.0_dp
                        write(*,*) "myexp 1.0", arg
                else
                        myexp = exp(-0.5*arg)
                        write(*,*) "myexp ",myexp,arg
                endif
                G = 1.0 / (sigma * sqrt(2.0 * pi)) * myexp 
        end function gaussian
        
end module gpmdcov_dos_mod
