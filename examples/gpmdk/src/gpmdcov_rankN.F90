r = K0*(q[n]_core-n_core)
do i = 1, Nranks

        do l = 1, nparts

        ! Compute a part for the preconditioner residual vector
        ! K0Res_core = K0_core*r

        enddo

        ! v_i = r/norm(r)
        ! v_i orthonormalized to all v_j, j = 1,2,..., i-1
        ! v_i (charge perturbation vector) -> VCoulomb

        do l = 1, nparts

        ! VCoulomb for core + halo of each subgraph
        ! DMPRT(H0,H1 = VCoulomb, ....) -> P1' (from part), 
        ! dPmu = beta*P0*(I-P0) (from part), 
        ! Core traces of P1', and dPmu, i.e. TrP1' and TrdPmu (fore each
        ! core)

        enddo

        ! mu_1 = -Sum TrP1' / Sum TrPmu over all graphs
         
        do l = 1, nparts

        ! P1 = P1' + mu_1*dPmu  (core + halo)
        ! Trace(P1*S) for each atom and core -> dq[n+lambda*v_i]/dlambda

        enddo

        ! r = dq[n+lambda*v_i]/d_lambda - v_i
        ! r = K0*r
        ! f_vi = r                
        
enddo        
