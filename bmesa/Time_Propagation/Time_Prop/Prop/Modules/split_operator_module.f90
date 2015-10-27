!***********************************************************************
! split_operator_module
!**begin prologue     split_operator_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to propagate
!***                  a wavefunction in time using the SO
!***                  algorithm.  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       Given a starting vector, the wavefunction is propagated 
!***                  until a certain convergence criterion is reached.
!***                  In imaginary time, the wavefunction is propagated 
!***                  until the eigenvalue has converged to a given accuracy. 
!***                  For real time, one just proceeds from step to step. 
!**references
!**modules needed     See USE statements below
!**end prologue       arnoldi_module
!***********************************************************************
                     MODULE split_operator_module
                     USE arnoldi_global
                     USE dvrprop_global
                     USE dvr_global
                     USE dvr_shared
                     USE regional_module
                     USE normalize_module
                     USE moment_module
                     USE auto_correlation_module
                     USE plot_module
                     USE spatial_wavefunction_module
                     USE initial_state_module
                     USE non_linear_potential_module
                     USE matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_sector_time_prop
             MODULE PROCEDURE so_sector_time_prop_d,                    &
                              so_sector_time_prop_z
                 END INTERFACE so_sector_time_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_propagation
             MODULE PROCEDURE so_propagation_d,                         &
                              so_propagation_z
                 END INTERFACE so_propagation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_exp
            MODULE PROCEDURE so_exp_d,                                  &
                             so_exp_z
                 END INTERFACE so_exp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_exp_diagonal_propagator
            MODULE PROCEDURE so_exp_diag_prop_d,                        &
                             so_exp_diag_prop_z
                 END INTERFACE so_exp_diagonal_propagator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_exp_diagonal_multiply
            MODULE PROCEDURE so_exp_diagonal_mul_d,                     &
                             so_exp_diagonal_mul_z
                 END INTERFACE so_exp_diagonal_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_exp_off_diagonal_multiply
            MODULE PROCEDURE so_exp_off_diagonal_mul_d,         &
                             so_exp_off_diagonal_mul_z      
                 END INTERFACE so_exp_off_diagonal_multiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck SO_Sector_Time_Prop_d
!***begin prologue     Split_Operator_Sector_Time_Prop_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           split-operator, propagation, RSP2, RSP4
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate the sector time propagators
!***                   based on Lie-Trotter-Suziki formulation.
!***                   A second and fourth order method have
!***                   been programmed.  the splitting is based
!***                   on the finite element DVR approach.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       SO_Sector_Time_Prop_d
  Subroutine SO_Sector_Time_Prop_d(d_exp)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                   :: d_exp
  REAL*4                                   :: secnds
  INTEGER                                  :: i, j
!
!        Calculate the  time dependent propagators
!
  DO i=1,spdim
     WRITE(iout,1) i
     p_fac = 1.d0
     n_prop =1
     IF(prop_order == 4) THEN
        p_fac = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
        n_prop = 2
     END IF
     DO j=1,num_reg(i)
        ALLOCATE(mat_reg(j,i)%exp_t_mat_d                        &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop))         
     END DO
     starting_reg = 1
     ending_reg = num_reg(i)
     n_reg = num_reg(i)
     write(iout,2)
     call regional_prop_d(i,deltat)
     DO j=1,num_reg(i)
        DEALLOCATE(mat_reg(j,i)%eigval_mat_d,                  &
                   mat_reg(j,i)%eigvec_mat_d)
     END DO
  END DO
1 FORMAT(/,15x,'Sector Time Propagtors Variable = ',i2)
2 FORMAT('Constructing Propagators for Sectors with Real '     &
            'Potentials')
3 FORMAT('Constructing Propagators for Sectors with Complex '  &
            'Potentials')
END SUBROUTINE SO_Sector_Time_Prop_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck SO_Sector_Time_Prop_z
!***begin prologue     SO_Sector_Time_Prop_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           split-operator, propagation, RSP2, RSP4
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate the sector time propagators
!***                   based on Lie-Trotter-Suziki formulation.
!***                   A second and fourth order method have
!***                   been programmed.  the splitting is based
!***                   on the finite element DVR approach.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       SO_Sector_Time_Prop_z
  Subroutine SO_Sector_Time_Prop_z(z_exp)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)               :: z_exp
  REAL*4                                   :: secnds
  INTEGER                                  :: i, j
!
!        Calculate the  time dependent propagators
!
  DO i=1,spdim
     WRITE(iout,1) i
     p_fac = 1.d0
     n_prop =1
     IF(prop_order == 4) THEN
        p_fac = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
        n_prop = 2
     END IF
     DO j=1,num_reg(i)
        ALLOCATE(mat_reg(j,i)%exp_t_mat_z                        &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop))         
     END DO
     IF(.not.absorb) THEN
        starting_reg = 1
        ending_reg = num_reg(i)
        n_reg = num_reg(i)
        WRITE(iout,2)
        call regional_prop_h(i,deltat)
     ELSE
        starting_reg = 1
        ending_reg = n_reg_real(i)
        n_reg = n_reg_real(i)
        write(iout,2)
        CALL regional_prop_h(i,deltat)
        starting_reg = n_reg_real(i) + 1
        ending_reg = num_reg(i)
        n_reg = ending_reg - starting_reg + 1
        write(iout,3)
        call regional_prop_z(i,.5d0*deltat)
     END IF
     DO j=1,n_reg_real(i)
        DEALLOCATE(mat_reg(j,i)%eigval_mat_d,                  &
                   mat_reg(j,i)%eigvec_mat_d)
     END DO
     IF(absorb) THEN
        DO j=n_reg_real(i)+1, num_reg(i)
           DEALLOCATE(mat_reg(j,i)%v_add_z,                    &
                      mat_reg(j,i)%eigval_mat_z,               &
                      mat_reg(j,i)%eigvec_mat_z_r,             &
                      mat_reg(j,i)%eigvec_mat_z_l)
        END DO
     END IF
  END DO
1 FORMAT(/,15x,'Sector Time Propagtors Variable = ',i2)
2 FORMAT('Constructing Propagators for Sectors with Real '     &
            'Potentials')
3 FORMAT('Constructing Propagators for Sectors with Complex '  &
            'Potentials')
END SUBROUTINE SO_Sector_Time_Prop_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_propagation_d
!**begin prologue     so_propagation_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, real-space, split-operator, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             FEDVR_Prop
!**purpose            time dependent schroedinger equation
!**                   using split-operator method with 
!**                   finite difference or dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       so_propagation_d
  SUBROUTINE so_propagation_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)              :: wave_function
  REAL*8, DIMENSION(n3d)              :: scratch_vector
  CHARACTER (LEN=2)                   :: itoc
  INTEGER                             :: i, t, tst
  INTEGER                             :: length, len
  REAL*8                              :: norm
  COMPLEX*16                          :: auto
  REAL*8                              :: total_t
!
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     IF (e_method == 'hamiltonian') THEN
!
!           Calculate Energy in propagation using Hamiltonian
!
         ALLOCATE(buf(spdim))
         DO i=1,spdim
            ALLOCATE( buf(i)%d(nphy(i)),              &
                      buf(i)%hbuf(nphy(i)*nphy(i)),   &
                      buf(i)%hibuf(2,nphy(i)*nphy(i)))
            CALL pack_h(i)
         END DO
     END IF
  END IF
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  eig_old=0.0d0
  CALL initial_vector(wave_function)
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
  
!        Calculate the time dependent perturbation.  It consists of a
!        space and a time part.
!
      v_tot = 0.d0
      CALL v_tim
!
!        In subroutine pert, the one-body potentials are added to
!        any more complicated interaction.  Note that, the one-body
!        potentials may contain the diagonal contributions from the one-body
!        kinetic energy.
!
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$. 
!
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!     the main real space product driving routine
!
      CALL so_exp(wave_function,scratch_vector)  
!
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( tst == 0 ) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function,t)
     END IF
     IF (e_method == 'exponential') then
         CALL gs_energy(scratch_vector,wave_function,e_cur)
     ELSE
         call check_norm(wave_function,norm)
         wave_function = wave_function/sqrt(norm) 
!         call h_v_d(wave_function,scratch_vector,1)
         call finite_element_h_v_d(wave_function,scratch_vector,1)
         call check_gs_energy(wave_function,scratch_vector,e_cur)
     END IF
     eig_tst=abs(e_cur-eig_old)
     write(iout,3) t, tim_pts(t+1), norm, e_cur, eig_tst
     IF( eig_tst <= con ) THEN
!
!        End the propagation
!
         CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
         WRITE(iout,1)
         WRITE(iout,4)
         WRITE(iout,1)
         RETURN
     ELSE
         eig_old=e_cur
     END IF  
     call iosys('write real solution to bec',n3d,                     &
                 wave_function,0,' ')
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function,norm)
  wave_function = wave_function/sqrt(norm)
  CALL plot_psi(wave_function,t-1) 
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     IF (e_method == 'hamiltonian') THEN
         DO i=1,spdim
            DEALLOCATE( buf(i)%d,buf(i)%hbuf,buf(i)%hibuf )
         END DO
         DEALLOCATE(buf)
     END IF
  END IF
!
1 FORMAT('***********************************************'            &
         '*************************')
2 FORMAT(/,20x,'Begin Split Operator Propagation')
3 FORMAT(/,5x,'Time Step     = ',i6,5x,'Time = ',f15.8,               &
         /,5x,'Normalization = ',e15.8,1x,'Energy    = ',e15.8,       &
         /5x,'Convergence = ',e15.8)
4 FORMAT(/,20x,'End Split Operator Propagation')
!
END SUBROUTINE so_propagation_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_propagation_z
!**begin prologue     so_propagation_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, real-space, split-operator, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             FEDVR_Prop
!**purpose            time dependent schroedinger equation
!**                   using split-operator method with 
!**                   finite difference or dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       so_propagation_z
  SUBROUTINE so_propagation_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)          :: wave_function
  COMPLEX*16, DIMENSION(n3d)          :: scratch_vector
  CHARACTER (LEN=2)                   :: itoc
  INTEGER                             :: i, t, tst
  INTEGER                             :: length, len, count
  REAL*8                              :: norm, rtemp
  COMPLEX*16                          :: auto
  REAL*8                              :: total_t
!
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  count=0
  CALL initial_vector(wave_function)
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
  
!        Calculate the time dependent perturbation.  It consists of a
!        space and a time part.
      v_tot = 0.d0
      CALL v_tim
!
!        In subroutine pert, the one-body potentials are added to
!        any more complicated interaction.  Note that, the one-body
!        potentials contain the diagonal contributions from the one-body
!        kinetic energy.
!
      CALL pert
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.  The first time through
!     the initial state is stored in soln which is needed later.
!
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL so_exp(wave_function,scratch_vector)   
      call iosys('write real solution to bec',2*n3d,wave_function,0,' ')
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( count == plot_step ) THEN
        count=0
        IF(i0stat == 'gaussian-pulse') THEN
           call calculate_moment(wave_function,tim_pts(t+1))
        END IF
!
!       Form the total solution and compute the autocorrelation
!       function.

        CALL spatial_psi(wave_function,scratch_vector)
        CALL plot_psi(scratch_vector,t)
        CALL auto_correlation_function(auto,scratch_vector,     &
                                       wave_function)
        rtemp = 1.d0 - auto * conjg(auto)
        call check_norm(wave_function,norm)
        write(iout,3) t, tim_pts(t+1), norm, rtemp
     END IF
     count=count+1
  END DO
  WRITE(iout,1)
  WRITE(iout,4)
  WRITE(iout,1)
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
!
!
1 FORMAT('***********************************************'              &
         '*************************')
2 FORMAT(/,20x,'Begin Split Operator Propagation')
3 FORMAT(/,5x,'Time Step           = ',i6,5x,'Time = ',f15.8,           &
         /,5x,'Normalization       = ', e15.8,5x,'Survival Probablity = ', e15.8)
4 FORMAT(/,20x,'End Split Operator Propagation')
END SUBROUTINE so_propagation_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_d.f
!***begin prologue     so_exp_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       so_exp_d
  SUBROUTINE so_exp_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)       :: wave_function
  REAL*8, DIMENSION(n3d)       :: scratch_vector
!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     CALL print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!
!            First Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_d)
      END IF
!
      CALL so_exp_diagonal_multiply(wave_function)
!
     IF(log_main(9)) then
         title='First Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Off diagonal scaling
!                t = p_fac * deltat / hbar
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
      IF(log_main(9)) then
         title='Off Diagonally Propagated Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Second Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      CALL so_exp_diagonal_multiply(wave_function)
!
      IF(log_main(9)) then
         title='Second Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
  ELSE IF (prop_order == 4) THEN
!
!            First Diagonal Scaling by Potential
!               t = p_fac * deltat * .5 / h_bar  
!
      tau_loc = p_fac * deltat * .5d0 / hbar
!
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!               t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Second Diagonal Scaling by potential
!               t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Third Diagonal Scaling
!               t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Fourth Diagonal Scaling
!               Same as the third.
!
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!     OK.  We are done.
! 
  END IF
END SUBROUTINE so_exp_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_z.f
!***begin prologue     so_exp_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       so_exp_3d_z
  SUBROUTINE so_exp_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)       :: wave_function
  COMPLEX*16, DIMENSION(n3d)       :: scratch_vector
!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     CALL print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!
!            First Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z) 
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
!
      CALL so_exp_diagonal_multiply(wave_function)
!
     IF(log_main(9)) then
        title='First Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Off diagonal scaling
!                t = p_fac * deltat / hbar
     prop_point = 1
     CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='Off Diagonally Propagated Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Second Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
     CALL so_exp_diagonal_multiply(wave_function)
!
     IF(log_main(9)) then
        title='Second Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
      write(iout,*) prop_order
  ELSE IF (prop_order == 4) THEN
!
!            First Diagonal Scaling by Potential
!               t = p_fac * deltat * .5 / h_bar  
!
      tau_loc = p_fac * deltat * .5d0 / hbar
!
      CALL so_exp_diagonal_propagator(exp_diag_z)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
      CALL so_exp_diagonal_multiply(wave_function)
      IF(log_main(9)) then
         title='First Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!               t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
      IF(log_main(9)) then
         title='First Off Diagonally Propagated Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Second Diagonal Scaling by potential
!               t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
      CALL so_exp_diagonal_multiply(wave_function)
      IF(log_main(9)) then
         title='Second  Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
      IF(log_main(9)) then
         title='Second Off Diagonally Propagated Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Third Diagonal Scaling
!               t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
      CALL so_exp_diagonal_multiply(wave_function)
      IF(log_main(9)) then
         title='Third  Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
      IF(log_main(9)) then
         title='Third Off Diagonally Propagated Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Fourth Diagonal Scaling
!               Same as the third.
!
      CALL so_exp_diagonal_multiply(wave_function)
      IF(log_main(9)) then
         title='Third  Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
      IF(log_main(9)) then
         title='Fourth Off  Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
      CALL so_exp_diagonal_multiply(wave_function)
      IF(log_main(9)) then
         title='Fifth Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
      IF(log_main(9)) then
         title='Fifth Off  Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
      CALL so_exp_diagonal_multiply(wave_function)
      IF(log_main(9)) then
         title='Fifth Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!     OK.  We are done.
! 
  END IF
END SUBROUTINE so_exp_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!deck so_exp_diag_prop_d
!***begin prologue     so_exp_diag_prop_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate exponential time-dependent scaling
!***                   for a local potential, v_tot.
!***description        In the split operator approach, local potentials,
!***                   which may be time dependent are split off before
!***                   treating the kinetic energy operators.  These only
!***                   scale the vectors and do not require the more complex
!***                   even/odd splitting of the kinetic energy.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       so_exp_diag_prop_d
  SUBROUTINE so_exp_diag_prop_d(exp_diag)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                  :: exp_diag
  exp_diag =   exp(- v_tot * tau_loc)  
END SUBROUTINE so_exp_diag_prop_d
!***********************************************************************
!***********************************************************************
!deck so_exp_diag_prop_z
!***begin prologue     so_exp_diag_prop_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate exponential time-dependent scaling
!***                   for a local potential, v_tot.
!***description        In the split operator approach, local potentials,
!***                   which may be time dependent are split off before
!***                   treating the kinetic energy operators.  These only
!***                   scale the vectors and do not require the more complex
!***                   even/odd splitting of the kinetic energy.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       so_exp_diag_prop_z
  SUBROUTINE so_exp_diag_prop_z(exp_diag)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)                  :: exp_diag
  exp_diag =   exp( - eye * v_tot * tau_loc)  
END SUBROUTINE so_exp_diag_prop_z
!***********************************************************************
!***********************************************************************
!deck so_exp_diagonal_mul_d
!***begin prologue     so_exp_diagonal_mul_d
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            diagonal scaling of propagated vector by exponential.
!***
!***description        
!                      
!                      
!***references
!***routines called    
!***end prologue       so_exp_diagonal_mul_d
!
  SUBROUTINE so_exp_diagonal_mul_d(vector_in_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                         :: vector_in_out
  vector_in_out(:) = exp_diag_d(:) * vector_in_out(:)  
END SUBROUTINE so_exp_diagonal_mul_d
!***********************************************************************
!***********************************************************************
!deck so_exp_diagonal_mul_z
!***begin prologue     exp_diagonal_mul_z
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            diagonal scaling of propagated vector
!***
!***description        
!                      
!                      
!***references
!***routines called    
!***end prologue       so_exp_diagonal_mul_z
!
  SUBROUTINE so_exp_diagonal_mul_z(vector_in_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)                       :: vector_in_out
  vector_in_out(:) = exp_diag_z(:) * vector_in_out(:)  
END SUBROUTINE so_exp_diagonal_mul_z
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_off_diagonal_mul
!***begin prologue     so_exp_off_diagonal_mul     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            compute the effect of the split operator
!                      off-diagonal, exponential propagator on a vector.
!
!***description        These routine perform the complex operation of
!***                   applying the exponential of a FEDVR Hamiltonian
!***                   on a vector.  The operation consists of applying,
!
!     exp(-iH_odd*t/(2*hbar)) * exp(-iH_even*t/(hbar)) * exp(-iH_odd*t/(2*hbar))
!
!                      to the vector.  H_even and H_odd are separately block
!                      diagonal but the blocks overlap.  Thus the operations must
!                      alternate their input vectors.
!
!                     H = H(1,1) + H(2,2) + H(3,3)
!
!     Note that vectors are stored as V(nphy(3),nphy(2),nphy(1)).  Thus matrix
!     vector multiply in two and three dimensions must be done carefully for
!     nphy(2) and nphy(3).  We have( repeated indices summed over), in 3D,
!    
! V(3,2,1) = [ H(1,1) * V(3,2,1) + H(2,2) * V(3,2,1) + H(3,3) * V(3,2,1) ]  
!                                   =
!            [ V(3,2,1) * H(1,1) + V(3,2,1) * H(2,2) + H(3,3) * V(3,2,1) ]
!
! Thus the first mutiply may be done as a matrix V(3*2,1) on H(1,1), 
! the second as V(3,2,1) * H(2,2) with an outer loop over index 1 and 
! the third as a simple matrix vector multiply, H(3,3) * V(3,2*1).
!
! In 2D, we have,
!
! V(2,1)   = [ H(1,1) * V(2,1) + H(2,2) * V(2,1)  ]  
!                                   =
!            [ V(2,1) * H(1,1) + H(2,2) * V(2,1) ]
!
!
!
!***references
!***routines called    
!***                  
!
!***end prologue       so_exp_off_diagonal_mul
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_d(wave_function,scratch_vector)
     IMPLICIT NONE
     REAL*8, DIMENSION(n3d)                   :: wave_function
     REAL*8, DIMENSION(n3d)                   :: scratch_vector
     IF (spdim == 1) THEN
         CALL exp_off_diag_m_v_d(wave_function,                      &
                                 scratch_vector,                     &
                                 nphy(1),1,1)
     ELSE IF ( spdim == 2) THEN
         CALL so_exp_off_diagonal_mul_2d_d(wave_function,            &
                                           scratch_vector)
     ELSE IF( spdim == 3) THEN
         CALL so_exp_off_diagonal_mul_3d_d(wave_function,            &
                                           scratch_vector)
     END IF
  END SUBROUTINE so_exp_off_diagonal_mul_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(n3d)                   :: wave_function
     COMPLEX*16, DIMENSION(n3d)                   :: scratch_vector
     IF (spdim == 1) THEN
         CALL exp_off_diag_m_v_z(wave_function,                      &
                                 scratch_vector,                     &
                                 nphy(1),1,1)
     ELSE IF ( spdim == 2) THEN
         CALL so_exp_off_diagonal_mul_2d_z(wave_function,            &
                                           scratch_vector)
     ELSE IF( spdim == 3) THEN
         CALL so_exp_off_diagonal_mul_3d_z(wave_function,            &
                                           scratch_vector)
     END IF
  END SUBROUTINE so_exp_off_diagonal_mul_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_2d_d(wave_function,scratch_vector)
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(2),nphy(1))           :: wave_function
     REAL*8, DIMENSION(nphy(2),nphy(1))           :: scratch_vector
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(2),
!            that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v_d(wave_function,                     &
                              scratch_vector,                    &
                              nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(1),
!            V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_v_m_2_d(wave_function,                   &
                                scratch_vector,                  &
                                nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_2d_d
  SUBROUTINE so_exp_off_diagonal_mul_3d_d(wave_function,scratch_vector)
     USE dvrprop_global
     USE dvr_shared
     USE dvr_global
     IMPLICIT NONE
     REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))   :: wave_function
     REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))   :: scratch_vector

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v_d(wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_3_d                                  &
                             (wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2),nphy(1),2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_2_d                                  &
                             (wave_function,                     &  
                              scratch_vector,                    &
                              nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck exp_off_diag_m_v_d
!***begin prologue     exp_off_diag_m_v_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = exp_t_mat * V
!
!------------------------------------------------------------------------------------
!***                   where exp_t_mat is the regional
!***                   propagator for a FEDVR or FD Hamiltonian.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D_D nj=ny*nx.

!***references
!***routines called    v_m_v_gen, v_m_v_2, v_m_v_3
!
!***end prologue       exp_off_diag_m_v_d                     
!
  SUBROUTINE exp_off_diag_m_v_d(v,                                             &
                                v_scr,                                         &
                                ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i, ntrips, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!  The sequence is called three times.  First for the odd matrices, then
!  for the even matrices and then again for the odd matrices.
! 
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips=3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
  DO trips=1,ntrips
     DO i = begin(trips), num_reg(index), 2
!

        IF ( nfun_reg(i,index) > 10 ) then
!
!                       General Code
!
           CALL v_so_mat_v_gen(v(locate(trips):,:),                               &
                               v_scr(locate(trips):,:),                           &
                               mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                               ni,nj,nfun_reg(i,index))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           CALL v_so_mat_v_2(v(locate(trips):,:),                                &
                             v_scr(locate(trips):,:),                            &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),       &
                             ni,nj)
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           CALL v_so_mat_v_3(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           CALL v_so_mat_v_4(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           CALL v_so_mat_v_5(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
           CALL v_so_mat_v_6(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           CALL v_so_mat_v_7(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           CALL v_so_mat_v_8(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           CALL v_so_mat_v_9(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           CALL v_so_mat_v_10(v(locate(trips):,:),                              &
                              v_scr(locate(trips):,:),                          &
                              mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                              ni,nj)

!        ELSE IF (nfun_reg(i,index) == 11) then
!
!                       Special case for Eleven by Eleven
!
!           CALL v_so_mat_v_11(v(locate(trips):,:),                              &
!                              v_scr(locate(trips):,:),                          &
!                              mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
!                              ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                       &
                                          +                                         &
                            nfun_reg(i+1,index)                                     &
                                          -                                         &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_m_v_d
!
!*deck exp_off_diag_v_m_2_d
!***begin prologue     exp_off_diag_v_m_2_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = V * exp_t_mat 
!
!------------------------------------------------------------------------------------
!***                   where exp_t_mat is the regional
!***                   propagator for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D_D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    v_v_m_gen, v_v_m_2, v_v_m_3
!***end prologue       exp_off_diag_v_m_2_d
!
  SUBROUTINE exp_off_diag_v_m_2_d(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index, ntrips
  REAL*8, DIMENSION(ni,nj)                 :: v
  REAL*8, DIMENSION(ni,nj)                 :: v_scr
  INTEGER                                  :: i, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips = 3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
!
  DO trips = 1, ntrips
     DO i=begin(trips), num_reg(index), 2
        IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
           CALL v_v_so_mat_gen(v(:,locate(trips):),                            &
                               v_scr(:,locate(trips):),                        &
                               mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),   &
                               ni,nj,nfun_reg(i,index)) 
!
        ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
           CALL v_v_so_mat_2(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                             ni,nj) 
       ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
           CALL v_v_so_mat_3(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
           CALL v_v_so_mat_4(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
           CALL v_v_so_mat_5(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
           CALL v_v_so_mat_6(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
           CALL v_v_so_mat_7(v(:,locate(trips):),                               &
                             v_scr(:,locate(trips):),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
           CALL v_v_so_mat_8(v(:,locate(trips):),                               &
                             v_scr(:,locate(trips):),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
           CALL v_v_so_mat_9(v(:,locate(trips):),                               &
                             v_scr(:,locate(trips):),                           &
                             mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),      &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
           CALL v_v_so_mat_10(v(:,locate(trips):),                              &
                              v_scr(:,locate(trips):),                          &
                              mat_reg(i,index)%exp_t_mat_d(:,:,prop_point),     &
                              ni,nj) 
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                     &
                                          +                                       &
                            nfun_reg(i+1,index)                                   &
                                          -                                       &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_2_d
!*deck exp_off_diag_v_m_3_d
!***begin prologue     exp_off_diag_v_m_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   exp_off_diag_v_m_2 in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       exp_off_diag_v_m_3_d
!
  SUBROUTINE exp_off_diag_v_m_3_d(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  REAL*8, DIMENSION(ni,nj,nk)         :: v
  REAL*8, DIMENSION(ni,nj,nk)         :: v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL exp_off_diag_v_m_2_d(v(:,:,k),                                     &
                               v_scr(:,:,k),                                 &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_3_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_off_diagonal_mul_z
!***begin prologue     so_exp_off_diagonal_mul_z     
!***date written       040607   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            compute the effect of the split operator
!                      off-diagonal, exponential propagator on a vector.
!
!***description
!***references
!***routines called    
!***                  
!
!***end prologue       so_exp_off_diagonal_mul_z
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE so_exp_off_diagonal_mul_1d_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(nphy(1))                   :: wave_function
     COMPLEX*16, DIMENSION(nphy(1))                   :: scratch_vector
     call exp_off_diag_m_v_z(wave_function,                      &
                             scratch_vector,                     &
                             nphy(1),1,1)
  END SUBROUTINE so_exp_off_diagonal_mul_1d_z
  SUBROUTINE so_exp_off_diagonal_mul_2d_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(nphy(2),nphy(1))           :: wave_function
     COMPLEX*16, DIMENSION(nphy(2),nphy(1))           :: scratch_vector
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(2),
!            that is H(2,2) * V(2,1) 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v_z(wave_function,                     &
                              scratch_vector,                    &
                              nphy(2),nphy(1),2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!            This code is for the general problem involving nphy(1),
!            V(2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_v_m_2_z(wave_function,                   &
                                scratch_vector,                  &
                                nphy(2),nphy(1),1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_2d_z
  SUBROUTINE so_exp_off_diagonal_mul_3d_z(wave_function,scratch_vector)
     IMPLICIT NONE
     COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))   :: wave_function
     COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))   :: scratch_vector

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            This code is for the general problem involving nphy(3),
!            H(3,3) * V(3,2,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call exp_off_diag_m_v_z(wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2)*nphy(1),3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general case involving the matrix multiply 
!        involving nphy(2), that is, V(3,2,1) * H(2,2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_3_z                                  &
                             (wave_function,                     &
                              scratch_vector,                    &
                              nphy(3),nphy(2),nphy(1),2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        This is for the general problem involving the matrix multiply 
!        for nphy(1), that is, V(3,2,1) * H(1,1)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call exp_off_diag_v_m_2_z                                  &
                             (wave_function,                     &
                              scratch_vector,                    &
                              nphy(3)*nphy(2),nphy(1),1)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END SUBROUTINE so_exp_off_diagonal_mul_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck exp_off_diag_m_v_d_z
!***begin prologue     exp_off_diag_m_v_d_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which 
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = exp_t_mat * V
!
!------------------------------------------------------------------------------------
!***                   where exp_t_mat is the regional
!***                   propagator for a FEDVR or FD Hamiltonian.
!
!***                   The parameter nj is a dummy and can take on values
!***                   consistent with any dimensional problem.  In a 1D
!***                   case nj=1, in 2D nj=nx and in 3D_D nj=ny*nx.

!***references
!***routines called    v_m_v_gen, v_m_v_2, v_m_v_3
!
!***end prologue       exp_off_diag_m_v_z                     
!
  SUBROUTINE exp_off_diag_m_v_z(v,                                             &
                                v_scr,                                         &
                                ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  INTEGER                                  :: i, ntrips, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
!  The sequence is called three times.  First for the odd matrices, then
!  for the even matrices and then again for the odd matrices.
! 
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips=3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
  DO trips=1,ntrips
     DO i = begin(trips), num_reg(index), 2
!        write(iout,*) i, nfun_reg(i,index), locate(trips)
!
        IF ( nfun_reg(i,index) > 10 ) then
!
!                       General Code
!
           CALL v_so_mat_v_gen(v(locate(trips):,:),                               &
                               v_scr(locate(trips):,:),                           &
                               mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                               ni,nj,nfun_reg(i,index))
        ELSE IF (nfun_reg(i,index) == 2) then
! 
!                       Special case for Two by Two
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           CALL v_so_mat_v_2(v(locate(trips):,:),                                &
                             v_scr(locate(trips):,:),                            &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),       &
                             ni,nj)
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Special case for Three by Three
!
!           write(iout,*) 'Input Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
           CALL v_so_mat_v_3(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
!           write(iout,*) 'Output Arrays'
!           write(iout,*) v(locate(trips):,:)
!           write(iout,*) v_scr(locate(trips):,:)
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Special case for Four by Four
!
           CALL v_so_mat_v_4(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Special case for Five by Five
!
           CALL v_so_mat_v_5(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Special case for Six by Six
!
           CALL v_so_mat_v_6(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Special case for Seven by Seven
!
           CALL v_so_mat_v_7(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Special case for Eight by Eight
!
           CALL v_so_mat_v_8(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Special case for Nine by Nine
!
           CALL v_so_mat_v_9(v(locate(trips):,:),                               &
                             v_scr(locate(trips):,:),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj)
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Special case for Ten by Ten
!
           CALL v_so_mat_v_10(v(locate(trips):,:),                              &
                              v_scr(locate(trips):,:),                          &
                              mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                              ni,nj)
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                       &
                                          +                                         &
                            nfun_reg(i+1,index)                                     &
                                          -                                         &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_m_v_z
!
!*deck exp_off_diag_v_m_2_z
!***begin prologue     exp_off_diag_v_m_2_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine is the driver routine which
!***                   calculates,
!------------------------------------------------------------------------------------
!
!                 V = V * exp_t_mat 
!
!------------------------------------------------------------------------------------
!***                   where exp_t_mat is the regional
!***                   propagator for a FEDVR or FD Hamiltonian. 
!***                   parts.  
!***                   The routine is not needed in 1D problems. 
!***                   In 2D, the ni index runs over the number
!                      of points in the y coordinate while in 3D_D,
!***                   it runs over the product nz*ny of the number 
!***                   of points in the z and y coordinates.
!
!***references
!***routines called    v_v_m_gen, v_v_m_2, v_v_m_3
!***end prologue       exp_off_diag_v_m_2_z
!
  SUBROUTINE exp_off_diag_v_m_2_z(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,index)
  IMPLICIT NONE
  INTEGER                                  :: ni, nj, index, ntrips
  COMPLEX*16, DIMENSION(ni,nj)             :: v
  COMPLEX*16, DIMENSION(ni,nj)             :: v_scr
  INTEGER                                  :: i, trips
  INTEGER, DIMENSION(3)                    :: locate, begin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  locate(1) = 1
  locate(2) = nfun_reg(1,index)
  locate(3) = 1
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  ntrips = 3
  IF(num_reg(index) == 1 ) THEN
     ntrips=1
  END IF
!
  DO trips = 1, ntrips
     DO i=begin(trips), num_reg(index), 2
        IF ( nfun_reg(i,index) > 10 ) then
!
!                        General Case
!
           CALL v_v_so_mat_gen(v(:,locate(trips):),                            &
                               v_scr(:,locate(trips):),                        &
                               mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),   &
                               ni,nj,nfun_reg(i,index)) 
!
        ELSE IF (nfun_reg(i,index) == 2) then
!
!                       Two by Two
! 
           CALL v_v_so_mat_2(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                             ni,nj) 
       ELSE IF (nfun_reg(i,index) == 3) then
!
!                       Three by Three
!
           CALL v_v_so_mat_3(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 4) then
!
!                       Four by Four
!
           CALL v_v_so_mat_4(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 5) then
!
!                       Five by Five
!
           CALL v_v_so_mat_5(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 6) then
!
!                       Six by Six
!
           CALL v_v_so_mat_6(v(:,locate(trips):),                              &
                             v_scr(:,locate(trips):),                          &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 7) then
!
!                       Seven by Seven
!
           CALL v_v_so_mat_7(v(:,locate(trips):),                               &
                             v_scr(:,locate(trips):),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 8) then
!
!                       Eight by Eight
!
           CALL v_v_so_mat_8(v(:,locate(trips):),                               &
                             v_scr(:,locate(trips):),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 9) then
!
!                       Nine by Nine
!
           CALL v_v_so_mat_9(v(:,locate(trips):),                               &
                             v_scr(:,locate(trips):),                           &
                             mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),      &
                             ni,nj) 
        ELSE IF (nfun_reg(i,index) == 10) then
!
!                       Ten by Ten
!
           CALL v_v_so_mat_10(v(:,locate(trips):),                              &
                              v_scr(:,locate(trips):),                          &
                              mat_reg(i,index)%exp_t_mat_z(:,:,prop_point),     &
                              ni,nj) 
        END IF
        IF (i /= num_reg(index) ) then
            locate(trips) = locate(trips) + nfun_reg(i,index)                     &
                                          +                                       &
                            nfun_reg(i+1,index)                                   &
                                          -                                       &
                                          2
        END IF
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_2_z
!*deck exp_off_diag_v_m_3_z
!***begin prologue     exp_off_diag_v_m_3_z     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   exp_off_diag_v_m_2 in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       exp_off_diag_v_m_3_z
!
  SUBROUTINE exp_off_diag_v_m_3_z(v,                                         &
                                  v_scr,                                     &
                                  ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  COMPLEX*16, DIMENSION(ni,nj,nk)     :: v
  COMPLEX*16, DIMENSION(ni,nj,nk)     :: v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL exp_off_diag_v_m_2_z(v(:,:,k),                                     &
                               v_scr(:,:,k),                                 &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_3_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END MODULE split_operator_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
