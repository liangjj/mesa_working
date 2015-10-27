!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     MODULE finite_element_propagator
!**begin prologue     finite_element_propagator
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       finite_element_propagator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE so_prop
             MODULE PROCEDURE so_prop_1d_d, so_prop_2d_d, so_prop_3d_d
                 END INTERFACE so_prop
                     INTERFACE arn_prop
             MODULE PROCEDURE arn_prop_1d_d, arn_prop_2d_d, arn_prop_3d_d
                 END INTERFACE arn_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_prop_1d_d
!**begin prologue     so_prop_1d_d
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
!**end prologue       so_prop_1d_d
  SUBROUTINE so_prop_1d_d(wave_function,scratch_vector)
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),nvec)         :: wave_function
  REAL*8, DIMENSION(nphy(1),nvec)         :: scratch_vector
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, t, tst
  INTEGER                                 :: length, len
  REAL*8                                  :: norm
  COMPLEX*16                              :: auto
  REAL*8                                  :: total_t
!
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      time(1)=secnds(0.0)
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
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.  The first time through
!     the initial state is stored in soln which is needed later.
!
      CALL cp_psi(wave_function(:,1),t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function(:,1))
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!
      CALL rsp(wave_function(:,1),scratch_vector(:,1))  
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)    
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     call check_norm(wave_function(:,1),norm)
     scratch_vector = wave_function/sqrt(norm) 
     DO i=1,nphy(1)
        scratch_vector(i,:) = scratch_vector(i,:)/sqrt(grid(1)%wt(i)) 
     END DO
     IF( tst == 0 ) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(scratch_vector(:,1),t)
     END IF
     IF (e_method == 'exponential') then
         CALL eigen_solve(scratch_vector(:,1),wave_function(:,1))
     ELSE
         call check_norm(wave_function(:,1),norm)
         wave_function = wave_function/sqrt(norm) 
         call h_v(wave_function,scratch_vector,nvec)
         call check_energy(wave_function(:,1),scratch_vector(:,1))
     END IF
     call iosys('write real solution to bec',n3d,                    &
                 wave_function,0,' ')
     time(6)=secnds(0.0)
     delta(5)=time(6)-time(5)
     total_t = delta(1) + delta(2) + delta(3) + delta(4) + delta(5)
     WRITE(iout,2) t, (delta(i), i=1,5), total_t
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function(:,1),norm)
  wave_function = wave_function/sqrt(norm)
  DO i=1,nphy(1)
     scratch_vector(i,:) = wave_function(i,:)/sqrt(grid(1)%wt(i)) 
  END DO
  CALL plot_psi(scratch_vector(:,1),t-1) 
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i8,    &
         /,10X,'time for linear perturbation      = ',f15.8, &
         /,10X,'time for right hand side          = ',f15.8, &
         /,10X,'time for non-linear perturbation  = ',f15.8, &
         /,10X,'time to propagate                 = ',f15.8, &
         /,10X,'time for plots, moments and norms = ',f15.8, &
         /,10X,'Total                             = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
3 FORMAT(/5x,'Survival Probablity',/,/5x,'Time = ', e15.8, &
          5x,'Probability = ', e15.8)
END SUBROUTINE so_prop_1d_d
!deck so_prop_2d_d
!**begin prologue     so_prop_2d_d
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
!**end prologue       so_prop_2d_d
  SUBROUTINE so_prop_2d_d(wave_function,scratch_vector)
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE eigenvalues
  USE psi_h_psi
  USE h_on_vector
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)    :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)    :: scratch_vector
  CHARACTER (LEN=2)                          :: itoc
  INTEGER                                    :: i, j, t, tst
  INTEGER                                    :: length, len
  REAL*8                                     :: norm, rtemp
  COMPLEX*16                                 :: auto
  REAL*8                                     :: total_t, wt1, wt2
!
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      time(1)=secnds(0.0)
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
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.  The first time through
!     the initial state is stored in soln which is needed later.
!
      CALL cp_psi(wave_function(:,:,1),t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function(:,:,1))
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL rsp(wave_function(:,:,1),scratch_vector(:,:,1))    
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)    
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     call check_norm(wave_function(:,:,1),norm,f_1)
     scratch_vector = wave_function/sqrt(norm) 
     DO i=1,nphy(2)
        wt2=1.d0/sqrt(grid(2)%wt(i)) 
        DO j=1,nphy(1)
           wt1=1.d0/sqrt(grid(1)%wt(j)) 
           scratch_vector(i,j,:) = scratch_vector(i,j,:) * wt1 * wt2 
        END DO
     END DO
     IF( tst == 0 ) THEN
        CALL plot_psi(wave_function(:,:,1),t)
     END IF
     IF (e_method == 'exponential') then
         CALL eigen_solve(scratch_vector(:,:,1),             &
                          wave_function(:,:,1),f_1)
     ELSE
         call check_norm(wave_function(:,:,1),norm,f_1)
         wave_function = wave_function/sqrt(norm) 
         call h_v(wave_function,scratch_vector,nvec)
         call check_energy(wave_function(:,:,1),             &
                           scratch_vector(:,:,1),f_1)
     END IF
     call iosys('write real solution to bec',n3d,            &
                 wave_function,0,' ')
     time(6)=secnds(0.0)
     delta(5)=time(6)-time(5)
     total_t = delta(1) + delta(2) + delta(3) + delta(4) + delta(5)
     WRITE(iout,2) t, (delta(i), i=1,5), total_t
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function(:,:,1),norm,f_1)
  wave_function = wave_function/sqrt(norm) 
  DO i=1,nphy(2)
     wt2=1.d0/sqrt(grid(2)%wt(i)) 
     DO j=1,nphy(1)
        wt1=1.d0/sqrt(grid(1)%wt(j)) 
        scratch_vector(i,j,:) = wave_function(i,j,:) * wt1 * wt2 
     END DO
  END DO
  CALL plot_psi(scratch_vector(:,:,1),t-1)
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i8,    &
         /,10X,'time for linear perturbation      = ',f15.8, &
         /,10X,'time for right hand side          = ',f15.8, &
         /,10X,'time for non-linear perturbation  = ',f15.8, &
         /,10X,'time to propagate                 = ',f15.8, &
         /,10X,'time for plots, moments and norms = ',f15.8, &
         /,10X,'Total                             = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
END SUBROUTINE so_prop_2d_d
!deck so_prop_3d_d
!**begin prologue     so_prop_3d_d
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
!**end prologue       so_prop_3d_d
  SUBROUTINE so_prop_3d_d(wave_function,scratch_vector)
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE eigenvalues
  USE psi_h_psi
  USE h_on_vector
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)       :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)       :: scratch_vector
  CHARACTER (LEN=2)                                     :: itoc
  INTEGER                                               :: i, j, k, t, tst
  INTEGER                                               :: iostat, length, len
  REAL*8                                                :: norm, rtemp
  COMPLEX*16                                            :: auto
  REAL*8                                                :: t_int, total_t
  REAL*8                                                :: wt1, wt2, wt3
!
!
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      time(1)=secnds(0.0)
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
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.  The first time through
!     the initial state is stored in soln which is needed later.
!
      CALL cp_psi(wave_function(:,:,:,1),t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function(:,:,:,1))
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL rsp(wave_function(:,:,:,1),scratch_vector(:,:,:,1)) 
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)    
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     call check_norm(wave_function(:,:,:,1),norm,f_1,f_2)
     scratch_vector = wave_function/sqrt(norm) 
     DO i=1,nphy(3)
        wt3=1.d0/sqrt(grid(3)%wt(i)) 
        DO j=1,nphy(2)
           wt2=1.d0/sqrt(grid(2)%wt(j)) 
           DO k=1,nphy(1)
              wt1=1.d0/sqrt(grid(1)%wt(k)) 
              scratch_vector(i,j,k,:) = scratch_vector(i,j,k,:) * wt1 * wt2 * wt3
           END DO
       END DO
     END DO
     IF( tst == 0 ) THEN
        CALL plot_psi(wave_function(:,:,:,1),t)
     END IF
     IF (e_method == 'exponential') then
         CALL eigen_solve(scratch_vector(:,:,:,1),             &
                          wave_function(:,:,:,1),f_1,f_2)
     ELSE
         call check_norm(wave_function(:,:,:,1),norm,f_1,f_2)
         wave_function = wave_function/sqrt(norm) 
         call h_v(wave_function,scratch_vector,nvec)
         call check_energy(wave_function(:,:,:,1),             &
                           scratch_vector(:,:,:,1),f_1,f_2)
     END IF
     call iosys('write real solution to bec',n3d,         &
                 wave_function,0,' ')
     time(6)=secnds(0.0)
     delta(5)=time(6)-time(5)
     total_t = delta(1) + delta(2) + delta(3) + delta(4) + delta(5)
     WRITE(iout,2) t, (delta(i), i=1,5), total_t
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function(:,:,:,1),norm,f_1,f_2)
  wave_function = wave_function/sqrt(norm) 
  DO i=1,nphy(3)    
     wt3=1.d0/sqrt(grid(3)%wt(i)) 
     DO j=1,nphy(2)
        wt2=1.d0/sqrt(grid(2)%wt(j)) 
        DO k=1,nphy(1)
           wt1=1.d0/sqrt(grid(1)%wt(k)) 
           scratch_vector(i,j,k,:) = wave_function(i,j,k,:) * wt1 * wt2 * wt3
        END DO
    END DO
  END DO
  CALL plot_psi(scratch_vector(:,:,:,1),t-1)
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i8,    &
         /,10X,'time for linear perturbation      = ',f15.8, &
         /,10X,'time for right hand side          = ',f15.8, &
         /,10X,'time for non-linear perturbation  = ',f15.8, &
         /,10X,'time to propagate                 = ',f15.8, &
         /,10X,'time for plots, moments and norms = ',f15.8, &
         /,10X,'Total                             = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
3 FORMAT(/5x,'Survival Probablity',/,/5x,'Time = ', e15.8, &
          5x,'Probability = ', e15.8)
END SUBROUTINE so_prop_3d_d
!deck arn_prop_1d_d
!**begin prologue     arn_prop_1d_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Arnoldi_Prop_f90
!**purpose            time dependent schroedinger equation
!**                   using Arnoldi method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       arn_prop_1d_d_real
  SUBROUTINE arn_prop_1d_d(wave_function,scratch_vector)
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),nvec)         :: wave_function
  REAL*8, DIMENSION(nphy(1),nvec)         :: scratch_vector
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, t, tst
  INTEGER                                 :: bigs, bigv
  INTEGER                                 :: length, len
  REAL*8                                  :: t_int, total_t
  REAL*8                                  :: norm
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
      deltat=t1-t0
      time(1)=secnds(0.0)
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      CALL v_tim
      CALL pert
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi(wave_function(:,1),t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function(:,1))
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Arnoldi code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
      CALL arnoldi(wave_function)
!  
!        compare the approximate and exact solutions where possible
     tst = t - plot_step * ( t /plot_step )  
     IF(tst == 0) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function(:,1),t)
     END IF
! Form the total solution and compute the autocorrelation
! function.
  
     time(5)=secnds(0.0)
     delta(4)=time(5)-time(4)
     total_t = delta(1) + delta(2) + delta(3) + delta(4)
     WRITE(iout,2) t, (delta(i), i=1,4), total_t
     call check_norm(wave_function(:,1),norm)
   END DO
!
!     End the propagation
!
   CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
   call check_norm(wave_function(:,1),norm)
  wave_function = wave_function/sqrt(norm)
  CALL plot_psi(wave_function(:,1),t-1) 
!
1 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time to pack the spatial Hamiltonian = '     &
               ,f15.8,/                                      &
         '***********************************************'   &
         '*************************')
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i4,    &
         /,10X,'time for linear perturbation      = ',f15.8,   &
         /,10X,'time for right hand side          = ',f15.8,   &
         /,10X,'time for non-linear perturbation  = ',f15.8,   &
         /,10X,'time to propagate                 = ',f15.8,/, &
         /,10X,'Total                             = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
3 FORMAT(/5x,'Survival Probablity',/,/15x,'Time', &
          15x,'Probability')
4 FORMAT(/10x,e15.8,5x,e15.8)
END SUBROUTINE arn_prop_1d_d
!deck arn_prop_2d_d
!**begin prologue     arn_prop_2d_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Arnoldi_Prop_f90
!**purpose            time dependent schroedinger equation
!**                   using Arnoldi method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       arn_prop_2d_d
  SUBROUTINE arn_prop_2d_d(wave_function,scratch_vector)
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)   :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)   :: scratch_vector
  CHARACTER (LEN=2)                         :: itoc
  INTEGER                                   :: i, t, tst
  INTEGER                                   :: bigs, bigv
  INTEGER                                   :: length, len
  REAL*8                                    :: t_int, total_t
  REAL*8                                    :: norm
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
      deltat=t1-t0
      time(1)=secnds(0.0)
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      CALL v_tim
      CALL pert
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi(wave_function(:,:,1),t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function(:,:,1))
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Arnoldi code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
      CALL arnoldi
!  
!        compare the approximate and exact solutions where possible
     tst = t - plot_step * ( t /plot_step )  
     IF(tst == 0) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function(:,:,1),t)
     END IF
! Form the total solution and compute the autocorrelation
! function.
  
     time(5)=secnds(0.0)
     delta(4)=time(5)-time(4)
     total_t = delta(1) + delta(2) + delta(3) + delta(4)
     WRITE(iout,2) t, (delta(i), i=1,4), total_t
     call check_norm(wave_function(:,:,1),norm,f_1)
   END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function(:,:,1),norm,f_1)
  wave_function = wave_function/sqrt(norm)
  CALL plot_psi(wave_function(:,:,1),t-1) 
!
1 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time to pack the spatial Hamiltonian = '     &
               ,f15.8,/                                      &
         '***********************************************'   &
         '*************************')
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i4,    &
         /,10X,'time for linear perturbation      = ',f15.8,   &
         /,10X,'time for right hand side          = ',f15.8,   &
         /,10X,'time for non-linear perturbation  = ',f15.8,   &
         /,10X,'time to propagate                 = ',f15.8,/, &
         /,10X,'Total                             = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
3 FORMAT(/5x,'Survival Probablity',/,/15x,'Time', &
          15x,'Probability')
4 FORMAT(/10x,e15.8,5x,e15.8)
END SUBROUTINE arn_prop_2d_d
!deck arn_prop_3d_d
!**begin prologue     arn_prop_3d_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Arnoldi_Prop_f90
!**purpose            time dependent schroedinger equation
!**                   using Arnoldi method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       arn_prop_3d_d
  SUBROUTINE arn_prop_3d_d(wave_function,scratch_vector)
  USE arnoldi_global_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE psi_h_psi
  USE h_on_vector
  USE eigenvalues
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)       &
                                                 :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)       &
                                                 :: scratch_vector
  CHARACTER (LEN=2)                              :: itoc
  INTEGER                                        :: i, t, tst
  INTEGER                                        :: bigs, bigv
  INTEGER                                        :: length, len
  REAL*8                                         :: t_int, total_t
  REAL*8                                         :: norm
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=tim_pts(t)
      t1=tim_pts(t+1)
      deltat=t1-t0
      time(1)=secnds(0.0)
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      CALL v_tim
      CALL pert
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
!  
!     Initialize the wavefunction at $t_{0}$ or 
!     read in its value from the disk at $t_{i-1}$.
!
      CALL cp_psi(wave_function(:,:,:,1),t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function(:,:,:,1))
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Arnoldi code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
      CALL arnoldi
!  
!        compare the approximate and exact solutions where possible
     tst = t - plot_step * ( t /plot_step )  
     IF(tst == 0) THEN
!
!       Form the total solution and compute the autocorrelation
!       function.
!
        CALL plot_psi(wave_function(:,:,:,1),t)
     END IF
! Form the total solution and compute the autocorrelation
! function.
  
     time(5)=secnds(0.0)
     delta(4)=time(5)-time(4)
     total_t = delta(1) + delta(2) + delta(3) + delta(4)
     WRITE(iout,2) t, (delta(i), i=1,4), total_t
     call check_norm(wave_function(:,:,:,1),norm,f_1,f_2)
   END DO
!
!     End the propagation
!
   CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
  call check_norm(wave_function(:,:,:,1),norm,f_1,f_2)
  wave_function = wave_function/sqrt(norm) 
        CALL plot_psi(wave_function(:,:,:,1),t-1)
!
1 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time to pack the spatial Hamiltonian = '     &
               ,f15.8,/                                      &
         '***********************************************'   &
         '*************************')
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i4,    &
         /,10X,'time for linear perturbation      = ',f15.8,   &
         /,10X,'time for right hand side          = ',f15.8,   &
         /,10X,'time for non-linear perturbation  = ',f15.8,   &
         /,10X,'time to propagate                 = ',f15.8,/, &
         /,10X,'Total                             = ',f15.8,/, &
         '***********************************************'     &
         '*************************')
3 FORMAT(/5x,'Survival Probablity',/,/15x,'Time', &
          15x,'Probability')
4 FORMAT(/10x,e15.8,5x,e15.8)
END SUBROUTINE arn_prop_3d_d
END MODULE finite_element_propagator
