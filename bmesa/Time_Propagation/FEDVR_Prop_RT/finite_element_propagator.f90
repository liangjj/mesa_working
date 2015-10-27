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
          MODULE PROCEDURE so_prop_1d_z, so_prop_2d_z, so_prop_3d_z
                  END INTERFACE so_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck so_prop_1d_z
!**begin prologue     so_prop_1d_z
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
!**end prologue       so_prop_1d_z
  SUBROUTINE so_prop_1d_z(wave_function,scratch_vector)
  USE dvrprop_global_rt
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)            :: wave_function
  REAL*8, DIMENSION(nphy(1),2)            :: scratch_vector
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, t, tst
  INTEGER                                 :: length, len, count
  REAL*8                                  :: norm, rtemp
  COMPLEX*16                              :: auto
  REAL*8                                  :: total_t
!
  count=0
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
      CALL cp_psi(wave_function,scratch_vector,t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL rsp(wave_function,scratch_vector)  
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)    
      call iosys('write real solution to bec',2*n3d,wave_function,0,' ')
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( count == plot_step ) THEN
        count=0
        IF(i0stat == 'gaussian-pulse') THEN
           call calc_moment(wave_function,tim_pts(t+1))
        ELSE IF(i0stat == 'perturbed-state-vector') THEN
!           CALL evolve(wave_function,scratch_vector,tim_pts(t+1))
        END IF
!       Form the total solution and compute the autocorrelation
!       function.
        DO i=1,nphy(1)
           scratch_vector(i,:) = wave_function(i,:) /sqrt(grid(1)%wt(i))
        END DO        
        CALL plot_psi(scratch_vector,t)
        CALL auto_correlation_function(auto,scratch_vector,     &
                                       wave_function)
        rtemp = 1.d0 - auto * conjg(auto)
        write(iout,3) tim_pts(t+1), rtemp
        call check_norm(wave_function,norm)
     END IF
     count=count+1
     time(6)=secnds(0.0)
     delta(5)=time(6)-time(5)
     total_t = delta(1) + delta(2) + delta(3) + delta(4) + delta(5)
     WRITE(iout,2) t, (delta(i), i=1,5), total_t
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i4,    &
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
END SUBROUTINE so_prop_1d_z
!deck so_prop_2d_z
!**begin prologue     so_prop_2d_z
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
!**end prologue       so_prop_2d_z
  SUBROUTINE so_prop_2d_z(wave_function,scratch_vector)
  USE dvrprop_global_rt
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)    :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1),2)    :: scratch_vector
  CHARACTER (LEN=2)                       :: itoc
  INTEGER                                 :: i, j, t, tst
  INTEGER                                 :: length, len
  REAL*8                                  :: norm, rtemp
  REAL*8                                  :: wt1, wt2
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
      CALL cp_psi(wave_function,scratch_vector,t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL rsp(wave_function,scratch_vector)    
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)    
      call iosys('write real solution to bec',2*n3d,wave_function,0,' ')
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( tst == 0 ) THEN
        IF(i0stat == 'gaussian-pulse') THEN
           call calc_moment(wave_function,tim_pts(t+1),f_2)
        END IF
!       Form the total solution and compute the autocorrelation
!       function.
        DO i=1,nphy(2)
           wt2=1.d0/sqrt(grid(2)%wt(i))            
           DO j=1,nphy(1)
              wt1=1.d0/sqrt(grid(1)%wt(j))            
              scratch_vector(i,j,:) = wave_function(i,j,:) * wt1 * wt2
           END DO        
        END DO
        CALL plot_psi(scratch_vector,t)
        CALL auto_correlation_function(auto,scratch_vector,     &
                                       wave_function,f_2)
        rtemp = 1.d0 - auto * conjg(auto)
        write(iout,3) tim_pts(t+1), rtemp
        call check_norm(wave_function,norm,f_1)
     END IF
     time(6)=secnds(0.0)
     delta(5)=time(6)-time(5)
     total_t = delta(1) + delta(2) + delta(3) + delta(4) + delta(5)
     WRITE(iout,2) t, (delta(i), i=1,5), total_t
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i4,    &
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
END SUBROUTINE so_prop_2d_z
!deck so_prop_3d_z
!**begin prologue     so_prop_3d_z
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
!**end prologue       so_prop_3d_z
  SUBROUTINE so_prop_3d_z(wave_function,scratch_vector)
  USE dvrprop_global_rt
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  USE non_linear_potential
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)    :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)    :: scratch_vector
  CHARACTER (LEN=2)                               :: itoc
  INTEGER                                         :: i, j, k, t, tst
  INTEGER                                         :: iostat, length, len
  REAL*8                                          :: norm, rtemp
  COMPLEX*16                                      :: auto
  REAL*8                                          :: t_int, total_t
  REAL*8                                          :: wt1, wt2, wt3
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
      CALL cp_psi(wave_function,scratch_vector,t,tim_pts(t))
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
!  
!     Calculate the non-linear potential if present.
!  
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl(wave_function)
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      CALL rsp(wave_function,scratch_vector)    
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)    
      call iosys('write real solution to bec',2*n3d,wave_function,0,' ')
      tst = t - plot_step * ( t /plot_step )
!
!  
!        compare the approximate and exact solutions where possible
!
     IF( tst == 0 ) THEN
        IF(i0stat == 'gaussian-pulse') THEN
           call calc_moment(wave_function,tim_pts(t+1),f_2,f_3)
        END IF
!       Form the total solution and compute the autocorrelation
!       function.
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
        CALL plot_psi(scratch_vector,t)
        CALL auto_correlation_function(auto,scratch_vector,     &
                                       wave_function,f_2,f_3)
        rtemp = 1.d0 - auto * conjg(auto)
        write(iout,3) tim_pts(t+1), rtemp
        call check_norm(wave_function,norm,f_1,f_2)
     END IF
     time(6)=secnds(0.0)
     delta(5)=time(6)-time(5)
     total_t = delta(1) + delta(2) + delta(3) + delta(4) + delta(5)
     WRITE(iout,2) t, (delta(i), i=1,5), total_t
  END DO
!
!     End the propagation
!
  CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
2 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time Summary for Interval         = ',i4,    &
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
END SUBROUTINE so_prop_3d_z
END MODULE finite_element_propagator
