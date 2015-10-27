!deck so_prop_3d
!**begin prologue     so_prop_3d
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
!**end prologue       so_prop_3d
  SUBROUTINE so_prop_3d(wave_function,scratch_vector)
  USE dvrprop_global_rt
  USE dvr_global
  USE dvr_shared
  USE normalize
  USE moment
  USE auto_correlation
  USE plot_wavefunction
  USE initial_state
  USE real_space_propagator
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)    :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)    :: scratch_vector
  CHARACTER (LEN=2)                               :: itoc
  INTEGER                                         :: i, t, tst
  INTEGER                                         :: iostat, length, len
  REAL*8                                          :: norm, rtemp
  COMPLEX*16                                      :: auto
  REAL*8                                          :: t_int, total_t
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
         CALL v_nl
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
!
!     Since the potential is, in general, time-dependent, the diagonal
!     propagator must be recomputed each time.
!
!      The input initial state is psi.  The output is chi
!
      IF(prop_order == 2 ) then
         CALL rsp(wave_function,scratch_vector)    
         time(5)=secnds(0.0)
         delta(4)=time(5)-time(4)    
      ELSE
          CALL lnkerr('quit.  error in propagation order')
      END IF
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
        CALL plot_psi(wave_function,t)
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
END SUBROUTINE so_prop_3d
