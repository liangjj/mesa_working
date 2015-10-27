!deck setup
!**begin prologue     setup
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
!**end prologue       setup
  SUBROUTINE setup
  USE dvrprop_global
  USE dvr_global,               ONLY : spdim
  USE dvr_shared,               ONLY : typke, grid
  USE moment
  IMPLICIT NONE
  INTEGER                           :: t, i, iostat
!
  ALLOCATE(tim_pts(ntreg+1))
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     key='FEDVR'
  ELSE IF (typke == 'fd' ) then
     key='FD'
  END IF
!
!     Set up a file to hold the wavefunction for plotting
!
  call iosys('open plot as scratch',0,0,0,' ')
  call iosys('create real wavefunction on plot',ntreg*2*n3d,0,0,' ')
  tim_pts(1)=t_init
  DO t=2,ntreg+1
     tim_pts(t)=tim_pts(t-1) + deltat
  END DO
  IF(plot) then
     DO i=1,spdim
        write(iplot(1),*) grid(i)%pt
     END DO
     OPEN (UNIT=iplot(2),FILE='time_points',             &
           ACCESS='sequential',FORM='formatted',         &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(3),FILE='real_wavefunctions',      &
           ACCESS='sequential',FORM='formatted',         &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(4),FILE='imaginary_wavefunction',  &
           ACCESS='sequential',FORM='formatted',         &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(5),FILE='correlation_function',    &
           ACCESS='sequential',FORM='formatted',         &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     write(iplot(2),*) (tim_pts(i),i=1,ntreg+1)
  END IF
  IF(i0stat == 'gaussian-pulse') THEN
     CALL moment_data
  END IF
  WRITE(iout,1) spdim, key, n3d
!
!
1 Format(/,1x,'Entering Propagation for a ',i1,' dimensional problem' &
         /1x, 'Spatial discretization is ',a8,1x,'Number of Points = ',i10)  
END SUBROUTINE setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
