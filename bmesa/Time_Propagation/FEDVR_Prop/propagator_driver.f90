!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE propagator_driver
!
!           Generic Interfac for propagator routines
!
  INTERFACE so_prop
     MODULE PROCEDURE so_prop_1d, so_prop_2d, so_prop_3d
  END INTERFACE so_prop
  CONTAINS
!deck so_prop_1d
!**begin prologue     so_prop_1d
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
!**end prologue       so_prop_1d
  SUBROUTINE so_prop_1d(psi,v_scr,n_d)
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                 :: n_d
  REAL*8                                  :: psi(n_d,2), v_scr(n_d,2)
  n3d=n_d
  ALLOCATE(sin_diag(n3d),cos_diag(n3d),v_tot(n3d),tim_pts(ntreg+1))
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     key='FEDVR'
  ELSE IF (typke == 'fd' ) then
     key='FD'
  END IF
  CALL open_plot
  DEALLOCATE(sin_diag,cos_diag,v_tot,tim_pts)
  END SUBROUTINE so_prop_1d
!deck so_prop_2d
!**begin prologue     so_prop_2d
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
!**end prologue       so_prop_2d
  SUBROUTINE so_prop_2d(psi,v_scr,n_d)
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                   :: n_d
  INTEGER                                 :: words
  REAL*8, DIMENSION(n_d(2),n_d(1),2)      :: psi
  REAL*8, DIMENSION(n_d(2),n_d(1),2)      :: v_scr
  n3d = n_d(2) * n_d(1)
  ALLOCATE(sin_diag(n3d),cos_diag(n3d),v_tot(n3d),tim_pts(ntreg+1))
  maxdim=max(n_d(2),n_d(1))
  words=max(n_d(1)*4,maxdim**2)
  ALLOCATE(fac(words))
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     key='FEDVR'
  ELSE IF (typke == 'fd' ) then
     key='FD'
  END IF
  CALL open_plot
  DEALLOCATE(sin_diag,cos_diag,v_tot,tim_pts)
  DEALLOCATE(fac)
  END SUBROUTINE so_prop_2d
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
  SUBROUTINE so_prop_3d(psi,v_scr,n_d)
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER, DIMENSION(3)                      :: n_d
  INTEGER, DIMENSION(2)                      :: words
  REAL*8, DIMENSION(n_d(3),n_d(2),n_d(1),2)  :: psi
  REAL*8, DIMENSION(n_d(3),n_d(2),n_d(1),2)  :: v_scr
  n3d = n_d(3) * n_d(2) * n_d(1)
  ALLOCATE(sin_diag(n3d), cos_diag(n3d), v_tot(n3d), tim_pts(ntreg+1))
  maxdim=max(n_d(3),n_d(2),n_d(1))
  words(1)=max(n_d(2)*n_d(1)*4,maxdim*maxdim*2)
  words(2)=max(n_d(1)*4,maxdim*2)
  ALLOCATE(fac(words(1)),v_1(words(2)))
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     key='FEDVR'
  ELSE IF (typke == 'fd' ) then
     key='FD'
  END IF
  CALL open_plot
  DEALLOCATE(sin_diag,cos_diag,v_tot,tim_pts)
  DEALLOCATE(fac,v_1)
  END SUBROUTINE so_prop_3d
!deck open_plot
!**begin prologue     open_plot
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, real-space, split-operator, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             FEDVR_Prop
!**purpose            open plotting routines.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       open_plot
  SUBROUTINE open_plot
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                 :: t, i, iostat
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
     OPEN (UNIT=iplot(2),FILE='time_points',     &
           ACCESS='sequential',FORM='formatted', &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(3),FILE='real_wavefunctions',     &
           ACCESS='sequential',FORM='formatted', &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(4),FILE='imaginary_wavefunction',     &
           ACCESS='sequential',FORM='formatted', &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(5),FILE='correlation_function',     &
           ACCESS='sequential',FORM='formatted', &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     write(iplot(2),*) (tim_pts(i),i=1,ntreg+1)
  END IF
  END SUBROUTINE open_plot
END MODULE propagator_driver

