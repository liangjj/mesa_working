!***********************************************************************
                           MODULE trial_vectors
                           INTERFACE trial
                    MODULE PROCEDURE trial_d, trial_z
                       END INTERFACE trial
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck trial_z
!**begin prologue     trial_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, trial
!***
!**author             schneider, b. i.(nsf)
!**source             trial
!**purpose            trial vectors for  Arnoldi method.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       trial_z
  SUBROUTINE trial_z(v_in)
  USE arnoldi_global_rt
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)            :: v_in
  INTEGER                               :: i, nout 
!
!  First vector is taken as solution from previous step.  The space is
!  completed with unit vectors and the entire set schmidt orthonormalized.
!
  vec(:,1) = v_in
  IF(ntrial > 1) THEN
     vec(:,2:ntrial)=(0.d0,0.d0)
      DO  i=2,ntrial
          vec(i,i)=(1.d0,0.d0)
      END DO
  END IF
  CALL cschmt(vec,thresh,n3d,1,ntrial,nout,.true.,.false.)
END SUBROUTINE trial_z
!deck trial_d
!**begin prologue     trial_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, trial
!***
!**author             schneider, b. i.(nsf)
!**source             trial
!**purpose            trial vectors for  Arnoldi method.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       trial_d
  SUBROUTINE trial_d(v_in)
  USE arnoldi_global_it
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)            :: v_in
  INTEGER                           :: i, nout 
!
!  First vector is taken as solution from previous step.  The space is
!  completed with unit vectors and the entire set schmidt orthonormalized.
!
  vec(:,1) = v_in
  IF(ntrial > 1) THEN
     vec(:,2:ntrial)=0.d0
      DO  i=2,ntrial
          vec(i,i)=1.d0
      END DO
  END IF
  CALL dschmt(vec,thresh,n3d,1,ntrial,nout,.true.,.false.)
END SUBROUTINE trial_d
END MODULE trial_vectors





















