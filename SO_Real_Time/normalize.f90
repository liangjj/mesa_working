!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE normalize
!**begin prologue     normalize
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       normalize
!
  IMPLICIT NONE
  CONTAINS
!deck check_norm_1d.f
!***begin prologue     check_norm_1d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for one dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_1d
  SUBROUTINE check_norm_1d(v)
  USE dvrprop_global
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: v
  REAL*8                                 :: ddot
  REAL*8                                 :: norm, first, last
  INTEGER                                :: i
!
  norm = ddot(nphy(1),v(1,1),1,v(1,1),1)               &
                         +                             &
         ddot(nphy(1),v(1,2),1,v(1,2),1)
  IF( typke == 'fd' ) then
      first = v(1,1) * v(1,1) + v(1,2) * v(1,2)
      last  = v(nphy(1),1) * v(nphy(1),1)              &
                          +                            &
              v(nphy(1),2) * v(nphy(1),2)
      norm = norm - .5d0 * ( first + last)
      norm = del * norm
  END IF
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_1d
!
!
!deck check_norm_2d.f
!***begin prologue     check_norm_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for two dimension 
!                      wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_2d
  SUBROUTINE check_norm_2d(v,f_2)
  USE dvrprop_global
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(1))             :: f_2
  REAL*8                                 :: ddot
  REAL*8                                 :: nrm, first, last
  INTEGER                                :: i, j
!
  nrm = 0.d0
  DO i=1,nphy(1)
     f_2(i) = ddot(nphy(2),v(1,i,1),1,v(1,i,1),1)          &
                              +                            &
              ddot(nphy(2),v(1,i,2),1,v(1,i,2),1)
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         first = v(1,i,1) * v(1,i,1)                       &
                          +                                &
                 v(1,i,2) * v(1,i,2)
         last  = v(nphy(2),i,1) * v(nphy(2),i,1)           &
                                +                          &
                 v(nphy(2),i,2) * v(nphy(2),i,2)
         f_2(i) = f_2(i) - .5d0 * ( first + last)
      END DO
      f_2(:) = del * f_2(:)
  END IF
  DO i=1,nphy(1)
     nrm = nrm + f_2(i)
  END DO
  IF( typke == 'fd' ) then
      nrm = nrm - .5d0 * ( f_2(1) + f_2(nphy(1)) )
      nrm = del * nrm
  END IF
  write(iout,1) nrm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_2d
!deck check_norm_3d.f
!***begin prologue     check_norm_3d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_3d
  SUBROUTINE check_norm_3d(v,f_3,v_3)
  USE dvrprop_global
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(2),nphy(1))             :: f_3
  REAL*8, DIMENSION(nphy(1))                     :: v_3
  REAL*8                                         :: ddot
  REAL*8                                         :: nrm, first, last
  INTEGER                                        :: i, j, k
!
  nrm = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        f_3(j,i) =  ddot(nphy(3),v(1,j,i,1),1,v(1,j,i,1),1)     &
                                       +                        &
                    ddot(nphy(3),v(1,j,i,2),1,v(1,j,i,2),1)
     END DO
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         DO j=1,nphy(2)
            first = v(1,j,i,1) * v(1,j,i,1)                     &
                               +                                &
                    v(1,j,i,2) * v(1,j,i,2)
            last  = v(nphy(3),j,i,1) * v(nphy(3),j,i,1)         &
                                    +                           &
                    v(nphy(3),j,i,2) * v(nphy(3),j,i,2)
            f_3(j,i) = f_3(j,i) - .5d0 * ( first + last)
         END DO
      END DO
      f_3(:,:) = del * f_3(:,:)
  END IF
  v_3(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        v_3(i) = v_3(i) + f_3(j,i)
     END DO
  END DO
  IF( typke == 'fd' ) then
      v_3(:) = v_3(:) - .5d0 * ( f_3(1,:)                       &
                                       +                        &
                                 f_3(nphy(2),:) )
      v_3(:) = del * v_3(:)
  END IF
  DO i =1,nphy(1)
     nrm = nrm + v_3(i)
  END DO
  IF( typke == 'fd' ) then
      nrm = nrm -.5d0 * ( v_3(1) + v_3(nphy(1)) )
      nrm = del * nrm
  END IF
  write(iout,1) nrm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_3d

END MODULE normalize

