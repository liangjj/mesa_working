!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           MODULE non_linear_potential_module
                           USE dvrprop_global
                           USE dvr_shared
                           USE dvr_global
!**begin prologue     non_linear_potential
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       non_linear_potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                 INTERFACE v_nl
                 MODULE PROCEDURE v_nl_d,                               &
                                  v_nl_z
                             END INTERFACE v_nl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_d.f
!**begin prologue     v_nl_d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_d
  SUBROUTINE v_nl_d(wave_function)
  IMPLICIT NONE   
  REAL*8, DIMENSION(n3d)             :: wave_function
  IF (spdim == 1 ) THEN
      CALL v_nl_1d_d(wave_function)
  ELSE IF (spdim == 2 ) THEN
      CALL v_nl_2d_d(wave_function)
  ELSE IF ( spdim == 3 ) THEN
      CALL v_nl_3d_d(wave_function)
  END IF
END SUBROUTINE v_nl_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_z.f
!**begin prologue     v_nl_z
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_z
  SUBROUTINE v_nl_z(wave_function)
  IMPLICIT NONE   
  COMPLEX*16, DIMENSION(n3d)             :: wave_function
  INTEGER                                :: i
  IF (spdim == 1 ) THEN
      CALL v_nl_1d_z(wave_function)
  ELSE IF (spdim == 2 ) THEN
      CALL v_nl_2d_z(wave_function)
  ELSE IF ( spdim == 3 ) THEN
      CALL v_nl_3d_z(wave_function)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_1d_d.f
!**begin prologue     v_nl_1d_d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_1d_d
  SUBROUTINE v_nl_1d_d(wave_function)
  IMPLICIT NONE   
  REAL*8, DIMENSION(nphy(1))             :: wave_function
  INTEGER                                :: i
  DO i=1,nphy(1)
     v_tot(i) = v_tot(i) + nl_coef * wave_function(i) * wave_function(i)    &
                                                      *                     &
                                       grid(1)%f(i,i) * grid(1)%f(i,i)
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_2d_d.f
!**begin prologue     v_nl_2d_d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_2d_d
  SUBROUTINE v_nl_2d_d(wave_function)
  IMPLICIT NONE   
  REAL*8, DIMENSION(nphy(2),nphy(1))             :: wave_function
  REAL*8                                         :: fi, fj
  INTEGER                                        :: i, j, count
  count=0
  DO j=1,nphy(2)
     fj = grid(2)%f(j,j) * grid(2)%f(j,j)
     DO i=1,nphy(1)
        fi = grid(1)%f(i,i) * grid(1)%f(i,i)
        count=count+1
        v_tot(count) = v_tot(count)    +                                    &
                       nl_coef * wave_function(j,i) * wave_function(j,i)    &
                                                    *                       &
                                                 fj * fi
     END DO
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_3d_d.f
!**begin prologue     v_nl_3d_d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_3d_d
  SUBROUTINE v_nl_3d_d(wave_function)
  IMPLICIT NONE   
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))             :: wave_function
  REAL*8                                                 :: fi, fj, fk
  INTEGER                                                :: i, j, k, count
  count=0
  DO k=1,nphy(3)
     fk = grid(3)%f(k,k) * grid(3)%f(k,k)
     DO j=1,nphy(2)
        fj = grid(2)%f(j,j) * grid(2)%f(j,j)
        DO i=1,nphy(1)
           fi = grid(1)%f(i,i) * grid(1)%f(i,i)
           count=count+1
           v_tot(count) = v_tot(count)    +                                        &
                          nl_coef * wave_function(k,j,i) * wave_function(k,j,i)    &
                                                         *                         &
                                                      fk * fj * fi
        END DO
     END DO
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_1d_z
!**begin prologue     v_nl_1d_z
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_1d_z
  SUBROUTINE v_nl_1d_z(wave_function)
  IMPLICIT NONE   
  COMPLEX*16, DIMENSION(nphy(1))         :: wave_function
  INTEGER                                :: i
  DO i=1,nphy(1)
     v_tot(i) = v_tot(i) + nl_coef * wave_function(i)              &
                                   * conjg(wave_function(i))       &
                                   *                               &
                                     grid(1)%f(i,i) * grid(1)%f(i,i)
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_1d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_2d_z
!**begin prologue     v_nl_2d_z
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_2d_z
  SUBROUTINE v_nl_2d_z(wave_function)
  IMPLICIT NONE   
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))         :: wave_function
  REAL*8                                         :: fi, fj
  INTEGER                                        :: i, j, count
  count=0
  DO j=1,nphy(2)
     fj = grid(2)%f(j,j) * grid(2)%f(j,j)
     DO i=1,nphy(1)
        fi = grid(1)%f(i,i) * grid(1)%f(i,i)
        count=count+1
        v_tot(count) = v_tot(count)    +                      &
                       nl_coef * wave_function(j,i)           &
                               * conjg(wave_function(j,i))    &
                               *                              &
                            fj * fi
     END DO
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_2d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck v_nl_3d_z
!**begin prologue     v_nl_3d_z
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_3d_z
  SUBROUTINE v_nl_3d_z(wave_function)
  IMPLICIT NONE   
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))         :: wave_function
  REAL*8                                                 :: fi, fj, fk
  INTEGER                                                :: i, j, k, count
  count=0
  DO k=1,nphy(3)
     fk = grid(3)%f(k,k) * grid(3)%f(k,k)
     DO j=1,nphy(2)
        fj = grid(2)%f(j,j) * grid(2)%f(j,j)
        DO i=1,nphy(1)
           fi = grid(1)%f(i,i) * grid(1)%f(i,i)
           count=count+1
           v_tot(count) = v_tot(count)    +                       &
                          nl_coef * wave_function(k,j,i)          &
                                  * conjg(wave_function(k,j,i))   &
                                  *                               &
                            fk * fj * fi
        END DO
     END DO
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,n3d,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE non_linear_potential_module
