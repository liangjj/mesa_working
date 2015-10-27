!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     MODULE spatial_wavefunction
!**begin prologue     spatial_wavefunction module
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       auto_correlation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     INTERFACE spatial_psi
         MODULE PROCEDURE spatial_psi_d,                                &
                          spatial_psi_z
                 END INTERFACE spatial_psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_d
!***begin prologue     spatial_psi_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_d
  SUBROUTINE spatial_psi_d(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                 :: wave_function, scratch_vector
!
  IF(spdim == 1 ) THEN
     CALL spatial_psi_1d_d(wave_function,scratch_vector)
  ELSE IF(spdim == 2 ) THEN
     CALL spatial_psi_2d_d(wave_function,scratch_vector)
  ELSE IF(spdim == 3 ) THEN
     CALL spatial_psi_3d_d(wave_function,scratch_vector)
  END IF
END SUBROUTINE spatial_psi_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_z
!***begin prologue     spatial_psi_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_z
  SUBROUTINE spatial_psi_z(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)       :: wave_function, scratch_vector
!
  IF(spdim == 1 ) THEN
     CALL spatial_psi_1d_z(wave_function,scratch_vector)
  ELSE IF(spdim == 2 ) THEN
     CALL spatial_psi_2d_z(wave_function,scratch_vector)
  ELSE IF(spdim == 3 ) THEN
     CALL spatial_psi_3d_z(wave_function,scratch_vector)
  END IF
END SUBROUTINE spatial_psi_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_1d_d
!***begin prologue     spatial_psi_1d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_1d_d
  SUBROUTINE spatial_psi_1d_d(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: wave_function,scratch_vector
!
  scratch_vector(:) = wave_function(:) /sqrt(grid(1)%wt(:))
END SUBROUTINE spatial_psi_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_1d_z
!***begin prologue     spatial_psi_1d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_1d_z
  SUBROUTINE spatial_psi_1d_z(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))         :: wave_function,scratch_vector
  scratch_vector(:) = wave_function(:) /sqrt(grid(1)%wt(:))
END SUBROUTINE spatial_psi_1d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_2d_d.f
!***begin prologue     spatial_psi_2d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_2d_d
  SUBROUTINE spatial_psi_2d_d(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: i, j 
  REAL*8, DIMENSION(nphy(2),nphy(1))     :: wave_function, scratch_vector
  REAL*8                                 :: sq_i, sq_j
  DO i=1,nphy(1)
       sq_i = 1.d0 / sqrt(grid(1)%wt(i))
       DO j=1,nphy(2)
          sq_j = 1.d0 / sqrt(grid(2)%wt(j))
          scratch_vector(j,i) = wave_function(j,i) * sq_i * sq_j
       END DO
  END DO
END SUBROUTINE spatial_psi_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_2d_z.f
!***begin prologue     spatial_psi_2d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_2d_z
  SUBROUTINE spatial_psi_2d_z(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: i, j 
  COMPLEX*16, DIMENSION(nphy(2),nphy(1)) :: wave_function, scratch_vector
  REAL*8                                 :: sq_i, sq_j
!
  DO i=1,nphy(1)
       sq_i = 1.d0 / sqrt(grid(1)%wt(i))
       DO j=1,nphy(2)
          sq_j = 1.d0 / sqrt(grid(2)%wt(j))
          scratch_vector(j,i) = wave_function(j,i) * sq_i * sq_j
       END DO
  END DO
END SUBROUTINE spatial_psi_2d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_3d_d.f
!***begin prologue     spatial_psi_3d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_3d_d
  SUBROUTINE spatial_psi_3d_d(wave_function,scratch_vector)
  USE arnoldi_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                        :: i, j, k 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))     :: wave_function, scratch_vector
  REAL*8                                         :: sq_i, sq_j, sq_k
  DO i=1,nphy(1)
       sq_i = 1.d0 / sqrt(grid(1)%wt(i))
       DO j=1,nphy(2)
          sq_j = 1.d0 / sqrt(grid(2)%wt(j))
          DO k=1,nphy(3)
             sq_k = 1.d0 / sqrt(grid(3)%wt(k))
             scratch_vector(k,j,i) = wave_function(k,j,i) * sq_i * sq_j * sq_k
          END DO
       END DO
  END DO
END SUBROUTINE spatial_psi_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck spatial_psi_3d_z.f
!***begin prologue     spatial_psi_3d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       spatial_psi_3d_z
  SUBROUTINE spatial_psi_3d_z(wave_function,scratch_vector)
  USE dvr_global
  USE dvr_shared  
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                        :: i, j, k 
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1)) :: wave_function, scratch_vector
  REAL*8                                         :: sq_i, sq_j, sq_k 
!
  DO i=1,nphy(1)
       sq_i = 1.d0 / sqrt(grid(1)%wt(i))
       DO j=1,nphy(2)
          sq_j = 1.d0 / sqrt(grid(2)%wt(j))
          DO k=1,nphy(3)
             sq_k = 1.d0 / sqrt(grid(3)%wt(k))
             scratch_vector(k,j,i) = wave_function(k,j,i) * sq_i * sq_j * sq_k
          END DO
       END DO
  END DO
END SUBROUTINE spatial_psi_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE spatial_wavefunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


