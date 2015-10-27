!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               MODULE eigenvalues
!**begin prologue     eigenvalues module
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       eigenvalues
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               INTERFACE eigen_solve
             MODULE PROCEDURE eigen_solve_1d_d, eigen_solve_2d_d,   &
                              eigen_solve_3d_d
                           END INTERFACE eigen_solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               CONTAINS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck eigen_solve_1d_d.f
!***begin prologue     eigen_solve_1d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           eigenvalues
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate eigenvalues from TD wavefunction 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       eigen_solve_1d_d
  SUBROUTINE eigen_solve_1d_d(vector_0,vector)
  USE arnoldi_global
  USE dvr_shared  
  USE normalize
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: vector_0, vector
  REAL*8                                 :: ddot, first, last
  REAL*8                                 :: numerator, denominator
  REAL*8                                 :: lambda, norm
!
  call iosys('read real solution from bec',nphy(1),vector_0,0,' ')
  numerator   = ddot(nphy(1),vector_0,1,vector,1)    
  denominator = ddot(nphy(1),vector,1,vector,1)       
  lambda = numerator / denominator
  lambda = log(lambda)/deltat
  write(iout,1) lambda
1 FORMAT(/,5x,'Energy from Exponential            = ',e20.12)
END SUBROUTINE eigen_solve_1d_d
!deck eigen_solve_2d_d.f
!***begin prologue     eigen_solve_2d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           eigenvalues
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate eigenvalues from TD wavefunction 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       eigen_solve_2d_d
  SUBROUTINE eigen_solve_2d_d(vector_0,vector,g_2)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                  :: i 
  REAL*8, DIMENSION(nphy(2),nphy(1))       :: vector_0, vector
  REAL*8, DIMENSION(nphy(1))               :: g_2
  REAL*8                                   :: ddot, first, last
  REAL*8                                   :: numerator, denominator
  REAL*8                                   :: lambda
!
  call iosys('read real solution from bec',   &
              nphy(2)*nphy(1),vector_0,0,' ')
  DO i=1,nphy(1)
     g_2(i) = ddot(nphy(2),vector_0(1,i),1,vector(1,i),1)
  END DO
  numerator = 0.d0
  DO i=1,nphy(1)
     numerator = numerator + g_2(i)
  END DO
  DO i=1,nphy(1)
     g_2(i) = ddot(nphy(2),vector(1,i),1,vector(1,i),1)
  END DO
  denominator = 0.d0
  DO i=1,nphy(1)
     denominator = denominator + g_2(i)
  END DO
  lambda = numerator / denominator
  lambda = log(lambda)/deltat
  write(iout,1) lambda
1 FORMAT(/,5x,'Energy from Exponential            = ',e20.12)
END SUBROUTINE eigen_solve_2d_d
!deck eigen_solve_3d_d.f
!***begin prologue     eigen_solve_3d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           eignevalues
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate eigenvalues from TD wavefunction 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       eigen_solve_3d_d
  SUBROUTINE eigen_solve_3d_d(vector_0,vector,g_2,g_3)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  INTEGER                                      :: i, j 
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))   :: vector_0, vector
  REAL*8, DIMENSION(nphy(1))                   :: g_2
  REAL*8, DIMENSION(nphy(2),nphy(1))           :: g_3
  REAL*8                                       :: ddot, first 
  REAL*8                                       :: last
  REAL*8                                       :: numerator
  REAL*8                                       :: denominator
  REAL*8                                       :: lambda
!
  call iosys('read real wavefunction from bec',             &
              nphy(3)*nphy(2)*nphy(1),vector_0,0,' ')
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_3(j,i) = ddot(nphy(3),vector_0(1,j,i),1,vector(1,j,i),1)
     END DO
  END DO
  g_2(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(i) = g_2(i) + g_3(j,i)
     END DO
  END DO
  numerator = 0.d0
  DO i=1,nphy(1)
     numerator = numerator + g_2(i)
  END DO
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_3(j,i) = ddot(nphy(3),vector(1,j,i),1,vector(1,j,i),1)
     END DO
  END DO
  g_2(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(i) = g_2(i) + g_3(j,i)
     END DO
  END DO
  denominator = 0.d0
  DO i=1,nphy(1)
     denominator = denominator + g_2(i)
  END DO
  lambda = numerator / denominator
  lambda = log(lambda)/deltat
  write(iout,1) lambda
1 FORMAT(/,5x,'Energy from Exponential            = ',e20.12)
END SUBROUTINE eigen_solve_3d_d
END MODULE eigenvalues

