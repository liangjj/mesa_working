!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE plot_wavefunction
!**begin prologue     plot_wavefunction
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       plot_wavefunction
!
  IMPLICIT NONE
  CONTAINS
!deck plot_all
!***begin prolgue      plot_all
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            plotting
!***
!***references
!***routines called
!***end prologue       plot_all
!
  SUBROUTINE plot_all(t)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: t, i
  CHARACTER (LEN=16)                     :: fptoc
  REAL*8                                 :: rtemp
!
!
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_wavefunction(psi)
  END IF
  IF (plot.and.t == ntreg) THEN
      IF(spdim==1) then
         CALL plot_1d(psi)
      ELSE IF(spdim==2) THEN
         CALL plot_2d(psi)
      ELSE
         CALL plot_3d(psi)
      END IF
  END IF
END SUBROUTINE plot_all
!deck print_wavefunction
!***begin prologue     print_wavefunction
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_wavefunction
!
  SUBROUTINE print_wavefunction(wave_function)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: i, j, k, count
  REAL*8, DIMENSION(n3d,2)               :: wave_function
!
  write(iout,1) title
  IF(spdim == 1 ) then
     write(iout,2)
      DO i=1,n3d
         write(iout,3) grid(1)%pt(i), wave_function(i,1), wave_function(i,2)
      END DO
  ELSE IF(spdim == 2) then
      count = 0
      write(iout,4)
      DO i = 1,nphy(1)
         DO j = 1, nphy(2)
            count = count + 1
            write(iout,5) grid(1)%pt(i), grid(2)%pt(j), wave_function(count,1), &
                                                        wave_function(count,2)
         END DO
      END DO
  ELSE IF(spdim == 3) then
      write(iout,6)
      count = 0
      DO i = 1,nphy(1)
         DO j = 1, nphy(2)
            DO k = 1, nphy(3)
               count = count + 1
               write(iout,7) grid(1)%pt(i), grid(2)%pt(j), grid(3)%pt(k),  &
                             wave_function(count,1), wave_function(count,2)
            END DO
         END DO
      END DO
  END IF
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'   Real Psi  ',5x,'Imag Psi ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6)
4 FORMAT('     x      ',5x,'     y      ',5x,'  Real Psi  ',5x,' Imag Psi ')
5 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
6 FORMAT('     x      ',5x,'     y      ',5x,'     z     ','         &
                                         Real Psi  ',5x,' Imag Psi ')
7 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_wavefunction
!deck plot_1d.f
!***begin prologue     plot_1d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_1d
  SUBROUTINE plot_1d(vector)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: i
  DO i = 1,nphy(1)
     rtemp = SQRT ( ( vector(i,1) * vector(i,1)        &
                                  +                    &
                      vector(i,2) * vector(i,2) ) )
     write(iplot(3),*) grid(1)%pt(i), vector(i,1)
     write(iplot(4),*) grid(1)%pt(i), vector(i,2)
     write(iplot(5),*) grid(1)%pt(i), rtemp
  END DO
END SUBROUTINE plot_1d
!deck plot_2d.f
!***begin prologue     plot_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_2d
  SUBROUTINE plot_2d(vector)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: i, j
  DO i = 1,nphy(2)
     DO j = 1,nphy(1)
        rtemp = SQRT ( ( vector(i,j,1) * vector(i,j,1)        &
                                       +                      &
                         vector(i,j,2) * vector(i,j,2) ) )
        write(iplot(3),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j,1)
        write(iplot(4),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j,2)
        write(iplot(5),*) grid(2)%pt(i), grid(1)%pt(j), rtemp
     END DO
  END DO
END SUBROUTINE plot_2d
!deck plot_3d.f
!***begin prologue     plot_3d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_3d
  SUBROUTINE plot_3d(vector)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vector
  REAL*8                                         :: rtemp
  INTEGER                                        :: i, j, k
  DO i = 1,nphy(3)
     DO j = 1,nphy(2)
        DO k = 1,nphy(1)
           rtemp = SQRT ( ( vector(i,j,k,1) * vector(i,j,k,1)      &
                                            +                      &
                            vector(i,j,k,2) * vector(i,j,k,2) ) )
          write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),          &
                            grid(1)%pt(k), vector(i,j,k,1)      
          write(iplot(4),*) grid(3)%pt(i), grid(2)%pt(j),          &
                            grid(3)%pt(k), vector(i,j,k,2)
          write(iplot(5),*) grid(3)%pt(i), grid(2)%pt(j),          &
                            grid(3)%pt(k), rtemp
        END DO
     END DO
  END DO
END SUBROUTINE plot_3d
END MODULE plot_wavefunction

