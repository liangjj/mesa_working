!deck plot_wavefunction
!***begin prologue     plot_wavefunction
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       plot_wavefunction
!
  SUBROUTINE plot_wavefunction(wave_function)
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
END SUBROUTINE plot_wavefunction
