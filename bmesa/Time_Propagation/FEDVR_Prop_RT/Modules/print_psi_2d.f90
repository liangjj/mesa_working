!deck print_psi_2d
!***begin prologue     print_psi_2d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_2d
!
  SUBROUTINE print_psi_2d(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: i, j
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i = 1,nphy(1)
     DO j = 1, nphy(2)
        write(iout,3) grid(1)%pt(i), grid(2)%pt(j),                 &
                      wave_function(j,i,1), wave_function(j,i,2)
     END DO
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'     y      ',5x,'  Real Psi  ',5x,' Imag Psi ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_2d