!deck print_psi_3d
!***begin prologue     print_psi_3d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_3d
!
  SUBROUTINE print_psi_3d(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                        :: i, j, k
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i = 1,nphy(1)
     DO j = 1, nphy(2)
        DO k = 1, nphy(3)
           write(iout,3) grid(1)%pt(i), grid(2)%pt(j), grid(3)%pt(k),  &
                         wave_function(k,j,i,1), wave_function(k,j,i,2)
        END DO
     END DO
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'     y      ',5x,'     z     ','         &
                                         Real Psi  ',5x,' Imag Psi ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_3d
