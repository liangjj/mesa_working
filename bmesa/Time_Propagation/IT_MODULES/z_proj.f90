!deck z_proj
!***begin prologue     z_proj
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            tabulate gaussian wavepacket
!***                   on grid
!***references
!***routines called
!***end prologue       z_proj
  SUBROUTINE z_proj(q,wt,alpha,sigma,x_0,z,n)
  USE io
  USE dvrprop_global_it,       ONLY    : title, prnton
  USE dvr_shared,              ONLY    : typke, system            
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8                                 :: alpha, sigma, x_0
  REAL*8, DIMENSION(n)                   :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z = EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
                        / ( 2.d0 * sigma * sigma))
  ELSE
     z = q * EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) )
  END IF
  z = sqrt(wt) * z 
END SUBROUTINE z_proj
