!deck z_proj.f
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
  SUBROUTINE z_proj(q,wt,p,alpha,sigma,x_0,beta,z,np)
  USE arnoldi_global,       ONLY    : iout, eye, title, &
                                      typke, system, prnton
  USE dvr_shared
!  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: np
  REAL*8, DIMENSION(np)                  :: q, wt
  REAL*8, DIMENSION(np,np)               :: p
  REAL*8                                 :: alpha, sigma, x_0, beta
  COMPLEX*16, DIMENSION(np)              :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z = EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) / ( 2.d0 * sigma * sigma)  &
              - eye * beta * (q-x_0) )
  ELSE
     z = q * EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
            / ( 2.d0 * sigma * sigma) - eye * beta * (q-x_0) )
  END IF
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     z = z * wt
     do i=1,np
        z(i) = z(i)*p(i,i)
     END DO
  END IF
END SUBROUTINE z_proj

