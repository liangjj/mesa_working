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
  SUBROUTINE z_proj(q,wt,p,alpha,sigma,x_0,beta,z,n)
  USE io
  USE dvrprop_global,       ONLY    : title, typke, system, prnton
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8, DIMENSION(n,n)                 :: p
  REAL*8                                 :: alpha, sigma, x_0, beta
  REAL*8, DIMENSION(n,2)                 :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z(:,1) = EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
                           / ( 2.d0 * sigma * sigma))
     z(:,2) = - z(:,1) * sin( beta * (q-x_0) )
     z(:,1) =   z(:,1) * cos( beta * (q-x_0) )
  ELSE
     z(:,1) = q * EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) )
     z(:,2) = - z(:,1) * sin( beta * (q-x_0) )
     z(:,1) =   z(:,1) * cos( beta * (q-x_0) )
  END IF
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     z(:,1) = wt(:) * z(:,1) 
     z(:,2) = wt(:) * z(:,2) 
     do i=1,n
        z(i,:) = p(i,i) * z(i,:) 
     END DO
  END IF
END SUBROUTINE z_proj
