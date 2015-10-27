!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE z_project
                              USE io
                              USE dvrprop_global,          ONLY    : title, prnton, eye
                              USE dvr_shared,              ONLY    : typke, system            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE z_proj
              MODULE PROCEDURE z_proj_d,                                &
                               z_proj_z                       
                          END INTERFACE z_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck z_proj_d
!***begin prologue     z_proj_d
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            tabulate gaussian wavepacket
!***                   on grid
!***references
!***routines called
!***end prologue       z_proj_d
  SUBROUTINE z_proj_d(q,wt,alpha,sigma,x_0,n_p,beta,z,n)
  IMPLICIT NONE
  INTEGER                                :: n, n_p
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8                                 :: alpha, sigma, x_0, beta
  REAL*8, DIMENSION(n)                   :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
                        / ( 2.d0 * sigma * sigma))
  ELSE
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) )
  END IF
  z = sqrt(wt) * z 
END SUBROUTINE z_proj_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck z_proj_z
!***begin prologue     z_proj_z
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            tabulate gaussian wavepacket
!***                   on grid
!***references
!***routines called
!***end prologue       z_proj_z
  SUBROUTINE z_proj_z(q,wt,alpha,sigma,x_0,n_p,beta,z,n)
  IMPLICIT NONE
  INTEGER                                :: n, n_p
  REAL*8, DIMENSION(n)                   :: q, wt
  REAL*8                                 :: alpha, sigma, x_0, beta
  COMPLEX*16, DIMENSION(n)               :: z
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 )  &
                        / ( 2.d0 * sigma * sigma))
     z = z * exp( - eye * beta * ( q - x_0) )
  ELSE
     z = q**n_p * EXP( - alpha * ( q - x_0 ) * ( q - x_0 ) )
     z = z * exp( - eye * beta * ( q - x_0) )
  END IF
  z = sqrt(wt) * z 
END SUBROUTINE z_proj_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE z_project
