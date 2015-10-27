!deck z_fill
!**begin prologue     z_fill
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            fill array
!**references
!**routines called
!**end prologue       z_fill
  SUBROUTINE z_fill(z1,z2,z3)
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: z1
  REAL*8, DIMENSION(nphy(2),2)           :: z2
  REAL*8, DIMENSION(nphy(3),2)           :: z3
  INTEGER                                :: i, j, k, count 
  IF(spdim == 2) then
     count=0
     DO  i=1,nphy(1)
         DO  j=1,nphy(2)
               count=count+1
               psi(count,1) = z1(i,1)*z2(j,1) - z1(i,2)*z2(j,2)
               psi(count,2) = z1(i,1)*z2(j,2) + z1(i,2)*z2(j,1)
         END DO
     END DO
  ELSE IF(spdim == 3) then
     count=0
     DO  i=1,nphy(1)
         DO  j=1,nphy(2)
             Do k=1,nphy(3)
                count=count+1
                psi(count,1) = z1(i,1)*z2(j,1)*z3(k,1) - &
                                  z1(i,2)*z2(j,2)*z3(k,1) - &
                                  z1(i,2)*z2(j,1)*z3(k,2) - &
                                  z1(i,1)*z2(j,2)*z3(k,2)
                psi(count,2) = z1(i,2)*z2(j,1)*z3(k,1) + &
                                  z1(i,1)*z2(j,2)*z3(k,1) + &
                                  z1(i,1)*z2(j,1)*z3(k,2) - &
                                  z1(i,2)*z2(j,2)*z3(k,2)
             END DO  
         END DO
     END DO
  END IF
END SUBROUTINE z_fill
