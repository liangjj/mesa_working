!deck zfil_308
!**begin prologue     zfil_308
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            fill array
!**references
!**routines called
!**end prologue       zfil_308
  SUBROUTINE zfil_308(z1,z2,z3)
  USE dvrprop_global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1),2)         :: z1
  COMPLEX*16, DIMENSION(nphy(2),2)         :: z2
  COMPLEX*16, DIMENSION(nphy(3),2)         :: z3
  INTEGER                                  :: i, j, k, count 

END SUBROUTINE zfil_308