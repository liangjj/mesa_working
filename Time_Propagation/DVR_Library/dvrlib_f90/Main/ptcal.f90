!deck ptcal
!***begin prologue     ptcal
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            determine the size of the gid with all constraints
!***description
!***references
!***routines called    
!***end prologue       ptcal
  SUBROUTINE PTCAL(nphy,nglobal,typwt)
  USE dvr_global, ONLY  : iout, nreg, npt, bcl, bcr  
  IMPLICIT NONE
  INTEGER nphy, nglobal, i
  CHARACTER (LEN = * )  typwt
!
  IF(typwt == 'hermite') THEN
     nphy = npt(1)
     nglobal=nphy
  ELSE IF(typwt == 'laguerre') THEN
     nphy = npt(1)
     nglobal=nphy
     IF(bcl == 0) THEN
        nphy=nphy-1
     END IF
  ELSE IF(typwt == 'fourier') THEN
     nphy=npt(1)
     nglobal=nphy
  ELSE
     nphy=0
     DO  i=1,nreg
!
!        number of internal functions
! 
          nphy = nphy + npt(i) - 2
     END DO
! 
!     add one bridge function between each interval.
!
     nphy = nphy + nreg - 1
! 
!     add the extreme left and extreme right points
!     we have not yet dropped any functions at the endpoints.
!
     nphy=nphy + 2
     nglobal=nphy
     IF(bcl == 0) THEN
        nphy=nphy-1
     END IF
     IF(bcr == 0) THEN
        nphy=nphy-1
     END IF
  END IF
END SUBROUTINE ptcal
