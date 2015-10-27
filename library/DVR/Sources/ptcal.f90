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
  SUBROUTINE PTCAL(nphy,nglobal)
  USE dvr_global, ONLY  : iout, nreg, npt, bcl, bcr  
  IMPLICIT NONE
  INTEGER                    :: nphy
  INTEGER                    :: nglobal
  INTEGER                    :: i
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
  Write(iout,1) nglobal, nphy
1 Format(/,20x,'Number of Global Grid Points   = ',i5,/,20x,
               'Number of Physical Grid Points = ',i5)
END SUBROUTINE ptcal
