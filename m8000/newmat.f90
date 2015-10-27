!deck @(#)newmat.f 1.1 9/7/91
!***begin prologue     newmat
!***date written       910304   (yymmdd)
!***revision date               (yymmdd)
!***keywords           m6010, link 6010, kohn data
!***author             schneider, barry (lanl)
!***source             m6010
!***purpose            re-arrange matrices
!***                   for kohn codes.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       newmat

  SUBROUTINE newmat(matin,matout,scr,type)
  USE mesa_global
  IMPLICIT NONE
  CHARACTER (LEN=*)                       :: type
  INTEGER                                 :: i, j, ii, jj
  REAL*8, DIMENSION(nbf,nbf)              :: matin                       
  REAL*8, DIMENSION(nbf,nmotot)           :: matout                       
  REAL*8, DIMENSION(nmotot,nmotot)        :: scr                       
  IF (TYPE == 'ao-mo') THEN
      matout(1:nbf,1:nmotot) = 0.d0
      DO  i=1,nmotot
          ii=list(i)
          DO  j=1,nbf
              matout(j,i)=matin(j,ii)
          END DO
      END DO
      call copy(matout,matin,nbf*nmotot)
  ELSE IF(TYPE == 'mo-mo') THEN
      scr(1:nmotot,1:nmotot) = 0.d0
      DO  i=1,nmotot
          ii=list(i)
          DO  j=1,nmotot
              jj=list(j)
              scr(i,j)=matin(ii,jj)
          END DO
      END DO
      call copy(scr,matin,nmo*nmo)
  ELSE
    CALL lnkerr('error in newmat')
  END IF
END SUBROUTINE newmat
