!deck nrmlze.f
!***begin prologue     nrmlze
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            normalize lobatto polynomials and their
!**                    first derivatives.
!***
!***references

!***routines called
!***end prologue       nrmlze

  SUBROUTINE nrmlze(p,dp,wt,n)
  IMPLICIT NONE
  INTEGER                      :: n
  INTEGER                      :: i, j
  REAL*8, DIMENSION(n,n)       :: p
  REAL*8, DIMENSION(n,n)       :: dp
  REAL*8, DIMENSION(n)         :: wt
  REAL*8  nrm
  DO  i=1,n
    nrm=1.d0/SQRT(wt(i))
    p(:,i)=nrm*p(:,i)
    dp(:,i)=nrm*dp(:,i)
!    CALL vscale(p(1,i),p(1,i),nrm,n)
!    CALL vscale(dp(1,i),dp(1,i),nrm,n)
  END DO
  RETURN
END SUBROUTINE nrmlze


