!deck abcoef
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-01-21  Time: 10:31:47

SUBROUTINE abcoef(a,b,a0,b0,c0,d0,eta,l,n,energy,level)
!***begin prologue     abcoef
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            expansion coefficients for coulomb functions
!***                   at small rho.
!***description
!***references         nbs mathematical handbook

!***routines called

!***end prologue       abcoef


REAL*8, INTENT(OUT)                      :: a(l+1:l+1+n)
REAL*8, INTENT(OUT)                      :: b(-l:-l+n)
REAL*8, INTENT(OUT)                      :: a0(0:n)
REAL*8, INTENT(OUT)                      :: b0(0:n)
REAL*8, INTENT(OUT)                      :: c0(0:n)
REAL*8, INTENT(OUT)                      :: d0(0:n)
REAL*8, INTENT(IN)                       :: eta
INTEGER, INTENT(IN)                      :: l
INTEGER, INTENT(IN)                      :: n
CHARACTER (LEN=*), INTENT(IN)            :: energy
INTEGER, INTENT(IN)                      :: level
IMPLICIT INTEGER (a-z)
REAL*8  twoeta, plfun, plf, zero, one, two, sqrt2
REAL*8 gamma, psi, argdum



DATA zero, one, two, sqrt2 /0.d0,1.d0,2.d0, 1.41421356237309505D+00/
COMMON /io/ inp, iout

IF (energy == 'positive') THEN
  twoeta=two*eta
!**********************************************************************c
!             generate the a coefficients needed for the regular       c
!                            series solution                           c
!**********************************************************************c
  a(l+1)=one
  a(l+2)=eta/(l+1)
  DO  i=l+3,l+1+n
    a(i)=twoeta*a(i-1)-a(i-2)
    a(i)=a(i)/((i+l)*(i-l-1))
  END DO
!**********************************************************************c
!              now for the b coefficients needed for the irregular     c
!              solution, things are more complicated since             c
!              they depend on the previous generation of the a's       c
!**********************************************************************c
  plf=plfun(eta,l)
  b(-l)=one
  DO  i=-l+1,-l+n
    IF (i == l+1) THEN
      b(i)=0.d0
    ELSE
      b(i)=twoeta*b(i-1)-b(i-2)-(i+i-1)*plf*a(i)
      b(i)=b(i)/((i-l-1)*(i+l))
    END IF
  END DO
  IF (level >= 2) THEN
    WRITE(iout,*)
    WRITE(iout,*) '          the a and b expansion '// 'coefficients'
    WRITE(iout,*)
    WRITE(iout,*) '      ia          a           ib'// '          b'
    DO  i=0,n
      ia=l+1+i
      ib=-l+i
      WRITE(iout,1) ia,a(ia),ib,b(ib)
    END DO
  END IF
ELSE IF (energy == 'negative') THEN
!**********************************************************************c
!             coefficients for series solution for negative energy     c
!             regular and irregular coulomb functions at small         c
!                               distances.                             c
!          * the indexing is less fancy here and follows the paper     c
!                     of Henry and Roundtree in CPC                    c
!          * there are errors in that paper which will ne noted        c
!          * the coefficients defined here are identical to the        c
!            paper cited and the errors have been fixed in the         c
!            definition of the small rho and large rho behavior        c
!                           of the functions                           c
!**********************************************************************c
  argdum=eta+l+1
  a0(0)=gamma(argdum)
  argdum=l+l+2
  a0(0)=a0(0)/gamma(argdum)
  DO  i=1,n
    a0(i)=(i+eta+l)*a0(i-1)/(i*(i+l+l+1))
  END DO
  argdum=1.d0
  d0(0)=psi(argdum)
  argdum=l+l+2
  d0(0)=d0(0)+psi(argdum)
  argdum=eta+l+1
  d0(0)=d0(0)-psi(argdum)
  DO  i=1,n
    d0(i)=d0(i-1)+1.d0/i+1.d0/(i+l+l+1.d0)-1.d0/(i+eta+l)
  END DO
  DO  i=0,n
    b0(i)=a0(i)*d0(i)
  END DO
  argdum=l+l+1.d0
  c0(0)=gamma(eta-l)*gamma(argdum)
  IF (l > 0) THEN
    DO  i=1,l+l
      c0(i)=(i+eta-l-1)*c0(i-1)/(i*(l+l+1-i))
    END DO
  END IF
  IF (level >= 2) THEN
    WRITE(iout,*)
    WRITE(iout,*) '          the a,b,c and d expansion '// 'coefficients'
    WRITE(iout,*)
    WRITE(iout,*) '      i          a                b'
    DO  i=0,n
      WRITE(iout,2) i,a0(i),b0(i)
    END DO
    WRITE(iout,*)
    WRITE(iout,*) '      i          c                d'
    twoel=l+l
    DO  i=0,n
      IF (i <= twoel) THEN
        WRITE(iout,2) i,c0(i),d0(i)
      ELSE
        WRITE(iout,3) i,d0(i)
      END IF
    END DO
  END IF
END IF
1 FORMAT(5X,i3,3X,e15.8,3X,i3,3X,e15.8)
2 FORMAT(5X,i3,3X,e15.8,3X,e15.8)
3 FORMAT(5X,i3,21X,e15.8)
RETURN
END SUBROUTINE abcoef















