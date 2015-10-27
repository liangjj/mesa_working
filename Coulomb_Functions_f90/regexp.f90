!deck regexp
  SUBROUTINE regexp(fl,dfl,r)
!***begin prologue     regexp
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for positive energy
!***                   regular coulomb function at small rho.
!***description
!***references         NBS handbook
!***routines called
!***end prologue       regexp
  IMPLICIT NONE
  REAL*8                      :: fl
  REAL*8                      :: dfl
  REAL*8                      :: r
  REAL*8                      :: cl
  REAL*8                      :: pre
  REAL*8                      :: add
  pre=one
  count=0
  l_1 = l_val + int_one
  DO i = l_1 , l_1 + size
     count=count+1
     add=a(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
     IF (ABS(add) < tol.AND.count >= 25) EXIT
         fl = fl + add
         dfl = dfl + i*add
         pre = pre*r
  END DO
  IF (level >= 3) THEN
      WRITE (iout,*) 'power series converged in ',count,' terms'
  END IF
  pre = r**l_val
  fl =  cl * fl * pre * r
  dfl = cl * dfl * pre
END SUBROUTINE regexp



