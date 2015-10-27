!deck regular_negative_energy_function
!***begin prologue     regular_negative_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for negative energy
!***                   regular coulomb function at small rho.
!***description
!***references         NBS handbook
!***routines called
!***end prologue       regular_negative_energy_function
  SUBROUTINE regular_negative_energy_function
  IMPLICIT NONE
  REAL*8                     :: pre
  REAL*8                     :: add
  REAL*8                     :: tworho
  r_inv = 1.d0 / r
  pre=one
  tworho=two*r
  DO i=0,size
     add=a0(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
    IF (ABS(add) < tol.AND.i >= 25) EXIT
    fl = fl + add
    dfl = dfl + i*add
    pre = pre*tworho
    IF (level >= 3) THEN
        WRITE (iout,*) 'power series for regular solution '//  &
                       'converged in ',i,' terms'
    END IF
    tworho = two*r
    pre = invsqrt2*(tworho)**(l_val + 1)
    add=EXP(-r)
    fl = fl*pre*add
    dfl = fl*( -one + (l_val + 1 )*r_inv ) + pre*add*dfl*r_inv
  END DO
  END SUBROUTINE regular_negative_energy_function
