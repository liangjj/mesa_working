!deck irregular_negative_energy_function
!***begin prologue     irregular_negative_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for negative energy
!***                   irregular coulomb function at small rho.
!***description
!***references
!***routines called
!***end prologue       irregular_negative_energy_function
  SUBROUTINE irregular_negative_energy_function
  IMPLICIT NONE
  REAL*8                                 :: prefac
  REAL*8                                 :: tworho
  REAL*8                                 :: mtworo
  REAL*8                                 :: pre
  REAL*8                                 :: add
  REAL*8                                 :: expfac
  REAL*8                                 :: dxpfac
  REAL*8                                 :: fac
  REAL*8                                 :: ifac
  REAL*8                                 :: lgfac
  REAL*8                                 :: fac1
  REAL*8                                 :: fac2
  REAL*8                                 :: fac3
  REAL*8                                 :: csum
  REAL*8                                 :: dcsum
  REAL*8                                 :: asum
  REAL*8                                 :: dasum
  REAL*8                                 :: bsum
  REAL*8                                 :: dbsum
  INTEGER                                :: i
  INTEGER                                :: twoel1
!**********************************************************************c
!      gl = -1./sqrt(2.) * series * exp(-rho) * (2.*rho)**(-l) /       c
!                    (gamma(eta+l+1)*gamma(eta-l))                     c
!                                                                      c
!      series = sum(n=0 to n=2*l) * c0(n) * (-2.*rho)**n               c
!               - (-2.*rho)**(2*l+1) * { ln(2.*rho) *                  c
!                                        sum(n=0 to n=infinity)        c
!                                            a0(n) * (2.*rho)**n -     c
!                                        sum(n=0 to n=infinity)        c
!                                            b0(n) * (2.*rho)**n }     c
!                                                                      c
!      * Henry has two misprints in his paper. the c0-sum has a        c
!        + 2.*rho and the last sum is added instead of subtracted.     c
!        these errors were discovered by comparing with the NBS        c
!        handbook and using the Henry normalization.                   c
!**********************************************************************c
  r_inv = 1.d0 / r
  twoel1= l_val + l_val + 1
  prefac = invsqrt2 / (gamma(eta+l_val+1)*gamma(eta-l_val))
!**********************************************************************c
!        exponential prefactor and its derivative                      c
!**********************************************************************c
  expfac = prefac*EXP(-r)
  dxpfac = -expfac
!**********************************************************************c
!             (2.*rho)**(-l) and derivative                            c
!             then combine with exponential factor                     c
!**********************************************************************c
  ifac = one
  IF (l_val /= 0) THEN
    ifac = (half*r_inv))**l_val
    expfac = ifac*expfac
    dxpfac = ifac*dxpfac-l*expfac*r_inv
  END IF
!**********************************************************************c
!                prefactors in front of sums                           c
!**********************************************************************c
  tworho = two*r
  fac = (-tworho)**twoel1
  lgfac = LOG(tworho)
!**********************************************************************c
!                     do the six sums                                  c
!**********************************************************************c
  csum = zero
  dcsum = zero
  pre = one
  mtworo = -tworho
  DO i=0,l_val+l
     csum = csum+c0(i)*pre
     dcsum = dcsum+i*pre
     pre = pre*mtworo
  END DO
  dcsum = dcsum*r_inv
  pre = one
  asum = zero
  dasum = zero
  DO i=0,size
     add = a0(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
     IF (ABS(add) < tol.AND.i >= 25) EXIT
     asum = asum+add
     dasum = dasum+i*add
     pre = pre*tworho
  END DO
  imax = i
  dasum = dasum*r_inv
  pre = one
  bsum = zero
  dbsum = zero
  DO i=0,size
     add = b0(i)*pre
     IF (ABS(add) < tol.AND.i >= 25) EXIT
     bsum = bsum+add
     dbsum = dbsum+i*add
     pre = pre*tworho
  END DO
  imax=MAX(imax,i)
  IF (level >= 3) THEN
      WRITE (iout,*) 'power series for irregular function '//  &
                     'converged in ',imax,' terms'
  END IF
  dbsum = dbsum*r_inv
!**********************************************************************c
!        combine all the factors to get the irregular function         c
!                                and                                   c
!                             its   derivative                         c
!**********************************************************************c
  fac1 = asum*lgfac - bsum
  fac2 = asum*r_inv + lgfac*dasum - dbsum
  fac3 = csum-fac*fac1
  gl = expfac*fac3
  dgl = expfac*(dcsum-fac*fac2-twoel1*fac*fac1*r_inv) + +dxpfac*fac3
  dgl = - dgl
  gl = -gl
  END SUBROUTINE irregular_negative_energy_function


