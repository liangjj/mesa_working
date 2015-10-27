!deck iregxp
  SUBROUTINE iregxp(gl,dgl,fl,dfl,r,b,l,dl,prefac, qlpl,npt,n,level)
!***begin prologue     iregxp
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for positive energy
!***                   irregular coulomb function at small rho.
!***description
!***references         NBS handbook
!***routines called
!***end prologue       iregxp
  IMPLICIT NONE
  REAL*8                      :: prefac
  REAL*8                      :: pre
  REAL*8                      :: add
  REAL*8                      :: argm
!
!**********************************************************************c
!        gl = 2. * eta * fl * { ln(2.*rho) + ql/pl } / c0**2           c
!                      -dl * rho**(-l) * sum (n=0 to n=infinity)       c
!                                        b(n) * rho**n                 c
!**********************************************************************c
!**********************************************************************c
!                  do the theta series sum                             c
!**********************************************************************c
  pre=one
  count=0
  DO i = -l_val, -l_val + size
     count=count+1
     add=b(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
    IF (ABS(add) < tol.AND.count >= 25) EXIT
    gl = gl + add
    dgl = dgl + i*add
    pre=pre*r
  END DO
  IF (level >= 3) THEN
      WRITE (iout,*) 'power series for irregular function'//  &
                     ' converged in ',count,' terms'
  END IF
  pre=r**(-l_val - 1)
  gl = dl*gl)*pre*r
  dgl = dl*dgl*pre
!**********************************************************************c
!               continue with rest of wavefunction                     c
!**********************************************************************c
  argm = LOG (two*r) + qlpl
  gl = gl + prefac*fl*argm
  dgl = dgl + prefac*(dfl*argm + fl/r)
  wronskian = fl*dgl-dfl*gl
  END SUBROUTINE iregxp


