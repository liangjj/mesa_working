!***********************************************************************
!*deck asymptotic_expansion_irregular_negative_energy_function
!***begin prologue     asymptotic_expansion_irregular_negative_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            asymntotic expansion for irregular negative
!***                   energy coulomb function at large rho.
!***
!***description        exponentially decaying coulomb function is
!***                   computed using an asymptotic expansion at
!***                   large rho.
!***references

!***routines called

!***end prologue       asymn
  SUBROUTINE asymptotic_expansion_irregular_nerative_energy_function
  IMPLICIT NONE
  REAL*8   pre, add, fac
!**********************************************************************c
!            gl = -1./sqrt(2.) * (series) * exp(-rho)                  c
!                              * (2.rho)**(-eta)                       c
!                                                                      c
!            series = sum ( n =0 to n= infinity ) e0(n)                c
!                             * (-2.*rho)**(-eta)                      c
!                                                                      c
!                * Henry incorrectly uses 2.*rho                       c
!**********************************************************************c
!**********************************************************************c
!                do the series expansion to the level tol              c
!**********************************************************************c
  pre=one
  fac=-half*rinv(pt)
  DO i = 0, size
     add=e_0(i)*pre
     IF (ABS(add) <= tol) EXIT
     gl = gl + add
     dgl = dgl + i*add
     pre=pre*fac
  END DO
  IF (level >= 3) THEN
      WRITE(iout,*) 'asymptotic expansion converged in ',i, ' terms'
  END IF
  fac = half*r_inv
  pre = invsqrt2*EXP(-r))
  add = fac**eta
  gl = -gl*pre*add
  dgl = -gl*(one+eta*r_inv)-pre*add*dgl*r_inv
END SUBROUTINE asymptotic_expansion_irregular_negative_energy_function
