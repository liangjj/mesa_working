!deck asymptotic_expansion_regular_positive_energy_function
!***begin prologue     asymptotic_expansion_regular_positive_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            asymptotic expansion for regular or irregular
!***                   positive coulomb function at large rho.
!***
!***description        f  = g * cos ( theta ) + f * sin ( theta )
!                       l                  l                   l
!***                   the leading term of f is 1.0
!***                   for the irregular function switch f and g and replace
!***                   the plus by a minus sign.
!***                   g, f and their derivatives are defined in the NBS
!***                   mathematical handbook.
!***references

!***routines called

!***end prologue       asymptotic_expansion_regular_positive_energy_function
  SUBROUTINE asymptotic_expansion_regular_positive_energy_function
  IMPLICIT NONE
REAL*8  pi, carg, sigmal
REAL*8  sn, cn,  thetal
REAL*8 f, g, fd, gd, fold, gold, fnew, gnew
REAL*8 foldd, goldd, fdnew, gdnew
REAL*8  zero, half, one, two
COMPLEX *16 cgamma, eye
!**********************************************************************c
!        calculate the coulomb phase shift and the sin and cos         c
!        of the argument of the asymptotic thetal                      c
!**********************************************************************c
  sigmal=carg(cgamma(l_val+one+eye*eta))
  IF(level >= 1) THEN
     WRITE(iout,1) sigmal
  END IF
  thetal = r - eta*LOG(two*r)) -l_val*pi*half+sigmal
  sn = SIN(thetal)
  cn = COS(thetal)
!**********************************************************************c
!                do the series expansion to the level tol              c
!**********************************************************************c
  fold = one
  gold = zero
  foldd = zero
  goldd = one - eta*r_inv
  f = fold
  g = gold
  fd = foldd
  gd = goldd
  DO i=0, sizde
     fnew =(a_k(i)*fold-b_k(i)*gold)*r_inv
    gnew = (a_k(i)*gold+b_k(i)*fold)*r_inv
    fdnew = (a_k(i)*foldd-b_k(i)*goldd-fnew)*r_inv
    gdnew = (a_k(i)*goldd+b_k(i)*foldd-gnew)*r_inv
    IF (ABS(fnew) < tol.AND.ABS(gnew) < tol) EXIT
    f = f+fnew
    g = g+gnew
    fd = fd+fdnew
    gd = gd+gdnew
    fold = fnew
    gold = gnew
    foldd = fdnew
    goldd = gdnew
  END DO
  IF (level >= 3) THEN
      WRITE(iout,*) 'asymptotic expansion converged in ',i, ' terms'
  END IF
  fl = f*sn+g*cn
  dfl = fd*sn+gd*cn
  gl = f*cn-g*sn
  dgl = fd*cn-gd*sn
  wron = fl*dgl - dfl*gl
1 FORMAT(/,5X,'coulomb phase shift',1X,e15.8,/)
END SUBROUTINE asymptotic_expansion_regular_positive_energy_function

