!deck rl_fun
!***begin prologue     rl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        r function for construction of irregular positive
!                       l energy coulomb functions
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called
!***end prologue       rl_fun
  REAL*8 FUNCTION rl_fun()
  IMPLICIT NONE
  REAL*8                    :: denominator
  COMPLEX*16                :: com_eta
  COMPLEX*16                :: eta_fac
  COMPLEX*16                :: num
  COMPLEX*16                :: cr_l
  INTEGER                   :: i
  INTEGER                   :: n_terms
  INTEGER                   :: l_fac
  INTEGER                   :: i_test
  IF (l_val == 0) THEN
      rl_fun = zero
  ELSE
  com_eta = eta_in*eye
  n_terms = l_val+l_val+1
  l_fac = n_terms
  eta_fac = com_eta-l_val
  num = one
  cr_l = num/l_fac
  DO i= 2 , n_terms
     l_fac = l_fac-1
     num = num * two * eta_fac
     denominator = l_fac*fact(i-1)
     cr_l =cr_l + num / denominator
     eta_fac = eta_fac + one
  END DO
  rl_fun = imag(cr_l)
  i_test = l_val - 2 * (l_val/2)
  IF (i_test == 0) THEN
      rl_fun = -rl_fun
  END IF
  rl_fun = rl_fun/fact(l_val + l_val)
END IF
END FUNCTION rl_fun
