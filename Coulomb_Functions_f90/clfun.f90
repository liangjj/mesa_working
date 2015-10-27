!deck  cl_fun
!***begin prologue     cl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            compute c  function for positive energy coulomb
!***                            l                           functions.
!***description
!***references         NBS handbook
!***routines called    cgamma(complex gamma function:clams)
!***end prologue       m2705
  REAL*8 Function cl_fun()
  IMPLICIT NONE
  COMPLEX *16                          :: arg
  COMPLEX *16                          :: cgamma
  cl_fun = (two**l_val)*EXP(-pi*eta_in*half)
  arg=l_val+1_val+eye*eta_in
  cl_fun=cl_fun*ABS(cgamma(arg))/fact(2*l+1)
END FUNCTION clfun
