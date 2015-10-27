!deck dl_fun
  REAL*8 FUNCTION dl_fun()
!***begin prologue     dl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        coulomb d  function
!                               l
!                      this function is trivially related to c
!                                                             l
!***references         NBS handbook

!***routines called    cl_fun(code)
!***end prologue       dl_fun
  IMPLICIT NONE
  dl_fun = one/( (l_val + l_val + 1)*cl_fun() )
  END FUNCTION dl_fun


