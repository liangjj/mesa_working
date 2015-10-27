!deck c0_sq_fn
  REAL*8 FUNCTION c0_sq_fn()
!***begin prologue     c0_sq_fn
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!                                         2
!***description        coulomb function c
!                                        0
!***references         NBS handbook

!***routines called

!***end prologue       c0_sq_fn
  IMPLICIT NONE
  REAL*8 two_pi_eta
  two_pi_eta = two * pi * eta
  c0_sq_fn = two_pi_eta/(EXP(two_pi_eta)-one)
  END FUNCTION c0_sq_fn
