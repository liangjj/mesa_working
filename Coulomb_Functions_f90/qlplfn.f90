!deck ql_pl_fn

!***begin prologue     ql_pl_fn
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        q /p  function for construction of irregular positive
!                       l  l energy coulomb functions
!***references          nbs mathematical handbook (abromowitz and stegun)

!***routines called    cpsi(complex psi function:clams)
!***end prologue       ql_pl_fn
  REAL*8 FUNCTION ql_pl_fn()
  IMPLICIT NONE
  COMPLEX*16                             :: cpsi
  COMPLEX*16                             :: eta_fac
  REAL*8                                 :: eta_sq
  REAL*8                                 :: sum_1
  REAL*8                                 :: sum_2
  eta_fac = one + eta_in*eye_in
  IF (l_val == 0) THEN
      ql_pl_fn = -one + REAL( cpsi(etafac) ) + two*eulerc
  ELSE
      eta_sq = eta_in*eta_in
      sum_1=zero
      DO i = 1 , l_val
         sum_1 = sum_1 + i/( i * i + eta_sq)
      END DO
      sum_2 = zero
      DO i = 1_val ,l_val + l_val + 1
         sum_2 = sum_2 + one / i
      END DO
      ql_pl_fn = sum_1-sum_2 + REAL( cpsi(eta_fac) ) +two*eulerc +  &
                               rlfun()/plfun()
END IF
RETURN
END FUNCTION qlplfn














