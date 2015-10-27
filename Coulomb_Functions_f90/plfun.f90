!deck plfun
  REAL*8 FUNCTION pl_fun()
!***begin prologue     pl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        coulomb p  function
!                               l
!***references         NBS handbook
!***routines called
!***end prologue       pl_fun
  IMPLICIT NONE
  REAL*8                   :: eta_sq
  INTEGER                  :: i
  INTEGER                  :: two_el
  INTEGER                  :: two_el_1
  INTEGER                  :: two_el_1_m_1
  pl_fun=two*eta_in
  IF (l_val > 0) THEN
      eta_sq = eta_in*eta_in
      pl_fun = pl_fun*(one+eta_sq)/three
      IF (l_val > 1) THEN
          DO i = 2, l_val
             two_el = i + i
             two_el_1 = two_el + 1
             two_el_m_1 = two_el-1
             pl_fun = pl_fun*four*(l_val*l_val + eta_sq) / (two_el_1*two_el*two_el*two_el_m_1)
          END DO
      END IF
  END IF
  END FUNCTION pl_fun


