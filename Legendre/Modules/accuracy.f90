!
MODULE accuracy
! Define floating-point working precision, idp
  IMPLICIT NONE
  INTEGER, PARAMETER        :: idp = SELECTED_REAL_KIND(15, 307)
END MODULE  accuracy
