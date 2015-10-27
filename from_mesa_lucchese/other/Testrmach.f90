PROGRAM Getr1mach
   IMPLICIT NONE
   INTEGER, PARAMETER :: XR = 8 

   REAL (KIND = XR) :: fac, newfac, newfacd
   INTEGER :: i, j, jmax, imax

   jmax = 2000
   imax = 2000

   DO i = -1024, 1024
      WRITE (UNIT = 6, FMT = "('2**', i5, '  ', e27.18)") i, 2.0_XR**i
   END DO

   fac = 1.0_XR
   DO i = 1, imax
      newfac = 2.0_XR**(-i)
      DO j = 1, jmax
         newfacd = newfac + 2.0_XR**(-(j+i))
         IF (newfac == newfacd) THEN
            EXIT
         END IF
      END DO
      IF (j > jmax) THEN
         STOP 'Bad jmax'
      END IF

      WRITE (UNIT = 6, FMT = "('i', i5, '  fac', e27.18, '  j', i8, '  dfac', e27.18)") -i, newfac, 1-j, 2.0_XR**(1-j) 
      IF (newfac == fac) EXIT
      fac = newfac
   END DO

   fac = 1.0_XR
   DO i = 1, imax
      newfac = 2.0_XR**(i)
      DO j = 1, jmax
         newfacd = newfac + 2.0_XR**(i-j)
         IF (newfac == newfacd) THEN
            EXIT
         END IF
      END DO
      IF (j > jmax) THEN
         STOP 'Bad jmax'
      END IF
      WRITE (UNIT = 6, FMT = "('i', i5, '  fac', e27.18, '  j', i8, '  dfac', e27.18)") i, newfac, 1-j, 2.0_XR**(1-j)
      IF (newfac == fac) EXIT
      fac = newfac
   END DO

   WRITE (UNIT = 6, FMT = "('log10(2)', e27.18)") LOG10(2.0_XR)

END PROGRAM Getr1mach

