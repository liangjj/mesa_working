PROGRAM MoreThan72
   IMPLICIT NONE
   CHARACTER (LEN = 132) :: Line
   INTEGER :: Count
   INTEGER :: ios

   Count = 0
   DO
      READ (UNIT = 5, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         EXIT
      END IF
      Count = Count+1
      IF (Line(1:1) /= 'c' .AND. Line(1:1) /= 'C') THEN
         IF (Line(73:) /= ' ') THEN
            WRITE (UNIT = 6, FMT = "('Long line ', i10, ' -', a)") Count, TRIM(Line)
         END IF
      END IF
   END DO

END PROGRAM MoreThan72
   
