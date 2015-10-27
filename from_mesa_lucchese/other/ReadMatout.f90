PROGRAM ReadMatout
   IMPLICIT NONE
   INTEGER, PARAMETER :: XR = 8
   INTEGER, PARAMETER :: NKeep = 10
   INTEGER :: ios, ios2
   CHARACTER (LEN = 256) :: LineIn, Line0, Line1
   INTEGER, DIMENSION(5) :: ColNum
   INTEGER :: ColUse, iCol
   INTEGER, DIMENSION(NKeep, 5) :: Row
   REAL (KIND = XR), DIMENSION(NKeep, 5) :: Val
   INTEGER :: RowIn
   REAL (KIND = XR), DIMENSION(5) :: ValIn, SumSq
   INTEGER :: j, k

   Line0 = ' '
   Line1 = ' '

   DO
      READ (UNIT = 5, FMT = "(a)", IOSTAT = ios) LineIn
      IF (ios /= 0) THEN
         EXIT
      END IF
      IF (LineIn(1:23) /= '     ------------------') THEN
         Line0 = Line1
         Line1 = LineIn
         CYCLE
      END IF
      READ (UNIT = Line1, FMT = "(5x, 10(5x, i7))") ColNum
      IF (ColNum(1) == 1) THEN
         WRITE (UNIT = 6, FMT = "('Results from ', a)") TRIM(ADJUSTL(Line0))
      END IF
      DO ColUse = 1, 5
         IF (ColNum(ColUse) == 0) THEN
            EXIT
         END IF
      END DO
      ColUse = ColUse-1
      IF (ColUse == 0) THEN
         STOP 'Zero Columns'
      END IF
      Val = 0.0_XR
      Row = 0
      SumSq = 0.0_XR
      DO
         ios2 = 0
         READ (UNIT = 5, FMT = "(a)", IOSTAT = ios) LineIn
         IF (ios /= 0) THEN
            EXIT
         END IF
         READ (UNIT = LineIn, FMT = "(i7,2x,5f12.7)", IOSTAT = ios2) RowIn, ValIn(1:ColUse)
         IF (ios2 /= 0 .OR. RowIn == 0) THEN
            EXIT
         END IF

! **     WRITE (UNIT = 6, FMT = "('Row In ', i7, 2x, 5f13.7)")  RowIn, ValIn(1:ColUse)

         SumSq(1:ColUse) = SumSq(1:ColUse) + ValIn(1:colUse)**2

! ** see if these are one of the top 10 terms
         DO iCol = 1, ColUse
            DO j = 1, NKeep
               IF (ABS(ValIn(iCol)) > ABS(Val(j, iCol))) THEN
! ** shift off and add in
                  DO k = NKeep-1, j, -1
                     Row(k+1, iCol) = Row(k, iCol)
                     Val(k+1, iCol) = Val(k, iCol)
                  END DO
                  Row(j, iCol) = RowIn
                  Val(j, iCol) = ValIn(iCol)
                  EXIT
               END IF
            END DO
         END DO
      END DO

      DO iCol = 1, ColUse
         WRITE (UNIT = 6, FMT = "('Most important terms for column ', i7, '  Sum of Squares ', f12.7)")&
              & ColNum(iCol), SumSq(iCol)
         DO j = 1, NKeep
            IF (Row(j, iCol) /= 0) THEN
               WRITE (UNIT = 6, FMT = "(i7, 2x, f12.7)") Row(j, iCol), Val(j, iCol)
            END IF
         END DO
      END DO
      IF (ios /= 0) THEN
! ** End of file
         EXIT
      END IF
      Line0=' '
      Line1=LineIn
   END DO
END PROGRAM ReadMatout
