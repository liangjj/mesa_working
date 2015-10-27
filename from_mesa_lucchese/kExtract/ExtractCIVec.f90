PROGRAM ExtractCIVec
   USE Extract_mod
   IMPLICIT NONE
   INTEGER (KIND = 4), PARAMETER :: UnitOutContract = 10, UnitIn = 11, UnitOut = 12
   INTEGER, PARAMETER :: NKeep = 10
   INTEGER, PARAMETER :: UnitScr = 13
   INTEGER (KIND = 4) :: ios, ios2
   CHARACTER (LEN = 256) :: LineIn, Line0, Line1
   INTEGER, DIMENSION(5) :: ColNum
   INTEGER :: ColUse, iCol
   INTEGER, DIMENSION(NKeep, 5) :: Row
   REAL (KIND = XR), DIMENSION(NKeep, 5) :: Val
   INTEGER :: RowIn
   REAL (KIND = XR), DIMENSION(5) :: ValIn, SumSq
   INTEGER :: j, k
   INTEGER :: ConfigCount, iConfig, MaxOrb
   CHARACTER (LEN = 8), PARAMETER :: Syms = "abcdefgh" 
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: ConfigOcc
   CHARACTER (LEN = 256), ALLOCATABLE, DIMENSION(:) :: ConfigLab
! ** ConfigOcc(is, iorb, isym, iconf)
! **  is = 1 for alpha and 2 for beta
! **  iorb = orbital number
! **  isym = symmetry type
! **  iconf = configuration number
   INTEGER :: iConf
   INTEGER :: OrbIn, OrbIn2, iOrb, OrbIn0, iOcc, iSym
   CHARACTER (LEN = 1) :: SymIn, Occ, SymIn0, SymIn2, Occ0, Occ2
   CHARACTER (LEN = 4) :: OccSym
   INTEGER, DIMENSION(2, 4) :: OccVal
   CHARACTER (LEN = 7) :: ThisTerm

   OPEN (UNIT = UnitIn, FILE = "mesa.out", POSITION = 'REWIND')
   OPEN (UNIT = UnitOut, FILE = "outExtractCIVec", POSITION = 'REWIND')

   WRITE (UNIT = UnitOut, FMT ="('Program ExtractCIVec - get CI vector description from mesa.out')")

   OccSym = " /\x"
   OccVal(:,1) = (/ 0,0 /)
   OccVal(:,2) = (/ 1,0 /)
   OccVal(:,3) = (/ 0,1 /)
   OccVal(:,4) = (/ 1,1 /)

   Line0 = ' '
   Line1 = ' '
   ConfigCount = 0
   OPEN (UNIT = UnitScr, STATUS = "SCRATCH", FORM = "FORMATTED")
   MainLoop: DO
      READ (UNIT = UnitIn, FMT = "(a)", IOSTAT = ios) LineIn
      IF (ios /= 0) THEN
         EXIT
      END IF
      IF (LineIn == '     configuration list') THEN
! ** Read in configuration list
         ConfigCount = 0
         MaxOrb = 0
         REWIND (UNIT = UnitScr)
         READ (UNIT = UnitIn, FMT = "(a)") LineIn
         READ (UNIT = UnitIn, FMT = "(a)") LineIn
         ConfigLoop:  DO
            READ (UNIT = LineIn, FMT = "(i6)") iConfig
            ConfigCount = ConfigCount+1
            WRITE (UNIT = UnitScr, FMT = "(a)") TRIM(LineIn)
            DO
               READ (UNIT = UnitIn, FMT = "(a)", IOSTAT = ios) LineIn
               IF (ios /= 0 .OR. LineIn == ' ') THEN
                  WRITE (UNIT = UnitScr, FMT = "()") 
                  EXIT ConfigLoop
               END IF
               IF (ConfigCount == 1) THEN
                  IF (LineIn(7:7) == ' ') THEN
                     CYCLE ConfigLoop
                  END IF
                  IF (LineIn(7:7) /= '-') THEN
                     READ (UNIT = LineIn, FMT = "(i6)") OrbIn
                  ELSE
                     READ (UNIT = LineIn, FMT = "(7x, i4)") OrbIn
                  END IF
                  MaxOrb = MAX(MaxOrb, OrbIn)

               ELSE IF (VERIFY(LineIn(9:9),Syms) /= 0) THEN
! ** this is the beginning of a new congfiguration
                  CYCLE ConfigLoop
               END IF
               WRITE (UNIT = UnitScr, FMT = "(a)") TRIM(LineIn)
            END DO
         END DO ConfigLoop
         IF (ConfigCount > 0) THEN
            WRITE (UNIT = UnitOut, FMT = "('Configurations found ', i10)") ConfigCount
            IF (ALLOCATED(ConfigOcc)) THEN
               DEALLOCATE (ConfigOcc, ConfigLab)
            END IF
            ALLOCATE (ConfigOcc(2, MaxOrb, 8, ConfigCount))
            ALLOCATE (ConfigLab(ConfigCount))
            ConfigOcc = 0
            REWIND (UNIT = UnitScr)
            READ (UNIT = UnitScr, FMT = "(a)") LineIn
            iConf = 0
            ConfigLoop2:  DO
               READ (UNIT = LineIn, FMT = "(i6)") iConfig
               iConf = iConf+1
               IF (iConf > 1) THEN
                  ConfigOcc(:,:,:,iConf) = ConfigOcc(:,:,:,1)
               END IF
               IF (iConf /= iConfig) STOP 'Error Bad configuration count'
               DO
                  READ (UNIT = UnitScr, FMT = "(a)", IOSTAT = ios) LineIn
                  IF (ios /= 0 .OR. LineIn == ' ') THEN
                     EXIT ConfigLoop2
                  END IF
                  IF (iConf == 1) THEN
                     IF (LineIn(7:7) == ' ') THEN
                        CYCLE ConfigLoop2
                     END IF
                     IF (LineIn(7:7) /= '-') THEN
                        READ (UNIT = LineIn, FMT = "(i6, a1, 2x, a1)") OrbIn, SymIn, Occ
                        ConfigOcc(:,OrbIn,INDEX(Syms, SymIn),iConf) = OccVal(:,INDEX(OccSym, Occ))
                     ELSE
                        READ (UNIT = LineIn, FMT = "(i6, 1x, i4, a1, 2x, a1)") OrbIn, OrbIn2, SymIn, Occ
                        DO iOrb = OrbIn, OrbIn2
                           ConfigOcc(:,iOrb,INDEX(Syms, SymIn),iConf) = OccVal(:,INDEX(OccSym, Occ))
                        END DO
                     END IF
                  ELSE IF (VERIFY(LineIn(9:9),Syms) /= 0) THEN
! ** this is the beginning of a new congfiguration
                     CYCLE ConfigLoop2
                  ELSE
                     READ (UNIT = LineIn, FMT = "(5x, i3, a1, 2x, a1, 3x, i3, a1, 2x, a1, 2x, i3, a1, 2x, a1)")&
                          & OrbIn0, SymIn0, Occ0, OrbIn, SymIn, Occ, OrbIn2, SymIn2, Occ2
                     ConfigOcc(:,OrbIn0,INDEX(Syms, SymIn0),iConf) = ConfigOcc(:,OrbIn0,INDEX(Syms, SymIn0),iConf)&
                          & -OccVal(:,INDEX(OccSym, Occ0))
                     ConfigOcc(:,OrbIn,INDEX(Syms, SymIn),iConf) = ConfigOcc(:,OrbIn,INDEX(Syms, SymIn),iConf)&
                          & +OccVal(:,INDEX(OccSym, Occ))
                     IF (OrbIn2 > 0) THEN
                        ConfigOcc(:,OrbIn2,INDEX(Syms, SymIn2),iConf) = ConfigOcc(:,OrbIn2,INDEX(Syms, SymIn2),iConf)&
                             & +OccVal(:,INDEX(OccSym, Occ2))
                     END IF
                  END IF
                  IF (ANY(ConfigOcc(:,:,:,iConf) > 1) .OR. ANY(ConfigOcc(:,:,:,iConf) < 0)) THEN
                     WRITE (UNIT = UnitOut, FMT = "('Error Bad configuration ', i8)") iConf
                     STOP
                  END IF
               END DO
            END DO ConfigLoop2
            DO iConf = 1, ConfigCount
               LineIn = ' '
               DO iSym = 1, 8
                  DO iOrb = 1, MaxOrb
                     IF (SUM(ConfigOcc(:,iOrb,iSym,iConf)) > 0) THEN
                        DO iOcc = 1, 4
                           IF (ALL(OccVal(:,iOcc) == ConfigOcc(:,iOrb,iSym,iConf))) THEN
                              EXIT
                           END IF
                        END DO
                        IF (iOcc > 4) STOP 'Error Bad iOcc'
                        WRITE (UNIT = ThisTerm, FMT = "(i5, a1, a1)") iOrb, Syms(iSym:iSym), OccSym(iOcc:iOcc)
                        LineIn = TRIM(LineIn)//' '//ADJUSTL(ThisTerm)
                     END IF
                  END DO
               END DO
               ConfigLab(iConf) = LineIn
            END DO
            WRITE (UNIT = UnitOut, FMT = "('Reference Configuration ', a)") TRIM(ConfigLab(1))
         END IF
         CYCLE MainLoop
      END IF
      IF (LineIn(1:23) /= '     ------------------') THEN
         Line0 = Line1
         Line1 = LineIn
         CYCLE
      END IF
      READ (UNIT = Line1, FMT = "(5x, 10(5x, i7))") ColNum
      IF (ColNum(1) == 1) THEN
         WRITE (UNIT = UnitOut, FMT = "('Results from ', a)") TRIM(ADJUSTL(Line0))
      END IF
      DO ColUse = 1, 5
         IF (ColNum(ColUse) == 0) THEN
            EXIT
         END IF
      END DO
      ColUse = ColUse-1
      IF (ColUse == 0) THEN
         STOP 'Error Zero Columns'
      END IF
      Val = 0.0_XR
      Row = 0
      SumSq = 0.0_XR
      DO
         ios2 = 0
         READ (UNIT = UnitIn, FMT = "(a)", IOSTAT = ios) LineIn
         IF (ios /= 0) THEN
            EXIT
         END IF
         READ (UNIT = LineIn, FMT = "(i7,2x,5f12.7)", IOSTAT = ios2) RowIn, ValIn(1:ColUse)
         IF (ios2 /= 0 .OR. RowIn == 0) THEN
            EXIT
         END IF

! **     WRITE (UNIT = UnitOut, FMT = "('Row In ', i7, 2x, 5f13.7)")  RowIn, ValIn(1:ColUse)

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
         WRITE (UNIT = UnitOut, FMT = "('Most important terms for column ', i7, '  Sum of Squares ', f12.7)")&
              & ColNum(iCol), SumSq(iCol)
         DO j = 1, NKeep
            IF (Row(j, iCol) /= 0) THEN
               IF (ConfigCount > 0 .AND. Row(j, iCol) <= ConfigCount) THEN
                  WRITE (UNIT = UnitOut, FMT = "(i7, 2x, f12.7, 2x, a)") Row(j, iCol), Val(j, iCol), &
                       & TRIM(ConfigLab(Row(j, iCol)))
               ELSE
                  WRITE (UNIT = UnitOut, FMT = "(i7, 2x, f12.7)") Row(j, iCol), Val(j, iCol)
               END IF
            END IF
         END DO
      END DO
      IF (ios /= 0) THEN
! ** End of file
         EXIT
      END IF
      Line0=' '
      Line1=LineIn
   END DO MainLoop
END PROGRAM ExtractCIVec
