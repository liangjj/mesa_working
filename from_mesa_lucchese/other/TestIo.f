      PROGRAM iotest
c **      INTEGER, parameter :: irecl = 32
c **       INTEGER, PARAMETER :: irecl = 131584
      DIMENSION ix(100), iy(100)

      irecl=32
      DO 20 i = 1, 100
         ix(i) = i
20    continue

      DO 60 nx = 1, 40
      OPEN (unit = 10, status = 'SCRATCH', form = 'UNFORMATTED',
     X     access = 'DIRECT', recl = irecl)

      WRITE (6, '('' irecl ='', i8)') irecl

      DO 30 i = 1, nx
         WRITE (6, '('' i ='', i8)') i
         WRITE (10, rec = i, err = 10) (ix(j), j = 1, i)
30    continue

      DO 50 i = 1, nx
         READ(10, rec = i, err = 10) (iy(j), j = 1, i)
         WRITE (6, '(''i ='', i5, /, (10i5))') i, (iy(j), j = 1, i)
         DO 40 j = 1, i
            IF (iy(j) .NE. ix(j)) THEN
               WRITE (6, '(''Bad read at i='', i10)') i
               stop
            end if
40       continue
50    continue
60    continue

      CLOSE (UNIT = 10)

      STOP

 10   CONTINUE
      WRITE (6, '('' Error occurred'')')
      STOP
      END

