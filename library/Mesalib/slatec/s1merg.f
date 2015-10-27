*deck s1merg
      subroutine s1merg (tcos, i1, m1, i2, m2, i3)
c***begin prologue  s1merg
c***subsidiary
c***purpose  merge two strings of ascending real numbers.
c***library   slatec
c***type      single precision (s1merg-s, d1merg-d, c1merg-c, i1merg-i)
c***author  (unknown)
c***description
c
c   this subroutine merges two ascending strings of numbers in the
c   array tcos.  the first string is of length m1 and starts at
c   tcos(i1+1).  the second string is of length m2 and starts at
c   tcos(i2+1).  the merged string goes into tcos(i3+1).
c
c***see also  genbun
c***routines called  scopy
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   901120  modified to use if-then-else.  previous spaghetti code did
c           not compile correctly with optimization on the ibm rs6000.
c           (rwc)
c   920130  code name changed from merge to s1merg.  (wrb)
c***end prologue  s1merg
      integer i1, i2, i3, m1, m2
      real tcos(*)
c
      integer j1, j2, j3
c
c***first executable statement  s1merg
      if (m1.eq.0 .and. m2.eq.0) return
c
      if (m1.eq.0 .and. m2.ne.0) then
         call scopy (m2, tcos(i2+1), 1, tcos(i3+1), 1)
         return
      endif
c
      if (m1.ne.0 .and. m2.eq.0) then
         call scopy (m1, tcos(i1+1), 1, tcos(i3+1), 1)
         return
      endif
c
      j1 = 1
      j2 = 1
      j3 = 1
c
   10 if (tcos(i1+j1) .le. tcos(i2+j2)) then
         tcos(i3+j3) = tcos(i1+j1)
         j1 = j1+1
         if (j1 .gt. m1) then
            call scopy (m2-j2+1, tcos(i2+j2), 1, tcos(i3+j3+1), 1)
            return
         endif
      else
         tcos(i3+j3) = tcos(i2+j2)
         j2 = j2+1
         if (j2 .gt. m2) then
            call scopy (m1-j1+1, tcos(i1+j1), 1, tcos(i3+j3+1), 1)
            return
         endif
      endif
      j3 = j3+1
      go to 10
      end
