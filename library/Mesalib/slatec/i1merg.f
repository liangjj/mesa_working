*deck i1merg
      subroutine i1merg (icos, i1, m1, i2, m2, i3)
c***begin prologue  i1merg
c***subsidiary
c***purpose  merge two strings of ascending integers.
c***library   slatec
c***type      integer (s1merg-s, d1merg-d, c1merg-c, i1merg-i)
c***author  boland, w. robert, (lanl)
c           clemens, reginald, (plk)
c***description
c
c   this subroutine merges two ascending strings of integers in the
c   array icos.  the first string is of length m1 and starts at
c   icos(i1+1).  the second string is of length m2 and starts at
c   icos(i2+1).  the merged string goes into icos(i3+1).
c
c***routines called  icopy
c***revision history  (yymmdd)
c   920202  date written
c***end prologue  i1merg
      integer i1, i2, i3, m1, m2
      real icos(*)
c
      integer j1, j2, j3
c
c***first executable statement  i1merg
      if (m1.eq.0 .and. m2.eq.0) return
c
      if (m1.eq.0 .and. m2.ne.0) then
         call icopy (m2, icos(i2+1), 1, icos(i3+1), 1)
         return
      endif
c
      if (m1.ne.0 .and. m2.eq.0) then
         call icopy (m1, icos(i1+1), 1, icos(i3+1), 1)
         return
      endif
c
      j1 = 1
      j2 = 1
      j3 = 1
c
   10 if (icos(i1+j1) .le. icos(i2+j2)) then
         icos(i3+j3) = icos(i1+j1)
         j1 = j1+1
         if (j1 .gt. m1) then
            call icopy (m2-j2+1, icos(i2+j2), 1, icos(i3+j3+1), 1)
            return
         endif
      else
         icos(i3+j3) = icos(i2+j2)
         j2 = j2+1
         if (j2 .gt. m2) then
            call icopy (m1-j1+1, icos(i1+j1), 1, icos(i3+j3+1), 1)
            return
         endif
      endif
      j3 = j3+1
      go to 10
      end
