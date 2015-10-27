*deck c1merg
      subroutine c1merg (tcos, i1, m1, i2, m2, i3)
c***begin prologue  c1merg
c***subsidiary
c***purpose  merge two strings of complex numbers.  each string is
c            ascending by the real part.
c***library   slatec
c***type      complex (s1merg-s, d1merg-d, c1merg-c, i1merg-i)
c***author  (unknown)
c***description
c
c   this subroutine merges two ascending strings of numbers in the
c   array tcos.  the first string is of length m1 and starts at
c   tcos(i1+1).  the second string is of length m2 and starts at
c   tcos(i2+1).  the merged string goes into tcos(i3+1).  the ordering
c   is on the real part.
c
c***see also  cmgnbn
c***routines called  ccopy
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   910408  modified to use if-then-else.  make it look like merge
c           which was modified earlier due to compiler problems on
c           the ibm rs6000.  (rwc)
c   920130  code name changed from cmpmrg to c1merg.  (wrb)
c***end prologue  c1merg
      integer i1, i2, i3, m1, m2
      complex tcos(*)
c
      integer j1, j2, j3
c
c***first executable statement  c1merg
      if (m1.eq.0 .and. m2.eq.0) return
c
      if (m1.eq.0 .and. m2.ne.0) then
         call ccopy (m2, tcos(i2+1), 1, tcos(i3+1), 1)
         return
      endif
c
      if (m1.ne.0 .and. m2.eq.0) then
         call ccopy (m1, tcos(i1+1), 1, tcos(i3+1), 1)
         return
      endif
c
      j1 = 1
      j2 = 1
      j3 = 1
c
   10 if (real(tcos(j1+i1)) .le. real(tcos(i2+j2))) then
         tcos(i3+j3) = tcos(i1+j1)
         j1 = j1+1
         if (j1 .gt. m1) then
            call ccopy (m2-j2+1, tcos(i2+j2), 1, tcos(i3+j3+1), 1)
            return
         endif
      else
         tcos(i3+j3) = tcos(i2+j2)
         j2 = j2+1
         if (j2 .gt. m2) then
            call ccopy (m1-j1+1, tcos(i1+j1), 1, tcos(i3+j3+1), 1)
            return
         endif
      endif
      j3 = j3+1
      go to 10
      end
