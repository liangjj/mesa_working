*deck ppsgf
      function ppsgf (x, iz, c, a, bh)
c***begin prologue  ppsgf
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (ppsgf-s)
c***author  (unknown)
c***see also  blktri
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  ppsgf
      dimension       a(*)       ,c(*)       ,bh(*)
c***first executable statement  ppsgf
      sum = 0.
      do 101 j=1,iz
         sum = sum-1./(x-bh(j))**2
  101 continue
      ppsgf = sum
      return
      end
