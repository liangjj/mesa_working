*deck ppgsf
      function ppgsf (x, iz, c, a, bh)
c***begin prologue  ppgsf
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      single precision (ppgsf-s)
c***author  (unknown)
c***see also  cblktr
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  ppgsf
      dimension       a(*)       ,c(*)       ,bh(*)
c***first executable statement  ppgsf
      sum = 0.
      do 101 j=1,iz
         sum = sum-1./(x-bh(j))**2
  101 continue
      ppgsf = sum
      return
      end
