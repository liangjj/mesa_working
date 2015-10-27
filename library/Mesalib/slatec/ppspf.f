*deck ppspf
      function ppspf (x, iz, c, a, bh)
c***begin prologue  ppspf
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (ppspf-s)
c***author  (unknown)
c***see also  blktri
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  ppspf
      dimension       a(*)       ,c(*)       ,bh(*)
c***first executable statement  ppspf
      sum = 0.
      do 101 j=1,iz
         sum = sum+1./(x-bh(j))
  101 continue
      ppspf = sum
      return
      end
