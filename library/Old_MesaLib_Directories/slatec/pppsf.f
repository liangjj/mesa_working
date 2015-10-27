*deck pppsf
      function pppsf (x, iz, c, a, bh)
c***begin prologue  pppsf
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      single precision (pppsf-s)
c***author  (unknown)
c***see also  cblktr
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  pppsf
      dimension       a(*)       ,c(*)       ,bh(*)
c***first executable statement  pppsf
      sum = 0.
      do 101 j=1,iz
         sum = sum+1./(x-bh(j))
  101 continue
      pppsf = sum
      return
      end
