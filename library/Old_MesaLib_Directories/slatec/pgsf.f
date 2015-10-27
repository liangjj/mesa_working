*deck pgsf
      function pgsf (x, iz, c, a, bh)
c***begin prologue  pgsf
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      single precision (pgsf-s)
c***author  (unknown)
c***see also  cblktr
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  pgsf
      dimension       a(*)       ,c(*)       ,bh(*)
c***first executable statement  pgsf
      fsg = 1.
      hsg = 1.
      do 101 j=1,iz
         dd = 1./(x-bh(j))
         fsg = fsg*a(j)*dd
         hsg = hsg*c(j)*dd
  101 continue
      if (mod(iz,2)) 103,102,103
  102 pgsf = 1.-fsg-hsg
      return
  103 pgsf = 1.+fsg+hsg
      return
      end
