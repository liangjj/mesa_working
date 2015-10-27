*deck cdiv
      subroutine cdiv (ar, ai, br, bi, cr, ci)
c***begin prologue  cdiv
c***subsidiary
c***purpose  compute the complex quotient of two complex numbers.
c***library   slatec
c***type      complex (cdiv-c)
c***author  (unknown)
c***description
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
c***see also  eisdoc
c***routines called  (none)
c***revision history  (yymmdd)
c   811101  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  cdiv
      real ar,ai,br,bi,cr,ci
c
      real s,ars,ais,brs,bis
c***first executable statement  cdiv
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
