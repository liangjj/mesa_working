      subroutine cdiv(ar,ai,br,bi,cr,ci)
c***begin prologue  cdiv
c***refer to  eisdoc
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c***routines called  (none)
c***end prologue  cdiv
      implicit real *8 (a-h,o-z)
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
