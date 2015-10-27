*deck bsrh
      function bsrh (xll, xrr, iz, c, a, bh, f, sgn)
c***begin prologue  bsrh
c***subsidiary
c***purpose  subsidiary to blktri
c***library   slatec
c***type      single precision (bcrh-s, bsrh-s)
c***author  (unknown)
c***see also  blktri
c***routines called  (none)
c***common blocks    cblkt
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  bsrh
      dimension       a(*)       ,c(*)       ,bh(*)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  bsrh
      xl = xll
      xr = xrr
      dx = .5*abs(xr-xl)
  101 x = .5*(xl+xr)
      if (sgn*f(x,iz,c,a,bh)) 103,105,102
  102 xr = x
      go to 104
  103 xl = x
  104 dx = .5*dx
      if (dx-cnv) 105,105,101
  105 bsrh = .5*(xl+xr)
      return
      end
