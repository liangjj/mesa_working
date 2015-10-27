*deck bcrh
      function bcrh (xll, xrr, iz, c, a, bh, f, sgn)
c***begin prologue  bcrh
c***subsidiary
c***purpose  subsidiary to cblktr
c***library   slatec
c***type      single precision (bcrh-s, bsrh-s)
c***author  (unknown)
c***see also  cblktr
c***routines called  (none)
c***common blocks    ccblk
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  bcrh
      dimension       a(*)       ,c(*)       ,bh(*)
      common /ccblk/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c***first executable statement  bcrh
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
  105 bcrh = .5*(xl+xr)
      return
      end
