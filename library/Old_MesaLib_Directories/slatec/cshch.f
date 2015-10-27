*deck cshch
      subroutine cshch (z, csh, cch)
c***begin prologue  cshch
c***subsidiary
c***purpose  subsidiary to cbesh and cbesk
c***library   slatec
c***type      all (cshch-a, zshch-a)
c***author  amos, d. e., (snl)
c***description
c
c     cshch computes the complex hyperbolic functions csh=sinh(x+i*y)
c     and cch=cosh(x+i*y), where i**2=-1.
c
c***see also  cbesh, cbesk
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cshch
      complex cch, csh, z
      real cchi, cchr, ch, cn, cshi, cshr, sh, sn, x, y
c***first executable statement  cshch
      x = real(z)
      y = aimag(z)
      sh = sinh(x)
      ch = cosh(x)
      sn = sin(y)
      cn = cos(y)
      cshr = sh*cn
      cshi = ch*sn
      csh = cmplx(cshr,cshi)
      cchr = ch*cn
      cchi = sh*sn
      cch = cmplx(cchr,cchi)
      return
      end
