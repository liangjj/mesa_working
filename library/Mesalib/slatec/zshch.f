*deck zshch
      subroutine zshch (zr, zi, cshr, cshi, cchr, cchi)
c***begin prologue  zshch
c***subsidiary
c***purpose  subsidiary to zbesh and zbesk
c***library   slatec
c***type      all (cshch-a, zshch-a)
c***author  amos, d. e., (snl)
c***description
c
c     zshch computes the complex hyperbolic functions csh=sinh(x+i*y)
c     and cch=cosh(x+i*y), where i**2=-1.
c
c***see also  zbesh, zbesk
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zshch
c
      double precision cchi, cchr, ch, cn, cshi, cshr, sh, sn, zi, zr
c***first executable statement  zshch
      sh = sinh(zr)
      ch = cosh(zr)
      sn = sin(zi)
      cn = cos(zi)
      cshr = sh*cn
      cshi = ch*sn
      cchr = ch*cn
      cchi = sh*sn
      return
      end
