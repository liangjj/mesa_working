*deck zabs
      double precision function zabs (zr, zi)
c***begin prologue  zabs
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi, zbesj, zbesk, zbesy, zairy and
c            zbiry
c***library   slatec
c***type      all (zabs-a)
c***author  amos, d. e., (snl)
c***description
c
c     zabs computes the absolute value or magnitude of a double
c     precision complex variable cmplx(zr,zi)
c
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbesy, zbiry
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zabs
      double precision zr, zi, u, v, q, s
c***first executable statement  zabs
      u = abs(zr)
      v = abs(zi)
      s = u + v
c-----------------------------------------------------------------------
c     s*1.0d0 makes an unnormalized underflow on cdc machines into a
c     true floating zero
c-----------------------------------------------------------------------
      s = s*1.0d+0
      if (s.eq.0.0d+0) go to 20
      if (u.gt.v) go to 10
      q = u/v
      zabs = v*sqrt(1.d+0+q*q)
      return
   10 q = v/u
      zabs = u*sqrt(1.d+0+q*q)
      return
   20 zabs = 0.0d+0
      return
      end
