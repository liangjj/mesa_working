*deck zexp
      subroutine zexp (ar, ai, br, bi)
c***begin prologue  zexp
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi, zbesj, zbesk, zbesy, zairy and
c            zbiry
c***library   slatec
c***type      all (zexp-a)
c***author  amos, d. e., (snl)
c***description
c
c     double precision complex exponential function b=exp(a)
c
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbesy, zbiry
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zexp
      double precision ar, ai, br, bi, zm, ca, cb
c***first executable statement  zexp
      zm = exp(ar)
      ca = zm*cos(ai)
      cb = zm*sin(ai)
      br = ca
      bi = cb
      return
      end
