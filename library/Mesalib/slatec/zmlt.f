*deck zmlt
      subroutine zmlt (ar, ai, br, bi, cr, ci)
c***begin prologue  zmlt
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi, zbesj, zbesk, zbesy, zairy and
c            zbiry
c***library   slatec
c***type      all (zmlt-a)
c***author  amos, d. e., (snl)
c***description
c
c     double precision complex multiply, c=a*b.
c
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbesy, zbiry
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zmlt
      double precision ar, ai, br, bi, cr, ci, ca, cb
c***first executable statement  zmlt
      ca = ar*br - ai*bi
      cb = ar*bi + ai*br
      cr = ca
      ci = cb
      return
      end
