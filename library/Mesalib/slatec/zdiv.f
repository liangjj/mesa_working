*deck zdiv
      subroutine zdiv (ar, ai, br, bi, cr, ci)
c***begin prologue  zdiv
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi, zbesj, zbesk, zbesy, zairy and
c            zbiry
c***library   slatec
c***type      all (zdiv-a)
c***author  amos, d. e., (snl)
c***description
c
c     double precision complex divide c=a/b.
c
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbesy, zbiry
c***routines called  zabs
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zdiv
      double precision ar, ai, br, bi, cr, ci, bm, ca, cb, cc, cd
      double precision zabs
      external zabs
c***first executable statement  zdiv
      bm = 1.0d0/zabs(br,bi)
      cc = br*bm
      cd = bi*bm
      ca = (ar*cc+ai*cd)*bm
      cb = (ai*cc-ar*cd)*bm
      cr = ca
      ci = cb
      return
      end
