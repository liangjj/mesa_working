*deck dintrv
      subroutine dintrv (xt, lxt, x, ilo, ileft, mflag)
c***begin prologue  dintrv
c***purpose  compute the largest integer ileft in 1 .le. ileft .le. lxt
c            such that xt(ileft) .le. x where xt(*) is a subdivision of
c            the x interval.
c***library   slatec
c***category  e3, k6
c***type      double precision (intrv-s, dintrv-d)
c***keywords  b-spline, data fitting, interpolation, splines
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract    **** a double precision routine ****
c         dintrv is the interv routine of the reference.
c
c         dintrv computes the largest integer ileft in 1 .le. ileft .le.
c         lxt such that xt(ileft) .le. x where xt(*) is a subdivision of
c         the x interval.  precisely,
c
c                      x .lt. xt(1)                1         -1
c         if  xt(i) .le. x .lt. xt(i+1)  then  ileft=i  , mflag=0
c           xt(lxt) .le. x                         lxt        1,
c
c         that is, when multiplicities are present in the break point
c         to the left of x, the largest index is taken for ileft.
c
c     description of arguments
c
c         input      xt,x are double precision
c          xt      - xt is a knot or break point vector of length lxt
c          lxt     - length of the xt vector
c          x       - argument
c          ilo     - an initialization parameter which must be set
c                    to 1 the first time the spline array xt is
c                    processed by dintrv.
c
c         output
c          ilo     - ilo contains information for efficient process-
c                    ing after the initial call and ilo must not be
c                    changed by the user.  distinct splines require
c                    distinct ilo parameters.
c          ileft   - largest integer satisfying xt(ileft) .le. x
c          mflag   - signals when x lies out of bounds
c
c     error conditions
c         none
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  (none)
c***revision history  (yymmdd)
c   800901  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dintrv
c
      integer ihi, ileft, ilo, istep, lxt, mflag, middle
      double precision x, xt
      dimension xt(*)
c***first executable statement  dintrv
      ihi = ilo + 1
      if (ihi.lt.lxt) go to 10
      if (x.ge.xt(lxt)) go to 110
      if (lxt.le.1) go to 90
      ilo = lxt - 1
      ihi = lxt
c
   10 if (x.ge.xt(ihi)) go to 40
      if (x.ge.xt(ilo)) go to 100
c
c *** now x .lt. xt(ihi) . find lower bound
      istep = 1
   20 ihi = ilo
      ilo = ihi - istep
      if (ilo.le.1) go to 30
      if (x.ge.xt(ilo)) go to 70
      istep = istep*2
      go to 20
   30 ilo = 1
      if (x.lt.xt(1)) go to 90
      go to 70
c *** now x .ge. xt(ilo) . find upper bound
   40 istep = 1
   50 ilo = ihi
      ihi = ilo + istep
      if (ihi.ge.lxt) go to 60
      if (x.lt.xt(ihi)) go to 70
      istep = istep*2
      go to 50
   60 if (x.ge.xt(lxt)) go to 110
      ihi = lxt
c
c *** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval
   70 middle = (ilo+ihi)/2
      if (middle.eq.ilo) go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1
      if (x.lt.xt(middle)) go to 80
      ilo = middle
      go to 70
   80 ihi = middle
      go to 70
c *** set output and return
   90 mflag = -1
      ileft = 1
      return
  100 mflag = 0
      ileft = ilo
      return
  110 mflag = 1
      ileft = lxt
      return
      end
