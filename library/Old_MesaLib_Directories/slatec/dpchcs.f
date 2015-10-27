*deck dpchcs
      subroutine dpchcs (switch, n, h, slope, d, incfd, ierr)
c***begin prologue  dpchcs
c***subsidiary
c***purpose  adjusts derivative values for dpchic
c***library   slatec (pchip)
c***type      double precision (pchcs-s, dpchcs-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c         dpchcs:  dpchic monotonicity switch derivative setter.
c
c     called by  dpchic  to adjust the values of d in the vicinity of a
c     switch in direction of monotonicity, to produce a more "visually
c     pleasing" curve than that given by  dpchim .
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  n, ierr
c        double precision  switch, h(n), slope(n), d(incfd,n)
c
c        call  dpchcs (switch, n, h, slope, d, incfd, ierr)
c
c   parameters:
c
c     switch -- (input) indicates the amount of control desired over
c           local excursions from data.
c
c     n -- (input) number of data points.  (assumes n.gt.2 .)
c
c     h -- (input) real*8 array of interval lengths.
c     slope -- (input) real*8 array of data slopes.
c           if the data are (x(i),y(i)), i=1(1)n, then these inputs are:
c                  h(i) =  x(i+1)-x(i),
c              slope(i) = (y(i+1)-y(i))/h(i),  i=1(1)n-1.
c
c     d -- (input) real*8 array of derivative values at the data points,
c           as determined by dpchci.
c          (output) derivatives in the vicinity of switches in direction
c           of monotonicity may be adjusted to produce a more "visually
c           pleasing" curve.
c           the value corresponding to x(i) is stored in
c                d(1+(i-1)*incfd),  i=1(1)n.
c           no other entries in d are changed.
c
c     incfd -- (input) increment between successive values in d.
c           this argument is provided primarily for 2-d applications.
c
c     ierr -- (output) error flag.  should be zero.
c           if negative, trouble in dpchsw.  (should never happen.)
c
c    -------
c    warning:  this routine does no validity-checking of arguments.
c    -------
c
c  fortran intrinsics used:  abs, max, min.
c
c***see also  dpchic
c***routines called  dpchst, dpchsw
c***revision history  (yymmdd)
c   820218  date written
c   820617  redesigned to (1) fix  problem with lack of continuity
c           approaching a flat-topped peak (2) be cleaner and
c           easier to verify.
c           eliminated subroutines pchsa and pchsx in the process.
c   820622  1. limited fact to not exceed one, so computed d is a
c             convex combination of dpchci value and dpchsd value.
c           2. changed fudge from 1 to 4 (based on experiments).
c   820623  moved pchsd to an inline function (eliminating mswtyp).
c   820805  converted to slatec library version.
c   870707  corrected conversion to double precision.
c   870813  minor cosmetic changes.
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  modified spacing in computation of dfloc.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated author section in prologue.  (wrb)
c   930503  improved purpose.  (fnf)
c***end prologue  dpchcs
c
c  programming notes:
c     1. the function  dpchst(arg1,arg2)  is assumed to return zero if
c        either argument is zero, +1 if they are of the same sign, and
c        -1 if they are of opposite sign.
c**end
c
c  declare arguments.
c
      integer  n, incfd, ierr
      double precision  switch, h(*), slope(*), d(incfd,*)
c
c  declare local variables.
c
      integer  i, indx, k, nless1
      double precision  del(3), dext, dfloc, dfmx, fact, fudge, one,
     *      slmax, wtave(2), zero
      save zero, one, fudge
      double precision  dpchst
c
c  define inline function for weighted average of slopes.
c
      double precision  dpchsd, s1, s2, h1, h2
      dpchsd(s1,s2,h1,h2) = (h2/(h1+h2))*s1 + (h1/(h1+h2))*s2
c
c  initialize.
c
      data  zero /0.d0/,  one/1.d0/
      data  fudge /4.d0/
c***first executable statement  dpchcs
      ierr = 0
      nless1 = n - 1
c
c  loop over segments.
c
      do 900  i = 2, nless1
         if ( dpchst(slope(i-1),slope(i)) )  100, 300, 900
c             --------------------------
c
  100    continue
c
c....... slope switches monotonicity at i-th point .....................
c
c           do not change d if 'up-down-up'.
            if (i .gt. 2)  then
               if ( dpchst(slope(i-2),slope(i)) .gt. zero)  go to 900
c                   --------------------------
            endif
            if (i .lt. nless1)  then
               if ( dpchst(slope(i+1),slope(i-1)) .gt. zero)  go to 900
c                   ----------------------------
            endif
c
c   ....... compute provisional value for d(1,i).
c
            dext = dpchsd (slope(i-1), slope(i), h(i-1), h(i))
c
c   ....... determine which interval contains the extremum.
c
            if ( dpchst(dext, slope(i-1)) )  200, 900, 250
c                -----------------------
c
  200       continue
c              dext and slope(i-1) have opposite signs --
c                        extremum is in (x(i-1),x(i)).
               k = i-1
c              set up to compute new values for d(1,i-1) and d(1,i).
               wtave(2) = dext
               if (k .gt. 1)
     *            wtave(1) = dpchsd (slope(k-1), slope(k), h(k-1), h(k))
               go to 400
c
  250       continue
c              dext and slope(i) have opposite signs --
c                        extremum is in (x(i),x(i+1)).
               k = i
c              set up to compute new values for d(1,i) and d(1,i+1).
               wtave(1) = dext
               if (k .lt. nless1)
     *            wtave(2) = dpchsd (slope(k), slope(k+1), h(k), h(k+1))
               go to 400
c
  300    continue
c
c....... at least one of slope(i-1) and slope(i) is zero --
c                     check for flat-topped peak .......................
c
            if (i .eq. nless1)  go to 900
            if ( dpchst(slope(i-1), slope(i+1)) .ge. zero)  go to 900
c                -----------------------------
c
c           we have flat-topped peak on (x(i),x(i+1)).
            k = i
c           set up to compute new values for d(1,i) and d(1,i+1).
            wtave(1) = dpchsd (slope(k-1), slope(k), h(k-1), h(k))
            wtave(2) = dpchsd (slope(k), slope(k+1), h(k), h(k+1))
c
  400    continue
c
c....... at this point we have determined that there will be an extremum
c        on (x(k),x(k+1)), where k=i or i-1, and have set array wtave--
c           wtave(1) is a weighted average of slope(k-1) and slope(k),
c                    if k.gt.1
c           wtave(2) is a weighted average of slope(k) and slope(k+1),
c                    if k.lt.n-1
c
         slmax = abs(slope(k))
         if (k .gt. 1)    slmax = max( slmax, abs(slope(k-1)) )
         if (k.lt.nless1) slmax = max( slmax, abs(slope(k+1)) )
c
         if (k .gt. 1)  del(1) = slope(k-1) / slmax
         del(2) = slope(k) / slmax
         if (k.lt.nless1)  del(3) = slope(k+1) / slmax
c
         if ((k.gt.1) .and. (k.lt.nless1))  then
c           normal case -- extremum is not in a boundary interval.
            fact = fudge* abs(del(3)*(del(1)-del(2))*(wtave(2)/slmax))
            d(1,k) = d(1,k) + min(fact,one)*(wtave(1) - d(1,k))
            fact = fudge* abs(del(1)*(del(3)-del(2))*(wtave(1)/slmax))
            d(1,k+1) = d(1,k+1) + min(fact,one)*(wtave(2) - d(1,k+1))
         else
c           special case k=1 (which can occur only if i=2) or
c                        k=nless1 (which can occur only if i=nless1).
            fact = fudge* abs(del(2))
            d(1,i) = min(fact,one) * wtave(i-k+1)
c              note that i-k+1 = 1 if k=i  (=nless1),
c                        i-k+1 = 2 if k=i-1(=1).
         endif
c
c
c....... adjust if necessary to limit excursions from data.
c
         if (switch .le. zero)  go to 900
c
         dfloc = h(k)*abs(slope(k))
         if (k .gt. 1)    dfloc = max( dfloc, h(k-1)*abs(slope(k-1)) )
         if (k.lt.nless1) dfloc = max( dfloc, h(k+1)*abs(slope(k+1)) )
         dfmx = switch*dfloc
         indx = i-k+1
c        indx = 1 if k=i, 2 if k=i-1.
c        ---------------------------------------------------------------
         call dpchsw(dfmx, indx, d(1,k), d(1,k+1), h(k), slope(k), ierr)
c        ---------------------------------------------------------------
         if (ierr .ne. 0)  return
c
c....... end of segment loop.
c
  900 continue
c
      return
c------------- last line of dpchcs follows -----------------------------
      end
