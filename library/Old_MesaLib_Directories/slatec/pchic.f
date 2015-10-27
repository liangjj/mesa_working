*deck pchic
      subroutine pchic (ic, vc, switch, n, x, f, d, incfd, wk, nwk,
     +   ierr)
c***begin prologue  pchic
c***purpose  set derivatives needed to determine a piecewise monotone
c            piecewise cubic hermite interpolant to given data.
c            user control is available over boundary conditions and/or
c            treatment of points where monotonicity switches direction.
c***library   slatec (pchip)
c***category  e1a
c***type      single precision (pchic-s, dpchic-d)
c***keywords  cubic hermite interpolation, monotone interpolation,
c             pchip, piecewise cubic interpolation,
c             shape-preserving interpolation
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c         pchic:  piecewise cubic hermite interpolation coefficients.
c
c     sets derivatives needed to determine a piecewise monotone piece-
c     wise cubic interpolant to the data given in x and f satisfying the
c     boundary conditions specified by ic and vc.
c
c     the treatment of points where monotonicity switches direction is
c     controlled by argument switch.
c
c     to facilitate two-dimensional applications, includes an increment
c     between successive values of the f- and d-arrays.
c
c     the resulting piecewise cubic hermite function may be evaluated
c     by pchfe or pchfd.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  ic(2), n, nwk, ierr
c        real  vc(2), switch, x(n), f(incfd,n), d(incfd,n), wk(nwk)
c
c        call  pchic (ic, vc, switch, n, x, f, d, incfd, wk, nwk, ierr)
c
c   parameters:
c
c     ic -- (input) integer array of length 2 specifying desired
c           boundary conditions:
c           ic(1) = ibeg, desired condition at beginning of data.
c           ic(2) = iend, desired condition at end of data.
c
c           ibeg = 0  for the default boundary condition (the same as
c                     used by pchim).
c           if ibeg.ne.0, then its sign indicates whether the boundary
c                     derivative is to be adjusted, if necessary, to be
c                     compatible with monotonicity:
c              ibeg.gt.0  if no adjustment is to be performed.
c              ibeg.lt.0  if the derivative is to be adjusted for
c                     monotonicity.
c
c           allowable values for the magnitude of ibeg are:
c           ibeg = 1  if first derivative at x(1) is given in vc(1).
c           ibeg = 2  if second derivative at x(1) is given in vc(1).
c           ibeg = 3  to use the 3-point difference formula for d(1).
c                     (reverts to the default b.c. if n.lt.3 .)
c           ibeg = 4  to use the 4-point difference formula for d(1).
c                     (reverts to the default b.c. if n.lt.4 .)
c           ibeg = 5  to set d(1) so that the second derivative is con-
c              tinuous at x(2). (reverts to the default b.c. if n.lt.4.)
c              this option is somewhat analogous to the "not a knot"
c              boundary condition provided by pchsp.
c
c          notes (ibeg):
c           1. an error return is taken if abs(ibeg).gt.5 .
c           2. only in case  ibeg.le.0  is it guaranteed that the
c              interpolant will be monotonic in the first interval.
c              if the returned value of d(1) lies between zero and
c              3*slope(1), the interpolant will be monotonic.  this
c              is **not** checked if ibeg.gt.0 .
c           3. if ibeg.lt.0 and d(1) had to be changed to achieve mono-
c              tonicity, a warning error is returned.
c
c           iend may take on the same values as ibeg, but applied to
c           derivative at x(n).  in case iend = 1 or 2, the value is
c           given in vc(2).
c
c          notes (iend):
c           1. an error return is taken if abs(iend).gt.5 .
c           2. only in case  iend.le.0  is it guaranteed that the
c              interpolant will be monotonic in the last interval.
c              if the returned value of d(1+(n-1)*incfd) lies between
c              zero and 3*slope(n-1), the interpolant will be monotonic.
c              this is **not** checked if iend.gt.0 .
c           3. if iend.lt.0 and d(1+(n-1)*incfd) had to be changed to
c              achieve monotonicity, a warning error is returned.
c
c     vc -- (input) real array of length 2 specifying desired boundary
c           values, as indicated above.
c           vc(1) need be set only if ic(1) = 1 or 2 .
c           vc(2) need be set only if ic(2) = 1 or 2 .
c
c     switch -- (input) indicates desired treatment of points where
c           direction of monotonicity switches:
c           set switch to zero if interpolant is required to be mono-
c           tonic in each interval, regardless of monotonicity of data.
c             notes:
c              1. this will cause d to be set to zero at all switch
c                 points, thus forcing extrema there.
c              2. the result of using this option with the default boun-
c                 dary conditions will be identical to using pchim, but
c                 will generally cost more compute time.
c                 this option is provided only to facilitate comparison
c                 of different switch and/or boundary conditions.
c           set switch nonzero to use a formula based on the 3-point
c              difference formula in the vicinity of switch points.
c           if switch is positive, the interpolant on each interval
c              containing an extremum is controlled to not deviate from
c              the data by more than switch*dfloc, where dfloc is the
c              maximum of the change of f on this interval and its two
c              immediate neighbors.
c           if switch is negative, no such control is to be imposed.
c
c     n -- (input) number of data points.  (error return if n.lt.2 .)
c
c     x -- (input) real array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.
c           (error return if not.)
c
c     f -- (input) real array of dependent variable values to be inter-
c           polated.  f(1+(i-1)*incfd) is value corresponding to x(i).
c
c     d -- (output) real array of derivative values at the data points.
c           these values will determine a monotone cubic hermite func-
c           tion on each subinterval on which the data are monotonic,
c           except possibly adjacent to switches in monotonicity.
c           the value corresponding to x(i) is stored in
c                d(1+(i-1)*incfd),  i=1(1)n.
c           no other entries in d are changed.
c
c     incfd -- (input) increment between successive values in f and d.
c           this argument is provided primarily for 2-d applications.
c           (error return if  incfd.lt.1 .)
c
c     wk -- (scratch) real array of working storage.  the user may wish
c           to know that the returned values are:
c              wk(i)     = h(i)     = x(i+1) - x(i) ;
c              wk(n-1+i) = slope(i) = (f(1,i+1) - f(1,i)) / h(i)
c           for  i = 1(1)n-1.
c
c     nwk -- (input) length of work array.
c           (error return if  nwk.lt.2*(n-1) .)
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           warning errors:
c              ierr = 1  if ibeg.lt.0 and d(1) had to be adjusted for
c                        monotonicity.
c              ierr = 2  if iend.lt.0 and d(1+(n-1)*incfd) had to be
c                        adjusted for monotonicity.
c              ierr = 3  if both of the above are true.
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c              ierr = -4  if abs(ibeg).gt.5 .
c              ierr = -5  if abs(iend).gt.5 .
c              ierr = -6  if both of the above are true.
c              ierr = -7  if nwk.lt.2*(n-1) .
c             (the d-array has not been changed in any of these cases.)
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c
c***references  1. f. n. fritsch, piecewise cubic hermite interpolation
c                 package, report ucrl-87285, lawrence livermore nation-
c                 al laboratory, july 1982.  [poster presented at the
c                 siam 30th anniversary meeting, 19-23 july 1982.]
c               2. f. n. fritsch and j. butland, a method for construc-
c                 ting local monotone piecewise cubic interpolants, siam
c                 journal on scientific and statistical computing 5, 2
c                 (june 1984), pp. 300-304.
c               3. f. n. fritsch and r. e. carlson, monotone piecewise
c                 cubic interpolation, siam journal on numerical ana-
c                 lysis 17, 2 (april 1980), pp. 238-246.
c***routines called  pchce, pchci, pchcs, xermsg
c***revision history  (yymmdd)
c   820218  date written
c   820804  converted to slatec library version.
c   870813  updated reference 2.
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890703  corrected category record.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920429  revised format and order of references.  (wrb,fnf)
c***end prologue  pchic
c  programming notes:
c
c     to produce a double precision version, simply:
c        a. change pchic to dpchic wherever it occurs,
c        b. change pchce to dpchce wherever it occurs,
c        c. change pchci to dpchci wherever it occurs,
c        d. change pchcs to dpchcs wherever it occurs,
c        e. change the real declarations to double precision, and
c        f. change the constant  zero  to double precision.
c
c  declare arguments.
c
      integer  ic(2), n, incfd, nwk, ierr
      real  vc(2), switch, x(*), f(incfd,*), d(incfd,*), wk(nwk)
c
c  declare local variables.
c
      integer  i, ibeg, iend, nless1
      real  zero
      save zero
      data  zero /0./
c
c  validity-check arguments.
c
c***first executable statement  pchic
      if ( n.lt.2 )  go to 5001
      if ( incfd.lt.1 )  go to 5002
      do 1  i = 2, n
         if ( x(i).le.x(i-1) )  go to 5003
    1 continue
c
      ibeg = ic(1)
      iend = ic(2)
      ierr = 0
      if (abs(ibeg) .gt. 5)  ierr = ierr - 1
      if (abs(iend) .gt. 5)  ierr = ierr - 2
      if (ierr .lt. 0)  go to 5004
c
c  function definition is ok -- go on.
c
      nless1 = n - 1
      if ( nwk .lt. 2*nless1 )  go to 5007
c
c  set up h and slope arrays.
c
      do 20  i = 1, nless1
         wk(i) = x(i+1) - x(i)
         wk(nless1+i) = (f(1,i+1) - f(1,i)) / wk(i)
   20 continue
c
c  special case n=2 -- use linear interpolation.
c
      if (nless1 .gt. 1)  go to 1000
      d(1,1) = wk(2)
      d(1,n) = wk(2)
      go to 3000
c
c  normal case  (n .ge. 3) .
c
 1000 continue
c
c  set interior derivatives and default end conditions.
c
c     --------------------------------------
      call pchci (n, wk(1), wk(n), d, incfd)
c     --------------------------------------
c
c  set derivatives at points where monotonicity switches direction.
c
      if (switch .eq. zero)  go to 3000
c     ----------------------------------------------------
      call pchcs (switch, n, wk(1), wk(n), d, incfd, ierr)
c     ----------------------------------------------------
      if (ierr .ne. 0)  go to 5008
c
c  set end conditions.
c
 3000 continue
      if ( (ibeg.eq.0) .and. (iend.eq.0) )  go to 5000
c     -------------------------------------------------------
      call pchce (ic, vc, n, x, wk(1), wk(n), d, incfd, ierr)
c     -------------------------------------------------------
      if (ierr .lt. 0)  go to 5009
c
c  normal return.
c
 5000 continue
      return
c
c  error returns.
c
 5001 continue
c     n.lt.2 return.
      ierr = -1
      call xermsg ('slatec', 'pchic',
     +   'number of data points less than two', ierr, 1)
      return
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'pchic', 'increment less than one', ierr,
     +   1)
      return
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'pchic', 'x-array not strictly increasing'
     +   , ierr, 1)
      return
c
 5004 continue
c     ic out of range return.
      ierr = ierr - 3
      call xermsg ('slatec', 'pchic', 'ic out of range', ierr, 1)
      return
c
 5007 continue
c     nwk .lt. 2*(n-1)  return.
      ierr = -7
      call xermsg ('slatec', 'pchic', 'work array too small', ierr, 1)
      return
c
 5008 continue
c     error return from pchcs.
      ierr = -8
      call xermsg ('slatec', 'pchic', 'error return from pchcs', ierr,
     +   1)
      return
c
 5009 continue
c     error return from pchce.
c   *** this case should never occur ***
      ierr = -9
      call xermsg ('slatec', 'pchic', 'error return from pchce', ierr,
     +   1)
      return
c------------- last line of pchic follows ------------------------------
      end
