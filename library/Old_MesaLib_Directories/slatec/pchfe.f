*deck pchfe
      subroutine pchfe (n, x, f, d, incfd, skip, ne, xe, fe, ierr)
c***begin prologue  pchfe
c***purpose  evaluate a piecewise cubic hermite function at an array of
c            points.  may be used by itself for hermite interpolation,
c            or as an evaluator for pchim or pchic.
c***library   slatec (pchip)
c***category  e3
c***type      single precision (pchfe-s, dpchfe-d)
c***keywords  cubic hermite evaluation, hermite interpolation, pchip,
c             piecewise cubic evaluation
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c          pchfe:  piecewise cubic hermite function evaluator
c
c     evaluates the cubic hermite function defined by  n, x, f, d  at
c     the points  xe(j), j=1(1)ne.
c
c     to provide compatibility with pchim and pchic, includes an
c     increment between successive values of the f- and d-arrays.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  n, ne, ierr
c        real  x(n), f(incfd,n), d(incfd,n), xe(ne), fe(ne)
c        logical  skip
c
c        call  pchfe (n, x, f, d, incfd, skip, ne, xe, fe, ierr)
c
c   parameters:
c
c     n -- (input) number of data points.  (error return if n.lt.2 .)
c
c     x -- (input) real array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.
c           (error return if not.)
c
c     f -- (input) real array of function values.  f(1+(i-1)*incfd) is
c           the value corresponding to x(i).
c
c     d -- (input) real array of derivative values.  d(1+(i-1)*incfd) is
c           the value corresponding to x(i).
c
c     incfd -- (input) increment between successive values in f and d.
c           (error return if  incfd.lt.1 .)
c
c     skip -- (input/output) logical variable which should be set to
c           .true. if the user wishes to skip checks for validity of
c           preceding parameters, or to .false. otherwise.
c           this will save time in case these checks have already
c           been performed (say, in pchim or pchic).
c           skip will be set to .true. on normal return.
c
c     ne -- (input) number of evaluation points.  (error return if
c           ne.lt.1 .)
c
c     xe -- (input) real array of points at which the function is to be
c           evaluated.
c
c          notes:
c           1. the evaluation will be most efficient if the elements
c              of xe are increasing relative to x;
c              that is,   xe(j) .ge. x(i)
c              implies    xe(k) .ge. x(i),  all k.ge.j .
c           2. if any of the xe are outside the interval [x(1),x(n)],
c              values are extrapolated from the nearest extreme cubic,
c              and a warning error is returned.
c
c     fe -- (output) real array of values of the cubic hermite function
c           defined by  n, x, f, d  at the points  xe.
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           warning error:
c              ierr.gt.0  means that extrapolation was performed at
c                 ierr points.
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c              ierr = -4  if ne.lt.1 .
c             (the fe-array has not been changed in any of these cases.)
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c
c***references  (none)
c***routines called  chfev, xermsg
c***revision history  (yymmdd)
c   811020  date written
c   820803  minor cosmetic changes for release 1.
c   870707  minor cosmetic changes to prologue.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  pchfe
c  programming notes:
c
c     1. to produce a double precision version, simply:
c        a. change pchfe to dpchfe, and chfev to dchfev, wherever they
c           occur,
c        b. change the real declaration to double precision,
c
c     2. most of the coding between the call to chfev and the end of
c        the ir-loop could be eliminated if it were permissible to
c        assume that xe is ordered relative to x.
c
c     3. chfev does not assume that x1 is less than x2.  thus, it would
c        be possible to write a version of pchfe that assumes a strict-
c        ly decreasing x-array by simply running the ir-loop backwards
c        (and reversing the order of appropriate tests).
c
c     4. the present code has a minor bug, which i have decided is not
c        worth the effort that would be required to fix it.
c        if xe contains points in [x(n-1),x(n)], followed by points .lt.
c        x(n-1), followed by points .gt.x(n), the extrapolation points
c        will be counted (at least) twice in the total returned in ierr.
c
c  declare arguments.
c
      integer  n, incfd, ne, ierr
      real  x(*), f(incfd,*), d(incfd,*), xe(*), fe(*)
      logical  skip
c
c  declare local variables.
c
      integer  i, ierc, ir, j, jfirst, next(2), nj
c
c  validity-check arguments.
c
c***first executable statement  pchfe
      if (skip)  go to 5
c
      if ( n.lt.2 )  go to 5001
      if ( incfd.lt.1 )  go to 5002
      do 1  i = 2, n
         if ( x(i).le.x(i-1) )  go to 5003
    1 continue
c
c  function definition is ok, go on.
c
    5 continue
      if ( ne.lt.1 )  go to 5004
      ierr = 0
      skip = .true.
c
c  loop over intervals.        (   interval index is  il = ir-1  . )
c                              ( interval is x(il).le.x.lt.x(ir) . )
      jfirst = 1
      ir = 2
   10 continue
c
c     skip out of loop if have processed all evaluation points.
c
         if (jfirst .gt. ne)  go to 5000
c
c     locate all points in interval.
c
         do 20  j = jfirst, ne
            if (xe(j) .ge. x(ir))  go to 30
   20    continue
         j = ne + 1
         go to 40
c
c     have located first point beyond interval.
c
   30    continue
         if (ir .eq. n)  j = ne + 1
c
   40    continue
         nj = j - jfirst
c
c     skip evaluation if no points in interval.
c
         if (nj .eq. 0)  go to 50
c
c     evaluate cubic at xe(i),  i = jfirst (1) j-1 .
c
c       ----------------------------------------------------------------
        call chfev (x(ir-1),x(ir), f(1,ir-1),f(1,ir), d(1,ir-1),d(1,ir),
     *              nj, xe(jfirst), fe(jfirst), next, ierc)
c       ----------------------------------------------------------------
         if (ierc .lt. 0)  go to 5005
c
         if (next(2) .eq. 0)  go to 42
c        if (next(2) .gt. 0)  then
c           in the current set of xe-points, there are next(2) to the
c           right of x(ir).
c
            if (ir .lt. n)  go to 41
c           if (ir .eq. n)  then
c              these are actually extrapolation points.
               ierr = ierr + next(2)
               go to 42
   41       continue
c           else
c              we should never have gotten here.
               go to 5005
c           endif
c        endif
   42    continue
c
         if (next(1) .eq. 0)  go to 49
c        if (next(1) .gt. 0)  then
c           in the current set of xe-points, there are next(1) to the
c           left of x(ir-1).
c
            if (ir .gt. 2)  go to 43
c           if (ir .eq. 2)  then
c              these are actually extrapolation points.
               ierr = ierr + next(1)
               go to 49
   43       continue
c           else
c              xe is not ordered relative to x, so must adjust
c              evaluation interval.
c
c              first, locate first point to left of x(ir-1).
               do 44  i = jfirst, j-1
                  if (xe(i) .lt. x(ir-1))  go to 45
   44          continue
c              note-- cannot drop through here unless there is an error
c                     in chfev.
               go to 5005
c
   45          continue
c              reset j.  (this will be the new jfirst.)
               j = i
c
c              now find out how far to back up in the x-array.
               do 46  i = 1, ir-1
                  if (xe(j) .lt. x(i)) go to 47
   46          continue
c              nb-- can never drop through here, since xe(j).lt.x(ir-1).
c
   47          continue
c              at this point, either  xe(j) .lt. x(1)
c                 or      x(i-1) .le. xe(j) .lt. x(i) .
c              reset ir, recognizing that it will be incremented before
c              cycling.
               ir = max(1, i-1)
c           endif
c        endif
   49    continue
c
         jfirst = j
c
c     end of ir-loop.
c
   50 continue
      ir = ir + 1
      if (ir .le. n)  go to 10
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
      call xermsg ('slatec', 'pchfe',
     +   'number of data points less than two', ierr, 1)
      return
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'pchfe', 'increment less than one', ierr,
     +   1)
      return
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'pchfe', 'x-array not strictly increasing'
     +   , ierr, 1)
      return
c
 5004 continue
c     ne.lt.1 return.
      ierr = -4
      call xermsg ('slatec', 'pchfe',
     +   'number of evaluation points less than one', ierr, 1)
      return
c
 5005 continue
c     error return from chfev.
c   *** this case should never occur ***
      ierr = -5
      call xermsg ('slatec', 'pchfe',
     +   'error return from chfev -- fatal', ierr, 2)
      return
c------------- last line of pchfe follows ------------------------------
      end
