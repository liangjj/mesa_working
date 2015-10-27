*deck dpchid
      double precision function dpchid (n, x, f, d, incfd, skip, ia, ib,
     +   ierr)
c***begin prologue  dpchid
c***purpose  evaluate the definite integral of a piecewise cubic
c            hermite function over an interval whose endpoints are data
c            points.
c***library   slatec (pchip)
c***category  e3, h2a1b2
c***type      double precision (pchid-s, dpchid-d)
c***keywords  cubic hermite interpolation, numerical integration, pchip,
c             quadrature
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c          dpchid:  piecewise cubic hermite integrator, data limits
c
c     evaluates the definite integral of the cubic hermite function
c     defined by  n, x, f, d  over the interval [x(ia), x(ib)].
c
c     to provide compatibility with dpchim and dpchic, includes an
c     increment between successive values of the f- and d-arrays.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  n, ia, ib, ierr
c        double precision  x(n), f(incfd,n), d(incfd,n)
c        logical  skip
c
c        value = dpchid (n, x, f, d, incfd, skip, ia, ib, ierr)
c
c   parameters:
c
c     value -- (output) value of the requested integral.
c
c     n -- (input) number of data points.  (error return if n.lt.2 .)
c
c     x -- (input) real*8 array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.
c           (error return if not.)
c
c     f -- (input) real*8 array of function values.  f(1+(i-1)*incfd) is
c           the value corresponding to x(i).
c
c     d -- (input) real*8 array of derivative values.  d(1+(i-1)*incfd)
c           is the value corresponding to x(i).
c
c     incfd -- (input) increment between successive values in f and d.
c           (error return if  incfd.lt.1 .)
c
c     skip -- (input/output) logical variable which should be set to
c           .true. if the user wishes to skip checks for validity of
c           preceding parameters, or to .false. otherwise.
c           this will save time in case these checks have already
c           been performed (say, in dpchim or dpchic).
c           skip will be set to .true. on return with ierr = 0 or -4.
c
c     ia,ib -- (input) indices in x-array for the limits of integration.
c           both must be in the range [1,n].  (error return if not.)
c           no restrictions on their relative values.
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c              ierr = -4  if ia or ib is out of range.
c                (value will be zero in any of these cases.)
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   820723  date written
c   820804  converted to slatec library version.
c   870707  corrected xerror calls for d.p. name(s).
c   870813  minor cosmetic changes.
c   890206  corrected xerror calls.
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890703  corrected category record.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   930504  corrected to set value=0 when ierr.ne.0.  (fnf)
c***end prologue  dpchid
c
c  programming notes:
c  1. this routine uses a special formula that is valid only for
c     integrals whose limits coincide with data values.  this is
c     mathematically equivalent to, but much more efficient than,
c     calls to dchfie.
c**end
c
c  declare arguments.
c
      integer  n, incfd, ia, ib, ierr
      double precision  x(*), f(incfd,*), d(incfd,*)
      logical  skip
c
c  declare local variables.
c
      integer  i, iup, low
      double precision  h, half, six, sum, value, zero
      save zero, half, six
c
c  initialize.
c
      data  zero /0.d0/,  half/.5d0/, six/6.d0/
c***first executable statement  dpchid
      value = zero
c
c  validity-check arguments.
c
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
      skip = .true.
      if ((ia.lt.1) .or. (ia.gt.n))  go to 5004
      if ((ib.lt.1) .or. (ib.gt.n))  go to 5004
      ierr = 0
c
c  compute integral value.
c
      if (ia .ne. ib)  then
         low = min(ia, ib)
         iup = max(ia, ib) - 1
         sum = zero
         do 10  i = low, iup
            h = x(i+1) - x(i)
            sum = sum + h*( (f(1,i) + f(1,i+1)) +
     *                      (d(1,i) - d(1,i+1))*(h/six) )
   10    continue
         value = half * sum
         if (ia .gt. ib)  value = -value
      endif
c
c  normal return.
c
 5000 continue
      dpchid = value
      return
c
c  error returns.
c
 5001 continue
c     n.lt.2 return.
      ierr = -1
      call xermsg ('slatec', 'dpchid',
     +   'number of data points less than two', ierr, 1)
      go to 5000
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'dpchid', 'increment less than one', ierr,
     +   1)
      go to 5000
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'dpchid',
     +   'x-array not strictly increasing', ierr, 1)
      go to 5000
c
 5004 continue
c     ia or ib out of range return.
      ierr = -4
      call xermsg ('slatec', 'dpchid', 'ia or ib out of range', ierr,
     +   1)
      go to 5000
c------------- last line of dpchid follows -----------------------------
      end
