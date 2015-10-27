*deck dpchia
      double precision function dpchia (n, x, f, d, incfd, skip, a, b,
     +   ierr)
c***begin prologue  dpchia
c***purpose  evaluate the definite integral of a piecewise cubic
c            hermite function over an arbitrary interval.
c***library   slatec (pchip)
c***category  e3, h2a1b2
c***type      double precision (pchia-s, dpchia-d)
c***keywords  cubic hermite interpolation, numerical integration, pchip,
c             quadrature
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c          dpchia:  piecewise cubic hermite integrator, arbitrary limits
c
c     evaluates the definite integral of the cubic hermite function
c     defined by  n, x, f, d  over the interval [a, b].
c
c     to provide compatibility with dpchim and dpchic, includes an
c     increment between successive values of the f- and d-arrays.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  n, ierr
c        double precision  x(n), f(incfd,n), d(incfd,n), a, b
c        double precision  value, dpchia
c        logical  skip
c
c        value = dpchia (n, x, f, d, incfd, skip, a, b, ierr)
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
c           skip will be set to .true. on return with ierr.ge.0 .
c
c     a,b -- (input) the limits of integration.
c           note:  there is no requirement that [a,b] be contained in
c                  [x(1),x(n)].  however, the resulting integral value
c                  will be highly suspect, if not.
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           warning errors:
c              ierr = 1  if  a  is outside the interval [x(1),x(n)].
c              ierr = 2  if  b  is outside the interval [x(1),x(n)].
c              ierr = 3  if both of the above are true.  (note that this
c                        means that either [a,b] contains data interval
c                        or the intervals do not intersect at all.)
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c                (value will be zero in any of these cases.)
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c              ierr = -4  in case of an error return from dpchid (which
c                         should never occur).
c
c***references  (none)
c***routines called  dchfie, dpchid, xermsg
c***revision history  (yymmdd)
c   820730  date written
c   820804  converted to slatec library version.
c   870707  corrected xerror calls for d.p. name(s).
c   870707  corrected conversion to double precision.
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
c   930503  corrected to set value=0 when ierr.lt.0.  (fnf)
c   930504  changed dchfiv to dchfie.  (fnf)
c***end prologue  dpchia
c
c  programming notes:
c  1. the error flag from dpchid is tested, because a logic flaw
c     could conceivably result in ierd=-4, which should be reported.
c**end
c
c  declare arguments.
c
      integer  n, incfd, ierr
      double precision  x(*), f(incfd,*), d(incfd,*), a, b
      logical  skip
c
c  declare local variables.
c
      integer  i, ia, ib, ierd, il, ir
      double precision  value, xa, xb, zero
      save zero
      double precision  dchfie, dpchid
c
c  initialize.
c
      data  zero /0.d0/
c***first executable statement  dpchia
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
      ierr = 0
      if ( (a.lt.x(1)) .or. (a.gt.x(n)) )  ierr = ierr + 1
      if ( (b.lt.x(1)) .or. (b.gt.x(n)) )  ierr = ierr + 2
c
c  compute integral value.
c
      if (a .ne. b)  then
         xa = min (a, b)
         xb = max (a, b)
         if (xb .le. x(2))  then
c           interval is to left of x(2), so use first cubic.
c                   ---------------------------------------
            value = dchfie (x(1),x(2), f(1,1),f(1,2),
     +                                 d(1,1),d(1,2), a, b)
c                   ---------------------------------------
         else if (xa .ge. x(n-1))  then
c           interval is to right of x(n-1), so use last cubic.
c                   ------------------------------------------
            value = dchfie(x(n-1),x(n), f(1,n-1),f(1,n),
     +                                  d(1,n-1),d(1,n), a, b)
c                   ------------------------------------------
         else
c           'normal' case -- xa.lt.xb, xa.lt.x(n-1), xb.gt.x(2).
c      ......locate ia and ib such that
c               x(ia-1).lt.xa.le.x(ia).le.x(ib).le.xb.le.x(ib+1)
            ia = 1
            do 10  i = 1, n-1
               if (xa .gt. x(i))  ia = i + 1
   10       continue
c             ia = 1 implies xa.lt.x(1) .  otherwise,
c             ia is largest index such that x(ia-1).lt.xa,.
c
            ib = n
            do 20  i = n, ia, -1
               if (xb .lt. x(i))  ib = i - 1
   20       continue
c             ib = n implies xb.gt.x(n) .  otherwise,
c             ib is smallest index such that xb.lt.x(ib+1) .
c
c     ......compute the integral.
            if (ib .lt. ia)  then
c              this means ib = ia-1 and
c                 (a,b) is a subset of (x(ib),x(ia)).
c                      -------------------------------------------
               value = dchfie (x(ib),x(ia), f(1,ib),f(1,ia),
     +                                      d(1,ib),d(1,ia), a, b)
c                      -------------------------------------------
            else
c
c              first compute integral over (x(ia),x(ib)).
c                (case (ib .eq. ia) is taken care of by initialization
c                 of value to zero.)
               if (ib .gt. ia)  then
c                         ---------------------------------------------
                  value = dpchid (n, x, f, d, incfd, skip, ia, ib, ierd)
c                         ---------------------------------------------
                  if (ierd .lt. 0)  go to 5004
               endif
c
c              then add on integral over (xa,x(ia)).
               if (xa .lt. x(ia))  then
                  il = max(1, ia-1)
                  ir = il + 1
c                                 -------------------------------------
                  value = value + dchfie (x(il),x(ir), f(1,il),f(1,ir),
     +                                      d(1,il),d(1,ir), xa, x(ia))
c                                 -------------------------------------
               endif
c
c              then add on integral over (x(ib),xb).
               if (xb .gt. x(ib))  then
                  ir = min (ib+1, n)
                  il = ir - 1
c                                 -------------------------------------
                  value = value + dchfie (x(il),x(ir), f(1,il),f(1,ir),
     +                                      d(1,il),d(1,ir), x(ib), xb)
c                                 -------------------------------------
               endif
c
c              finally, adjust sign if necessary.
               if (a .gt. b)  value = -value
            endif
         endif
      endif
c
c  normal return.
c
 5000 continue
      dpchia = value
      return
c
c  error returns.
c
 5001 continue
c     n.lt.2 return.
      ierr = -1
      call xermsg ('slatec', 'dpchia',
     +   'number of data points less than two', ierr, 1)
      go to 5000
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'dpchia', 'increment less than one', ierr,
     +   1)
      go to 5000
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'dpchia',
     +   'x-array not strictly increasing', ierr, 1)
      go to 5000
c
 5004 continue
c     trouble in dpchid.  (should never occur.)
      ierr = -4
      call xermsg ('slatec', 'dpchia', 'trouble in dpchid', ierr, 1)
      go to 5000
c------------- last line of dpchia follows -----------------------------
      end
