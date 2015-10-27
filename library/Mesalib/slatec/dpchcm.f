*deck dpchcm
      subroutine dpchcm (n, x, f, d, incfd, skip, ismon, ierr)
c***begin prologue  dpchcm
c***purpose  check a cubic hermite function for monotonicity.
c***library   slatec (pchip)
c***category  e3
c***type      double precision (pchcm-s, dpchcm-d)
c***keywords  cubic hermite interpolation, monotone interpolation,
c             pchip, piecewise cubic interpolation, utility routine
c***author  fritsch, f. n., (llnl)
c             computing & mathematics research division
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c *usage:
c
c        parameter  (incfd = ...)
c        integer  n, ismon(n), ierr
c        double precision  x(n), f(incfd,n), d(incfd,n)
c        logical  skip
c
c        call  dpchcm (n, x, f, d, incfd, skip, ismon, ierr)
c
c *arguments:
c
c     n:in  is the number of data points.  (error return if n.lt.2 .)
c
c     x:in  is a real*8 array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.
c           (error return if not.)
c
c     f:in  is a real*8 array of function values.  f(1+(i-1)*incfd) is
c           the value corresponding to x(i).
c
c     d:in  is a real*8 array of derivative values.  d(1+(i-1)*incfd) is
c           is the value corresponding to x(i).
c
c     incfd:in  is the increment between successive values in f and d.
c           (error return if  incfd.lt.1 .)
c
c     skip:inout  is a logical variable which should be set to
c           .true. if the user wishes to skip checks for validity of
c           preceding parameters, or to .false. otherwise.
c           this will save time in case these checks have already
c           been performed.
c           skip will be set to .true. on normal return.
c
c     ismon:out  is an integer array indicating on which intervals the
c           pch function defined by  n, x, f, d  is monotonic.
c           for data interval [x(i),x(i+1)],
c             ismon(i) = -3  if function is probably decreasing;
c             ismon(i) = -1  if function is strictly decreasing;
c             ismon(i) =  0  if function is constant;
c             ismon(i) =  1  if function is strictly increasing;
c             ismon(i) =  2  if function is non-monotonic;
c             ismon(i) =  3  if function is probably increasing.
c                if abs(ismon)=3, this means that the d-values are near
c                the boundary of the monotonicity region.  a small
c                increase produces non-monotonicity; decrease, strict
c                monotonicity.
c           the above applies to i=1(1)n-1.  ismon(n) indicates whether
c              the entire function is monotonic on [x(1),x(n)].
c
c     ierr:out  is an error flag.
c           normal return:
c              ierr = 0  (no errors).
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c          (the ismon-array has not been changed in any of these cases.)
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c
c *description:
c
c          dpchcm:  piecewise cubic hermite -- check monotonicity.
c
c     checks the piecewise cubic hermite function defined by  n,x,f,d
c     for monotonicity.
c
c     to provide compatibility with dpchim and dpchic, includes an
c     increment between successive values of the f- and d-arrays.
c
c *cautions:
c     this provides the same capability as old dpchmc, except that a
c     new output value, -3, was added february 1989.  (formerly, -3
c     and +3 were lumped together in the single value 3.)  codes that
c     flag nonmonotonicity by "if (ismon.eq.2)" need not be changed.
c     codes that check via "if (ismon.ge.3)" should change the test to
c     "if (iabs(ismon).ge.3)".  codes that declare monotonicity via
c     "if (ismon.le.1)" should change to "if (iabs(ismon).le.1)".
c
c***references  f. n. fritsch and r. e. carlson, monotone piecewise
c                 cubic interpolation, siam journal on numerical ana-
c                 lysis 17, 2 (april 1980), pp. 238-246.
c***routines called  dchfcm, xermsg
c***revision history  (yymmdd)
c   820518  date written
c   820804  converted to slatec library version.
c   831201  reversed order of subscripts of f and d, so that the
c           routine will work properly when incfd.gt.1 .  (bug!)
c   870707  corrected xerror calls for d.p. name(s).
c   890206  corrected xerror calls.
c   890209  added possible ismon value of -3 and modified code so
c           that 1,3,-1 produces ismon(n)=2, rather than 3.
c   890306  added caution about changed output.
c   890407  changed name from dpchmc to dpchcm, as requested at the
c           march 1989 slatec cml meeting, and made a few other
c           minor modifications necessitated by this change.
c   890407  converted to new slatec format.
c   890407  modified description to ldoc format.
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920429  revised format and order of references.  (wrb,fnf)
c***end prologue  dpchcm
c
c  fortran intrinsics used:  isign.
c  other routines used:  chfcm, xermsg.
c
c ----------------------------------------------------------------------
c
c  programming notes:
c
c     an alternate organization would have separate loops for computing
c     ismon(i), i=1,...,nseg, and for the computation of ismon(n).  the
c     first loop can be readily parallelized, since the nseg calls to
c     chfcm are independent.  the second loop can be cut short if
c     ismon(n) is ever equal to 2, for it cannot be changed further.
c
c     to produce a single precision version, simply:
c        a. change dpchcm to pchcm wherever it occurs,
c        b. change dchfcm to chfcm wherever it occurs, and
c        c. change the double precision declarations to real.
c
c  declare arguments.
c
      integer n, incfd, ismon(n), ierr
      double precision  x(n), f(incfd,n), d(incfd,n)
      logical  skip
c
c  declare local variables.
c
      integer i, nseg
      double precision  delta
      integer dchfcm
c
c  validity-check arguments.
c
c***first executable statement  dpchcm
      if (skip)  go to 5
c
      if ( n.lt.2 )  go to 5001
      if ( incfd.lt.1 )  go to 5002
      do 1  i = 2, n
         if ( x(i).le.x(i-1) )  go to 5003
    1 continue
      skip = .true.
c
c  function definition is ok -- go on.
c
    5 continue
      nseg = n - 1
      do 90  i = 1, nseg
         delta = (f(1,i+1)-f(1,i))/(x(i+1)-x(i))
c                   -------------------------------
         ismon(i) = dchfcm (d(1,i), d(1,i+1), delta)
c                   -------------------------------
         if (i .eq. 1)  then
            ismon(n) = ismon(1)
         else
c           need to figure out cumulative monotonicity from following
c           "multiplication table":
c
c                    +        i s m o n (i)
c                     +  -3  -1   0   1   3   2
c                      +------------------------+
c               i   -3 i -3  -3  -3   2   2   2 i
c               s   -1 i -3  -1  -1   2   2   2 i
c               m    0 i -3  -1   0   1   3   2 i
c               o    1 i  2   2   1   1   3   2 i
c               n    3 i  2   2   3   3   3   2 i
c              (n)   2 i  2   2   2   2   2   2 i
c                      +------------------------+
c           note that the 2 row and column are out of order so as not
c           to obscure the symmetry in the rest of the table.
c
c           no change needed if equal or constant on this interval or
c           already declared nonmonotonic.
            if ( (ismon(i).ne.ismon(n)) .and. (ismon(i).ne.0)
     .                                  .and. (ismon(n).ne.2) )  then
               if ( (ismon(i).eq.2) .or. (ismon(n).eq.0) )  then
                  ismon(n) =  ismon(i)
               else if (ismon(i)*ismon(n) .lt. 0)  then
c                 this interval has opposite sense from curve so far.
                  ismon(n) = 2
               else
c                 at this point, both are nonzero with same sign, and
c                 we have already eliminated case both +-1.
                  ismon(n) = isign (3, ismon(n))
               endif
            endif
         endif
   90 continue
c
c  normal return.
c
      ierr = 0
      return
c
c  error returns.
c
 5001 continue
c     n.lt.2 return.
      ierr = -1
      call xermsg ('slatec', 'dpchcm',
     +   'number of data points less than two', ierr, 1)
      return
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'dpchcm', 'increment less than one', ierr,
     +   1)
      return
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'dpchcm',
     +   'x-array not strictly increasing', ierr, 1)
      return
c------------- last line of dpchcm follows -----------------------------
      end
