*deck dpchim
      subroutine dpchim (n, x, f, d, incfd, ierr)
c***begin prologue  dpchim
c***purpose  set derivatives needed to determine a monotone piecewise
c            cubic hermite interpolant to given data.  boundary values
c            are provided which are compatible with monotonicity.  the
c            interpolant will have an extremum at each point where mono-
c            tonicity switches direction.  (see dpchic if user control
c            is desired over boundary or switch conditions.)
c***library   slatec (pchip)
c***category  e1a
c***type      double precision (pchim-s, dpchim-d)
c***keywords  cubic hermite interpolation, monotone interpolation,
c             pchip, piecewise cubic interpolation
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c          dpchim:  piecewise cubic hermite interpolation to
c                  monotone data.
c
c     sets derivatives needed to determine a monotone piecewise cubic
c     hermite interpolant to the data given in x and f.
c
c     default boundary conditions are provided which are compatible
c     with monotonicity.  (see dpchic if user control of boundary con-
c     ditions is desired.)
c
c     if the data are only piecewise monotonic, the interpolant will
c     have an extremum at each point where monotonicity switches direc-
c     tion.  (see dpchic if user control is desired in such cases.)
c
c     to facilitate two-dimensional applications, includes an increment
c     between successive values of the f- and d-arrays.
c
c     the resulting piecewise cubic hermite function may be evaluated
c     by dpchfe or dpchfd.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        parameter  (incfd = ...)
c        integer  n, ierr
c        double precision  x(n), f(incfd,n), d(incfd,n)
c
c        call  dpchim (n, x, f, d, incfd, ierr)
c
c   parameters:
c
c     n -- (input) number of data points.  (error return if n.lt.2 .)
c           if n=2, simply does linear interpolation.
c
c     x -- (input) real*8 array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.
c           (error return if not.)
c
c     f -- (input) real*8 array of dependent variable values to be
c           interpolated.  f(1+(i-1)*incfd) is value corresponding to
c           x(i).  dpchim is designed for monotonic data, but it will
c           work for any f-array.  it will force extrema at points where
c           monotonicity switches direction.  if some other treatment of
c           switch points is desired, dpchic should be used instead.
c                                     -----
c     d -- (output) real*8 array of derivative values at the data
c           points.  if the data are monotonic, these values will
c           determine a monotone cubic hermite function.
c           the value corresponding to x(i) is stored in
c                d(1+(i-1)*incfd),  i=1(1)n.
c           no other entries in d are changed.
c
c     incfd -- (input) increment between successive values in f and d.
c           this argument is provided primarily for 2-d applications.
c           (error return if  incfd.lt.1 .)
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           warning error:
c              ierr.gt.0  means that ierr switches in the direction
c                 of monotonicity were detected.
c           "recoverable" errors:
c              ierr = -1  if n.lt.2 .
c              ierr = -2  if incfd.lt.1 .
c              ierr = -3  if the x-array is not strictly increasing.
c             (the d-array has not been changed in any of these cases.)
c               note:  the above errors are checked in the order listed,
c                   and following arguments have **not** been validated.
c
c***references  1. f. n. fritsch and j. butland, a method for construc-
c                 ting local monotone piecewise cubic interpolants, siam
c                 journal on scientific and statistical computing 5, 2
c                 (june 1984), pp. 300-304.
c               2. f. n. fritsch and r. e. carlson, monotone piecewise
c                 cubic interpolation, siam journal on numerical ana-
c                 lysis 17, 2 (april 1980), pp. 238-246.
c***routines called  dpchst, xermsg
c***revision history  (yymmdd)
c   811103  date written
c   820201  1. introduced  dpchst  to reduce possible over/under-
c             flow problems.
c           2. rearranged derivative formula for same reason.
c   820602  1. modified end conditions to be continuous functions
c             of data when monotonicity switches in next interval.
c           2. modified formulas so end conditions are less prone
c             of over/underflow problems.
c   820803  minor cosmetic changes for release 1.
c   870707  corrected xerror calls for d.p. name(s).
c   870813  updated reference 1.
c   890206  corrected xerror calls.
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890703  corrected category record.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920429  revised format and order of references.  (wrb,fnf)
c***end prologue  dpchim
c  programming notes:
c
c     1. the function  dpchst(arg1,arg2)  is assumed to return zero if
c        either argument is zero, +1 if they are of the same sign, and
c        -1 if they are of opposite sign.
c     2. to produce a single precision version, simply:
c        a. change dpchim to pchim wherever it occurs,
c        b. change dpchst to pchst wherever it occurs,
c        c. change all references to the fortran intrinsics to their
c           single precision equivalents,
c        d. change the double precision declarations to real, and
c        e. change the constants zero and three to single precision.
c
c  declare arguments.
c
      integer  n, incfd, ierr
      double precision  x(*), f(incfd,*), d(incfd,*)
c
c  declare local variables.
c
      integer  i, nless1
      double precision  del1, del2, dmax, dmin, drat1, drat2, dsave,
     *      h1, h2, hsum, hsumt3, three, w1, w2, zero
      save zero, three
      double precision  dpchst
      data  zero /0.d0/, three/3.d0/
c
c  validity-check arguments.
c
c***first executable statement  dpchim
      if ( n.lt.2 )  go to 5001
      if ( incfd.lt.1 )  go to 5002
      do 1  i = 2, n
         if ( x(i).le.x(i-1) )  go to 5003
    1 continue
c
c  function definition is ok, go on.
c
      ierr = 0
      nless1 = n - 1
      h1 = x(2) - x(1)
      del1 = (f(1,2) - f(1,1))/h1
      dsave = del1
c
c  special case n=2 -- use linear interpolation.
c
      if (nless1 .gt. 1)  go to 10
      d(1,1) = del1
      d(1,n) = del1
      go to 5000
c
c  normal case  (n .ge. 3).
c
   10 continue
      h2 = x(3) - x(2)
      del2 = (f(1,3) - f(1,2))/h2
c
c  set d(1) via non-centered three-point formula, adjusted to be
c     shape-preserving.
c
      hsum = h1 + h2
      w1 = (h1 + hsum)/hsum
      w2 = -h1/hsum
      d(1,1) = w1*del1 + w2*del2
      if ( dpchst(d(1,1),del1) .le. zero)  then
         d(1,1) = zero
      else if ( dpchst(del1,del2) .lt. zero)  then
c        need do this check only if monotonicity switches.
         dmax = three*del1
         if (abs(d(1,1)) .gt. abs(dmax))  d(1,1) = dmax
      endif
c
c  loop through interior points.
c
      do 50  i = 2, nless1
         if (i .eq. 2)  go to 40
c
         h1 = h2
         h2 = x(i+1) - x(i)
         hsum = h1 + h2
         del1 = del2
         del2 = (f(1,i+1) - f(1,i))/h2
   40    continue
c
c        set d(i)=0 unless data are strictly monotonic.
c
         d(1,i) = zero
         if ( dpchst(del1,del2) )  42, 41, 45
c
c        count number of changes in direction of monotonicity.
c
   41    continue
         if (del2 .eq. zero)  go to 50
         if ( dpchst(dsave,del2) .lt. zero)  ierr = ierr + 1
         dsave = del2
         go to 50
c
   42    continue
         ierr = ierr + 1
         dsave = del2
         go to 50
c
c        use brodlie modification of butland formula.
c
   45    continue
         hsumt3 = hsum+hsum+hsum
         w1 = (hsum + h1)/hsumt3
         w2 = (hsum + h2)/hsumt3
         dmax = max( abs(del1), abs(del2) )
         dmin = min( abs(del1), abs(del2) )
         drat1 = del1/dmax
         drat2 = del2/dmax
         d(1,i) = dmin/(w1*drat1 + w2*drat2)
c
   50 continue
c
c  set d(n) via non-centered three-point formula, adjusted to be
c     shape-preserving.
c
      w1 = -h2/hsum
      w2 = (h2 + hsum)/hsum
      d(1,n) = w1*del1 + w2*del2
      if ( dpchst(d(1,n),del2) .le. zero)  then
         d(1,n) = zero
      else if ( dpchst(del1,del2) .lt. zero)  then
c        need do this check only if monotonicity switches.
         dmax = three*del2
         if (abs(d(1,n)) .gt. abs(dmax))  d(1,n) = dmax
      endif
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
      call xermsg ('slatec', 'dpchim',
     +   'number of data points less than two', ierr, 1)
      return
c
 5002 continue
c     incfd.lt.1 return.
      ierr = -2
      call xermsg ('slatec', 'dpchim', 'increment less than one', ierr,
     +   1)
      return
c
 5003 continue
c     x-array not strictly increasing.
      ierr = -3
      call xermsg ('slatec', 'dpchim',
     +   'x-array not strictly increasing', ierr, 1)
      return
c------------- last line of dpchim follows -----------------------------
      end
