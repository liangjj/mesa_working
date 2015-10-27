*deck dpchbs
      subroutine dpchbs (n, x, f, d, incfd, knotyp, nknots, t, bcoef,
     +   ndim, kord, ierr)
c***begin prologue  dpchbs
c***purpose  piecewise cubic hermite to b-spline converter.
c***library   slatec (pchip)
c***category  e3
c***type      double precision (pchbs-s, dpchbs-d)
c***keywords  b-splines, conversion, cubic hermite interpolation,
c             piecewise cubic interpolation
c***author  fritsch, f. n., (llnl)
c             computing and mathematics research division
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c *usage:
c
c        integer  n, incfd, knotyp, nknots, ndim, kord, ierr
c        parameter  (incfd = ...)
c        double precision  x(nmax), f(incfd,nmax), d(incfd,nmax),
c       *      t(2*nmax+4), bcoef(2*nmax)
c
c        call dpchbs (n, x, f, d, incfd, knotyp, nknots, t, bcoef,
c       *             ndim, kord, ierr)
c
c *arguments:
c
c     n:in  is the number of data points, n.ge.2 .  (not checked)
c
c     x:in  is the real array of independent variable values.  the
c           elements of x must be strictly increasing:
c                x(i-1) .lt. x(i),  i = 2(1)n.   (not checked)
c           nmax, the dimension of x, must be .ge.n.
c
c     f:in  is the real array of dependent variable values.
c           f(1+(i-1)*incfd) is the value corresponding to x(i).
c           nmax, the second dimension of f, must be .ge.n.
c
c     d:in  is the real array of derivative values at the data points.
c           d(1+(i-1)*incfd) is the value corresponding to x(i).
c           nmax, the second dimension of d, must be .ge.n.
c
c     incfd:in  is the increment between successive values in f and d.
c           this argument is provided primarily for 2-d applications.
c           it may have the value 1 for one-dimensional applications,
c           in which case f and d may be singly-subscripted arrays.
c
c     knotyp:in  is a flag to control the knot sequence.
c           the knot sequence t is normally computed from x by putting
c           a double knot at each x and setting the end knot pairs
c           according to the value of knotyp:
c              knotyp = 0:  quadruple knots at x(1) and x(n).  (default)
c              knotyp = 1:  replicate lengths of extreme subintervals:
c                           t( 1 ) = t( 2 ) = x(1) - (x(2)-x(1))  ;
c                           t(m+4) = t(m+3) = x(n) + (x(n)-x(n-1)).
c              knotyp = 2:  periodic placement of boundary knots:
c                           t( 1 ) = t( 2 ) = x(1) - (x(n)-x(n-1));
c                           t(m+4) = t(m+3) = x(n) + (x(2)-x(1))  .
c              here m=ndim=2*n.
c           if the input value of knotyp is negative, however, it is
c           assumed that nknots and t were set in a previous call.
c           this option is provided for improved efficiency when used
c           in a parametric setting.
c
c     nknots:inout  is the number of knots.
c           if knotyp.ge.0, then nknots will be set to ndim+4.
c           if knotyp.lt.0, then nknots is an input variable, and an
c              error return will be taken if it is not equal to ndim+4.
c
c     t:inout  is the array of 2*n+4 knots for the b-representation.
c           if knotyp.ge.0, t will be returned by dpchbs with the
c              interior double knots equal to the x-values and the
c              boundary knots set as indicated above.
c           if knotyp.lt.0, it is assumed that t was set by a
c              previous call to dpchbs.  (this routine does **not**
c              verify that t forms a legitimate knot sequence.)
c
c     bcoef:out  is the array of 2*n b-spline coefficients.
c
c     ndim:out  is the dimension of the b-spline space.  (set to 2*n.)
c
c     kord:out  is the order of the b-spline.  (set to 4.)
c
c     ierr:out  is an error flag.
c           normal return:
c              ierr = 0  (no errors).
c           "recoverable" errors:
c              ierr = -4  if knotyp.gt.2 .
c              ierr = -5  if knotyp.lt.0 and nknots.ne.(2*n+4).
c
c *description:
c     dpchbs computes the b-spline representation of the pch function
c     determined by n,x,f,d.  to be compatible with the rest of pchip,
c     dpchbs includes incfd, the increment between successive values of
c     the f- and d-arrays.
c
c     the output is the b-representation for the function:  nknots, t,
c     bcoef, ndim, kord.
c
c *caution:
c     since it is assumed that the input pch function has been
c     computed by one of the other routines in the package pchip,
c     input arguments n, x, incfd are **not** checked for validity.
c
c *restrictions/assumptions:
c     1. n.ge.2 .  (not checked)
c     2. x(i).lt.x(i+1), i=1,...,n .  (not checked)
c     3. incfd.gt.0 .  (not checked)
c     4. knotyp.le.2 .  (error return if not)
c    *5. nknots = ndim+4 = 2*n+4 .  (error return if not)
c    *6. t(2*k+1) = t(2*k) = x(k), k=1,...,n .  (not checked)
c
c       * indicates this applies only if knotyp.lt.0 .
c
c *portability:
c     argument incfd is used only to cause the compiler to generate
c     efficient code for the subscript expressions (1+(i-1)*incfd) .
c     the normal usage, in which dpchbs is called with one-dimensional
c     arrays f and d, is probably non-fortran 77, in the strict sense,
c     but it works on all systems on which dpchbs has been tested.
c
c *see also:
c     pchic, pchim, or pchsp can be used to determine an interpolating
c        pch function from a set of data.
c     the b-spline routine dbvalu can be used to evaluate the
c        b-representation that is output by dpchbs.
c        (see bspdoc for more information.)
c
c***references  f. n. fritsch, "representations for parametric cubic
c                 splines," computer aided geometric design 6 (1989),
c                 pp.79-82.
c***routines called  dpchkt, xermsg
c***revision history  (yymmdd)
c   870701  date written
c   900405  converted fortran to upper case.
c   900405  removed requirement that x be dimensioned n+1.
c   900406  modified to make pchkt a subsidiary routine to simplify
c           usage.  in the process, added argument incfd to be com-
c           patible with the rest of pchip.
c   900410  converted prologue to slatec 4.0 format.
c   900410  added calls to xermsg and changed constant 3. to 3 to
c           reduce single/double differences.
c   900411  added reference.
c   900430  produced double precision version.
c   900501  corrected declarations.
c   930317  minor cosmetic changes.  (fnf)
c   930514  corrected problems with dimensioning of arguments and
c           clarified description.  (fnf)
c   930604  removed  nknots from dpchkt call list.  (fnf)
c***end prologue  dpchbs
c
c*internal notes:
c
c**end
c
c  declare arguments.
c
      integer  n, incfd, knotyp, nknots, ndim, kord, ierr
      double precision  x(*), f(incfd,*), d(incfd,*), t(*), bcoef(*)
c
c  declare local variables.
c
      integer  k, kk
      double precision  dov3, hnew, hold
      character*8  libnam, subnam
c***first executable statement  dpchbs
c
c  initialize.
c
      ndim = 2*n
      kord = 4
      ierr = 0
      libnam = 'slatec'
      subnam = 'dpchbs'
c
c  check argument validity.  set up knot sequence if ok.
c
      if ( knotyp.gt.2 )  then
         ierr = -1
         call xermsg (libnam, subnam, 'knotyp greater than 2', ierr, 1)
         return
      endif
      if ( knotyp.lt.0 )  then
         if ( nknots.ne.ndim+4 )  then
            ierr = -2
            call xermsg (libnam, subnam,
     *                    'knotyp.lt.0 and nknots.ne.(2*n+4)', ierr, 1)
            return
         endif
      else
c          set up knot sequence.
         nknots = ndim + 4
         call dpchkt (n, x, knotyp, t)
      endif
c
c  compute b-spline coefficients.
c
      hnew = t(3) - t(1)
      do 40  k = 1, n
         kk = 2*k
         hold = hnew
c          the following requires mixed mode arithmetic.
         dov3 = d(1,k)/3
         bcoef(kk-1) = f(1,k) - hold*dov3
c          the following assumes t(2*k+1) = x(k).
         hnew = t(kk+3) - t(kk+1)
         bcoef(kk) = f(1,k) + hnew*dov3
   40 continue
c
c  terminate.
c
      return
c------------- last line of dpchbs follows -----------------------------
      end
