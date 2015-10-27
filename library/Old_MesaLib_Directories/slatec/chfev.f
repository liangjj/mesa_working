*deck chfev
      subroutine chfev (x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr)
c***begin prologue  chfev
c***purpose  evaluate a cubic polynomial given in hermite form at an
c            array of points.  while designed for use by pchfe, it may
c            be useful directly as an evaluator for a piecewise cubic
c            hermite function in applications, such as graphing, where
c            the interval is known in advance.
c***library   slatec (pchip)
c***category  e3
c***type      single precision (chfev-s, dchfev-d)
c***keywords  cubic hermite evaluation, cubic polynomial evaluation,
c             pchip
c***author  fritsch, f. n., (llnl)
c             lawrence livermore national laboratory
c             p.o. box 808  (l-316)
c             livermore, ca  94550
c             fts 532-4275, (510) 422-4275
c***description
c
c          chfev:  cubic hermite function evaluator
c
c     evaluates the cubic polynomial determined by function values
c     f1,f2 and derivatives d1,d2 on interval (x1,x2) at the points
c     xe(j), j=1(1)ne.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        integer  ne, next(2), ierr
c        real  x1, x2, f1, f2, d1, d2, xe(ne), fe(ne)
c
c        call  chfev (x1,x2, f1,f2, d1,d2, ne, xe, fe, next, ierr)
c
c   parameters:
c
c     x1,x2 -- (input) endpoints of interval of definition of cubic.
c           (error return if  x1.eq.x2 .)
c
c     f1,f2 -- (input) values of function at x1 and x2, respectively.
c
c     d1,d2 -- (input) values of derivative at x1 and x2, respectively.
c
c     ne -- (input) number of evaluation points.  (error return if
c           ne.lt.1 .)
c
c     xe -- (input) real array of points at which the function is to be
c           evaluated.  if any of the xe are outside the interval
c           [x1,x2], a warning error is returned in next.
c
c     fe -- (output) real array of values of the cubic function defined
c           by  x1,x2, f1,f2, d1,d2  at the points  xe.
c
c     next -- (output) integer array indicating number of extrapolation
c           points:
c            next(1) = number of evaluation points to left of interval.
c            next(2) = number of evaluation points to right of interval.
c
c     ierr -- (output) error flag.
c           normal return:
c              ierr = 0  (no errors).
c           "recoverable" errors:
c              ierr = -1  if ne.lt.1 .
c              ierr = -2  if x1.eq.x2 .
c                (the fe-array has not been changed in either case.)
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   811019  date written
c   820803  minor cosmetic changes for release 1.
c   890411  added save statements (vers. 3.2).
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890703  corrected category record.  (wrb)
c   890703  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  chfev
c  programming notes:
c
c     to produce a double precision version, simply:
c        a. change chfev to dchfev wherever it occurs,
c        b. change the real declaration to double precision, and
c        c. change the constant zero to double precision.
c
c  declare arguments.
c
      integer  ne, next(2), ierr
      real  x1, x2, f1, f2, d1, d2, xe(*), fe(*)
c
c  declare local variables.
c
      integer  i
      real  c2, c3, del1, del2, delta, h, x, xmi, xma, zero
      save zero
      data  zero /0./
c
c  validity-check arguments.
c
c***first executable statement  chfev
      if (ne .lt. 1)  go to 5001
      h = x2 - x1
      if (h .eq. zero)  go to 5002
c
c  initialize.
c
      ierr = 0
      next(1) = 0
      next(2) = 0
      xmi = min(zero, h)
      xma = max(zero, h)
c
c  compute cubic coefficients (expanded about x1).
c
      delta = (f2 - f1)/h
      del1 = (d1 - delta)/h
      del2 = (d2 - delta)/h
c                                           (delta is no longer needed.)
      c2 = -(del1+del1 + del2)
      c3 = (del1 + del2)/h
c                               (h, del1 and del2 are no longer needed.)
c
c  evaluation loop.
c
      do 500  i = 1, ne
         x = xe(i) - x1
         fe(i) = f1 + x*(d1 + x*(c2 + x*c3))
c          count extrapolation points.
         if ( x.lt.xmi )  next(1) = next(1) + 1
         if ( x.gt.xma )  next(2) = next(2) + 1
c        (note redundancy--if either condition is true, other is false.)
  500 continue
c
c  normal return.
c
      return
c
c  error returns.
c
 5001 continue
c     ne.lt.1 return.
      ierr = -1
      call xermsg ('slatec', 'chfev',
     +   'number of evaluation points less than one', ierr, 1)
      return
c
 5002 continue
c     x1.eq.x2 return.
      ierr = -2
      call xermsg ('slatec', 'chfev', 'interval endpoints equal', ierr,
     +   1)
      return
c------------- last line of chfev follows ------------------------------
      end
