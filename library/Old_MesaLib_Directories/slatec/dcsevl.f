*deck dcsevl
      double precision function dcsevl (x, cs, n)
c***begin prologue  dcsevl
c***purpose  evaluate a chebyshev series.
c***library   slatec (fnlib)
c***category  c3a2
c***type      double precision (csevl-s, dcsevl-d)
c***keywords  chebyshev series, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c  evaluate the n-term chebyshev series cs at x.  adapted from
c  a method presented in the paper by broucke referenced below.
c
c       input arguments --
c  x    value at which the series is to be evaluated.
c  cs   array of n terms of a chebyshev series.  in evaluating
c       cs, only half the first coefficient is summed.
c  n    number of terms in array cs.
c
c***references  r. broucke, ten subroutines for the manipulation of
c                 chebyshev series, algorithm 446, communications of
c                 the a.c.m. 16, (1973) pp. 254-256.
c               l. fox and i. b. parker, chebyshev polynomials in
c                 numerical analysis, oxford university press, 1968,
c                 page 56.
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900329  prologued revised extensively and code rewritten to allow
c           x to be slightly outside interval (-1,+1).  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dcsevl
      double precision b0, b1, b2, cs(*), onepl, twox, x, d1mach
      logical first
      save first, onepl
      data first /.true./
c***first executable statement  dcsevl
      if (first) onepl = 1.0d0 + d1mach(4)
      first = .false.
      if (n .lt. 1) call xermsg ('slatec', 'dcsevl',
     +   'number of terms .le. 0', 2, 2)
      if (n .gt. 1000) call xermsg ('slatec', 'dcsevl',
     +   'number of terms .gt. 1000', 3, 2)
      if (abs(x) .gt. onepl) call xermsg ('slatec', 'dcsevl',
     +   'x outside the interval (-1,+1)', 1, 1)
c
      b1 = 0.0d0
      b0 = 0.0d0
      twox = 2.0d0*x
      do 10 i = 1,n
         b2 = b1
         b1 = b0
         ni = n + 1 - i
         b0 = twox*b1 - b2 + cs(ni)
   10 continue
c
      dcsevl = 0.5d0*(b0-b2)
c
      return
      end
