*deck @(#)dcsevl.f	5.1  4/18/95
      function dcsevl (x, a, n)
      real*8 dcsevl
c
c evaluate the n-term chebyshev series a at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c             input arguments --
c x      dble prec value at which the series is to be evaluated.
c a      dble prec array of n terms of a chebyshev series.  in eval-
c        uating a, only half the first coef is summed.
c n      number of terms in array a.
c
      real*8 a(n), x, twox, b0, b1, b2
c
      if (n.lt.1) call lnkerr('dcsevl: number of terms le 0')
      if (n.gt.1000) call lnkerr('dcsevl: number of terms gt 1000')
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call lnkerr (
     1  'dcsevl: x outside (-1,+1)')
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end
