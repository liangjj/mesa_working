*deck dbesks
      subroutine dbesks (xnu, x, nin, bk)
c***begin prologue  dbesks
c***purpose  compute a sequence of modified bessel functions of the
c            third kind of fractional order.
c***library   slatec (fnlib)
c***category  c10b3
c***type      double precision (besks-s, dbesks-d)
c***keywords  fnlib, fractional order, modified bessel function,
c             sequence of bessel functions, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c dbesks computes a sequence of modified bessel functions of the third
c kind of order xnu + i at x, where x .gt. 0, xnu lies in (-1,1),
c and i = 0, 1, ... , nin - 1, if nin is positive and i = 0, 1, ... ,
c nin + 1, if nin is negative.  on return, the vector bk(.) contains
c the results at x for order starting at xnu.  xnu, x, and bk are
c double precision.  nin is an integer.
c
c***references  (none)
c***routines called  d1mach, dbskes, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesks
      double precision xnu, x, bk(*), expxi, xmax, d1mach
      save xmax
      data xmax / 0.d0 /
c***first executable statement  dbesks
      if (xmax.eq.0.d0) xmax = -log (d1mach(1))
c
      if (x .gt. xmax) call xermsg ('slatec', 'dbesks',
     +   'x so big bessel k underflows', 1, 2)
c
      call dbskes (xnu, x, nin, bk)
c
      expxi = exp (-x)
      n = abs (nin)
      do 20 i=1,n
        bk(i) = expxi * bk(i)
 20   continue
c
      return
      end
