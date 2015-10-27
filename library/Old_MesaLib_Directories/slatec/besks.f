*deck besks
      subroutine besks (xnu, x, nin, bk)
c***begin prologue  besks
c***purpose  compute a sequence of modified bessel functions of the
c            third kind of fractional order.
c***library   slatec (fnlib)
c***category  c10b3
c***type      single precision (besks-s, dbesks-d)
c***keywords  fnlib, fractional order, modified bessel function,
c             sequence of bessel functions, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c besks computes a sequence of modified bessel functions of the third
c kind of order xnu + i at x, where x .gt. 0, xnu lies in (-1,1),
c and i = 0, 1, ... , nin - 1, if nin is positive and i = 0, 1, ... ,
c nin + 1, if nin is negative.  on return, the vector bk(.) contains
c the results at x for order starting at xnu.
c
c***references  (none)
c***routines called  beskes, r1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besks
      dimension bk(*)
      save xmax
      data xmax / 0.0 /
c***first executable statement  besks
      if (xmax.eq.0.0) xmax = -log (r1mach(1))
c
      if (x .gt. xmax) call xermsg ('slatec', 'besks',
     +   'x so big bessel k underflows', 1, 2)
c
      call beskes (xnu, x, nin, bk)
c
      expxi = exp (-x)
      n = abs (nin)
      do 20 i=1,n
        bk(i) = expxi * bk(i)
 20   continue
c
      return
      end
