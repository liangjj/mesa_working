*deck alngam
      function alngam (x)
c***begin prologue  alngam
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      single precision (alngam-s, dlngam-d, clngam-c)
c***keywords  absolute value, complete gamma function, fnlib, logarithm,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c alngam(x) computes the logarithm of the absolute value of the
c gamma function at x.
c
c***references  (none)
c***routines called  gamma, r1mach, r9lgmc, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900727  added external statement.  (wrb)
c***end prologue  alngam
      logical first
      external gamma
      save sq2pil, sqpi2l, pi, xmax, dxrel, first
      data sq2pil / 0.9189385332 0467274e0/
      data sqpi2l / 0.2257913526 4472743e0/
      data pi     / 3.1415926535 8979324e0/
      data first  /.true./
c***first executable statement  alngam
      if (first) then
         xmax = r1mach(2)/log(r1mach(2))
         dxrel = sqrt (r1mach(4))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.10.0) go to 20
c
c log (abs (gamma(x))) for  abs(x) .le. 10.0
c
      alngam = log (abs (gamma(x)))
      return
c
c log (abs (gamma(x))) for abs(x) .gt. 10.0
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'alngam',
     +   'abs(x) so big alngam overflows', 2, 2)
c
      if (x.gt.0.) alngam = sq2pil + (x-0.5)*log(x) - x + r9lgmc(y)
      if (x.gt.0.) return
c
      sinpiy = abs (sin(pi*y))
      if (sinpiy .eq. 0.) call xermsg ('slatec', 'alngam',
     +   'x is a negative integer', 3, 2)
c
      if (abs((x-aint(x-0.5))/x) .lt. dxrel) call xermsg ('slatec',
     +   'alngam', 'answer lt half precision because x too near ' //
     +   'negative integer', 1, 1)
c
      alngam = sqpi2l + (x-0.5)*log(y) - x - log(sinpiy) - r9lgmc(y)
      return
c
      end
