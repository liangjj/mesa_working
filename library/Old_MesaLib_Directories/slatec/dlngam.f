*deck dlngam
      double precision function dlngam (x)
c***begin prologue  dlngam
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (alngam-s, dlngam-d, clngam-c)
c***keywords  absolute value, complete gamma function, fnlib, logarithm,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dlngam(x) calculates the double precision logarithm of the
c absolute value of the gamma function for double precision
c argument x.
c
c***references  (none)
c***routines called  d1mach, d9lgmc, dgamma, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dlngam
      double precision x, dxrel, pi, sinpiy, sqpi2l, sq2pil, xmax,
     1  y, dgamma, d9lgmc, d1mach, temp
      logical first
      external dgamma
      save sq2pil, sqpi2l, pi, xmax, dxrel, first
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data sqpi2l / +.2257913526 4472743236 3097614947 441 d+0    /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data first /.true./
c***first executable statement  dlngam
      if (first) then
         temp = 1.d0/log(d1mach(2))
         xmax = temp*d1mach(2)
         dxrel = sqrt(d1mach(4))
      endif
      first = .false.
c
      y = abs (x)
      if (y.gt.10.d0) go to 20
c
c log (abs (dgamma(x)) ) for abs(x) .le. 10.0
c
      dlngam = log (abs (dgamma(x)) )
      return
c
c log ( abs (dgamma(x)) ) for abs(x) .gt. 10.0
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'dlngam',
     +   'abs(x) so big dlngam overflows', 2, 2)
c
      if (x.gt.0.d0) dlngam = sq2pil + (x-0.5d0)*log(x) - x + d9lgmc(y)
      if (x.gt.0.d0) return
c
      sinpiy = abs (sin(pi*y))
      if (sinpiy .eq. 0.d0) call xermsg ('slatec', 'dlngam',
     +   'x is a negative integer', 3, 2)
c
      if (abs((x-aint(x-0.5d0))/x) .lt. dxrel) call xermsg ('slatec',
     +   'dlngam',
     +   'answer lt half precision because x too near negative integer',
     +   1, 1)
c
      dlngam = sqpi2l + (x-0.5d0)*log(y) - x - log(sinpiy) - d9lgmc(y)
      return
c
      end
