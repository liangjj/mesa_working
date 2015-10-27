*deck fac
      function fac (n)
c***begin prologue  fac
c***purpose  compute the factorial function.
c***library   slatec (fnlib)
c***category  c1
c***type      single precision (fac-s, dfac-d)
c***keywords  factorial, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c fac(n) evaluates the factorial function of n.  fac is single
c precision.  n must be an integer between 0 and 25 inclusive.
c
c***references  (none)
c***routines called  gamlim, r9lgmc, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  fac
      dimension facn(26)
      save facn, sq2pil, nmax
      data facn( 1) / 1.0e0 /
      data facn( 2) / 1.0e0 /
      data facn( 3) / 2.0e0 /
      data facn( 4) / 6.0e0 /
      data facn( 5) / 24.0e0 /
      data facn( 6) / 120.0e0 /
      data facn( 7) / 720.0e0 /
      data facn( 8) / 5040.0e0 /
      data facn( 9) / 40320.0e0 /
      data facn(10) / 362880.0e0 /
      data facn(11) / 3628800.0e0 /
      data facn(12) / 39916800.0e0 /
      data facn(13) / 479001600.0e0 /
      data facn(14) / 6227020800.0e0 /
      data facn(15) / 87178291200.0e0 /
      data facn(16) / 1307674368000.0e0 /
      data facn(17) / 20922789888000.0e0 /
      data facn(18) / 355687428096000.0e0 /
      data facn(19) / 6402373705728000.0e0 /
      data facn(20) /  .12164510040883200e18 /
      data facn(21) /  .24329020081766400e19 /
      data facn(22) /  .51090942171709440e20 /
      data facn(23) /  .11240007277776077e22 /
      data facn(24) /  .25852016738884977e23 /
      data facn(25) /  .62044840173323944e24 /
      data facn(26) /  .15511210043330986e26 /
      data sq2pil / 0.9189385332 0467274e0/
      data nmax / 0 /
c***first executable statement  fac
      if (nmax.ne.0) go to 10
      call gamlim (xmin, xmax)
      nmax = xmax - 1.
c
 10   if (n .lt. 0) call xermsg ('slatec', 'fac',
     +   'factorial of negative integer undefined', 1, 2)
c
      if (n.le.25) fac = facn(n+1)
      if (n.le.25) return
c
      if (n .gt. nmax) call xermsg ('slatec', 'fac',
     +   'n so big factorial(n) overflows', 2, 2)
c
      x = n + 1
      fac = exp ( (x-0.5)*log(x) - x + sq2pil + r9lgmc(x) )
c
      return
      end
