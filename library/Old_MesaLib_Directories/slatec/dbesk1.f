*deck dbesk1
      double precision function dbesk1 (x)
c***begin prologue  dbesk1
c***purpose  compute the modified (hyperbolic) bessel function of the
c            third kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besk1-s, dbesk1-d)
c***keywords  fnlib, hyperbolic bessel function,
c             modified bessel function, order one, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c dbesk1(x) calculates the double precision modified (hyperbolic)
c bessel function of the third kind of order one for double precision
c argument x.  the argument must be large enough that the result does
c not overflow and small enough that the result does not underflow.
c
c series for bk1        on the interval  0.          to  4.00000e+00
c                                        with weighted error   9.16e-32
c                                         log weighted error  31.04
c                               significant figures required  30.61
c                                    decimal places required  31.64
c
c***references  (none)
c***routines called  d1mach, dbesi1, dbsk1e, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesk1
      double precision x, bk1cs(16), xmax, xmaxt, xmin, xsml, y,
     1  d1mach, dcsevl, dbesi1, dbsk1e
      logical first
      save bk1cs, ntk1, xmin, xsml, xmax, first
      data bk1cs(  1) / +.2530022733 8947770532 5311208685 33 d-1     /
      data bk1cs(  2) / -.3531559607 7654487566 7238316918 01 d+0     /
      data bk1cs(  3) / -.1226111808 2265714823 4790679300 42 d+0     /
      data bk1cs(  4) / -.6975723859 6398643501 8129202960 83 d-2     /
      data bk1cs(  5) / -.1730288957 5130520630 1765073689 79 d-3     /
      data bk1cs(  6) / -.2433406141 5659682349 6007350301 64 d-5     /
      data bk1cs(  7) / -.2213387630 7347258558 3152525451 26 d-7     /
      data bk1cs(  8) / -.1411488392 6335277610 9583302126 08 d-9     /
      data bk1cs(  9) / -.6666901694 1993290060 8537512643 73 d-12    /
      data bk1cs( 10) / -.2427449850 5193659339 2631968648 53 d-14    /
      data bk1cs( 11) / -.7023863479 3862875971 7837971200 00 d-17    /
      data bk1cs( 12) / -.1654327515 5100994675 4910293333 33 d-19    /
      data bk1cs( 13) / -.3233834745 9944491991 8933333333 33 d-22    /
      data bk1cs( 14) / -.5331275052 9265274999 4666666666 66 d-25    /
      data bk1cs( 15) / -.7513040716 2157226666 6666666666 66 d-28    /
      data bk1cs( 16) / -.9155085717 6541866666 6666666666 66 d-31    /
      data first /.true./
c***first executable statement  dbesk1
      if (first) then
         ntk1 = initds (bk1cs, 16, 0.1*real(d1mach(3)))
         xmin = exp(max(log(d1mach(1)), -log(d1mach(2))) + 0.01d0)
         xsml = sqrt(4.0d0*d1mach(3))
         xmaxt = -log(d1mach(1))
         xmax = xmaxt - 0.5d0*xmaxt*log(xmaxt)/(xmaxt+0.5d0)
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'dbesk1',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.0d0) go to 20
c
      if (x .lt. xmin) call xermsg ('slatec', 'dbesk1',
     +   'x so small k1 overflows', 3, 2)
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesk1 = log(0.5d0*x)*dbesi1(x) + (0.75d0 + dcsevl (.5d0*y-1.d0,
     1  bk1cs, ntk1))/x
      return
c
 20   dbesk1 = 0.d0
      if (x .gt. xmax) call xermsg ('slatec', 'dbesk1',
     +   'x so big k1 underflows', 1, 1)
      if (x.gt.xmax) return
c
      dbesk1 = exp(-x) * dbsk1e(x)
c
      return
      end
