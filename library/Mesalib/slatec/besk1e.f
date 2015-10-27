*deck besk1e
      function besk1e (x)
c***begin prologue  besk1e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the third kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besk1e-s, dbsk1e-d)
c***keywords  exponentially scaled, fnlib, hyperbolic bessel function,
c             modified bessel function, order one, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c besk1e(x) computes the exponentially scaled modified (hyperbolic)
c bessel function of third kind of order one for real argument
c x .gt. 0.0, i.e., exp(x)*k1(x).
c
c series for bk1        on the interval  0.          to  4.00000d+00
c                                        with weighted error   7.02e-18
c                                         log weighted error  17.15
c                               significant figures required  16.73
c                                    decimal places required  17.67
c
c series for ak1        on the interval  1.25000d-01 to  5.00000d-01
c                                        with weighted error   6.06e-17
c                                         log weighted error  16.22
c                               significant figures required  15.41
c                                    decimal places required  16.83
c
c series for ak12       on the interval  0.          to  1.25000d-01
c                                        with weighted error   2.58e-17
c                                         log weighted error  16.59
c                               significant figures required  15.22
c                                    decimal places required  17.16
c
c***references  (none)
c***routines called  besi1, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besk1e
      dimension bk1cs(11), ak1cs(17), ak12cs(14)
      logical first
      save bk1cs, ak1cs, ak12cs, ntk1, ntak1, ntak12, xmin, xsml,
     1 first
      data bk1cs( 1) /    .0253002273 389477705e0 /
      data bk1cs( 2) /   -.3531559607 76544876e0 /
      data bk1cs( 3) /   -.1226111808 22657148e0 /
      data bk1cs( 4) /   -.0069757238 596398643e0 /
      data bk1cs( 5) /   -.0001730288 957513052e0 /
      data bk1cs( 6) /   -.0000024334 061415659e0 /
      data bk1cs( 7) /   -.0000000221 338763073e0 /
      data bk1cs( 8) /   -.0000000001 411488392e0 /
      data bk1cs( 9) /   -.0000000000 006666901e0 /
      data bk1cs(10) /   -.0000000000 000024274e0 /
      data bk1cs(11) /   -.0000000000 000000070e0 /
      data ak1cs( 1) /    .2744313406 973883e0 /
      data ak1cs( 2) /    .0757198995 3199368e0 /
      data ak1cs( 3) /   -.0014410515 5647540e0 /
      data ak1cs( 4) /    .0000665011 6955125e0 /
      data ak1cs( 5) /   -.0000043699 8470952e0 /
      data ak1cs( 6) /    .0000003540 2774997e0 /
      data ak1cs( 7) /   -.0000000331 1163779e0 /
      data ak1cs( 8) /    .0000000034 4597758e0 /
      data ak1cs( 9) /   -.0000000003 8989323e0 /
      data ak1cs(10) /    .0000000000 4720819e0 /
      data ak1cs(11) /   -.0000000000 0604783e0 /
      data ak1cs(12) /    .0000000000 0081284e0 /
      data ak1cs(13) /   -.0000000000 0011386e0 /
      data ak1cs(14) /    .0000000000 0001654e0 /
      data ak1cs(15) /   -.0000000000 0000248e0 /
      data ak1cs(16) /    .0000000000 0000038e0 /
      data ak1cs(17) /   -.0000000000 0000006e0 /
      data ak12cs( 1) /    .0637930834 3739001e0 /
      data ak12cs( 2) /    .0283288781 3049721e0 /
      data ak12cs( 3) /   -.0002475370 6739052e0 /
      data ak12cs( 4) /    .0000057719 7245160e0 /
      data ak12cs( 5) /   -.0000002068 9392195e0 /
      data ak12cs( 6) /    .0000000097 3998344e0 /
      data ak12cs( 7) /   -.0000000005 5853361e0 /
      data ak12cs( 8) /    .0000000000 3732996e0 /
      data ak12cs( 9) /   -.0000000000 0282505e0 /
      data ak12cs(10) /    .0000000000 0023720e0 /
      data ak12cs(11) /   -.0000000000 0002176e0 /
      data ak12cs(12) /    .0000000000 0000215e0 /
      data ak12cs(13) /   -.0000000000 0000022e0 /
      data ak12cs(14) /    .0000000000 0000002e0 /
      data first /.true./
c***first executable statement  besk1e
      if (first) then
         ntk1 = inits (bk1cs, 11, 0.1*r1mach(3))
         ntak1 = inits (ak1cs, 17, 0.1*r1mach(3))
         ntak12 = inits (ak12cs, 14, 0.1*r1mach(3))
c
         xmin = exp (max(log(r1mach(1)), -log(r1mach(2))) + .01)
         xsml = sqrt (4.0*r1mach(3))
      endif
      first = .false.
c
      if (x .le. 0.) call xermsg ('slatec', 'besk1e',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.0) go to 20
c
      if (x .lt. xmin) call xermsg ('slatec', 'besk1e',
     +   'x so small k1 overflows', 3, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besk1e = exp(x) * (log(0.5*x)*besi1(x) +
     1  (0.75 + csevl (.5*y-1., bk1cs, ntk1))/x )
      return
c
 20   if (x.le.8.) besk1e = (1.25 + csevl ((16./x-5.)/3., ak1cs, ntak1))
     1  / sqrt(x)
      if (x.gt.8.) besk1e = (1.25 + csevl (16./x-1., ak12cs, ntak12))
     1  / sqrt(x)
c
      return
      end
