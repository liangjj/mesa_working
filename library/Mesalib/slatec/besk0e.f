*deck besk0e
      function besk0e (x)
c***begin prologue  besk0e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the third kind of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besk0e-s, dbsk0e-d)
c***keywords  exponentially scaled, fnlib, hyperbolic bessel function,
c             modified bessel function, order zero, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c besk0e(x) computes the exponentially scaled modified (hyperbolic)
c bessel function of third kind of order zero for real argument
c x .gt. 0.0, i.e., exp(x)*k0(x).
c
c series for bk0        on the interval  0.          to  4.00000d+00
c                                        with weighted error   3.57e-19
c                                         log weighted error  18.45
c                               significant figures required  17.99
c                                    decimal places required  18.97
c
c series for ak0        on the interval  1.25000d-01 to  5.00000d-01
c                                        with weighted error   5.34e-17
c                                         log weighted error  16.27
c                               significant figures required  14.92
c                                    decimal places required  16.89
c
c series for ak02       on the interval  0.          to  1.25000d-01
c                                        with weighted error   2.34e-17
c                                         log weighted error  16.63
c                               significant figures required  14.67
c                                    decimal places required  17.20
c
c***references  (none)
c***routines called  besi0, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besk0e
      dimension bk0cs(11), ak0cs(17), ak02cs(14)
      logical first
      save bk0cs, ak0cs, ak02cs, ntk0, ntak0, ntak02, xsml, first
      data bk0cs( 1) /   -.0353273932 3390276872e0 /
      data bk0cs( 2) /    .3442898999 246284869e0 /
      data bk0cs( 3) /    .0359799365 1536150163e0 /
      data bk0cs( 4) /    .0012646154 1144692592e0 /
      data bk0cs( 5) /    .0000228621 2103119451e0 /
      data bk0cs( 6) /    .0000002534 7910790261e0 /
      data bk0cs( 7) /    .0000000019 0451637722e0 /
      data bk0cs( 8) /    .0000000000 1034969525e0 /
      data bk0cs( 9) /    .0000000000 0004259816e0 /
      data bk0cs(10) /    .0000000000 0000013744e0 /
      data bk0cs(11) /    .0000000000 0000000035e0 /
      data ak0cs( 1) /   -.0764394790 3327941e0 /
      data ak0cs( 2) /   -.0223565260 5699819e0 /
      data ak0cs( 3) /    .0007734181 1546938e0 /
      data ak0cs( 4) /   -.0000428100 6688886e0 /
      data ak0cs( 5) /    .0000030817 0017386e0 /
      data ak0cs( 6) /   -.0000002639 3672220e0 /
      data ak0cs( 7) /    .0000000256 3713036e0 /
      data ak0cs( 8) /   -.0000000027 4270554e0 /
      data ak0cs( 9) /    .0000000003 1694296e0 /
      data ak0cs(10) /   -.0000000000 3902353e0 /
      data ak0cs(11) /    .0000000000 0506804e0 /
      data ak0cs(12) /   -.0000000000 0068895e0 /
      data ak0cs(13) /    .0000000000 0009744e0 /
      data ak0cs(14) /   -.0000000000 0001427e0 /
      data ak0cs(15) /    .0000000000 0000215e0 /
      data ak0cs(16) /   -.0000000000 0000033e0 /
      data ak0cs(17) /    .0000000000 0000005e0 /
      data ak02cs( 1) /   -.0120186982 6307592e0 /
      data ak02cs( 2) /   -.0091748526 9102569e0 /
      data ak02cs( 3) /    .0001444550 9317750e0 /
      data ak02cs( 4) /   -.0000040136 1417543e0 /
      data ak02cs( 5) /    .0000001567 8318108e0 /
      data ak02cs( 6) /   -.0000000077 7011043e0 /
      data ak02cs( 7) /    .0000000004 6111825e0 /
      data ak02cs( 8) /   -.0000000000 3158592e0 /
      data ak02cs( 9) /    .0000000000 0243501e0 /
      data ak02cs(10) /   -.0000000000 0020743e0 /
      data ak02cs(11) /    .0000000000 0001925e0 /
      data ak02cs(12) /   -.0000000000 0000192e0 /
      data ak02cs(13) /    .0000000000 0000020e0 /
      data ak02cs(14) /   -.0000000000 0000002e0 /
      data first /.true./
c***first executable statement  besk0e
      if (first) then
         ntk0 = inits (bk0cs, 11, 0.1*r1mach(3))
         ntak0 = inits (ak0cs, 17, 0.1*r1mach(3))
         ntak02 = inits (ak02cs, 14, 0.1*r1mach(3))
         xsml = sqrt (4.0*r1mach(3))
      endif
      first = .false.
c
      if (x .le. 0.) call xermsg ('slatec', 'besk0e',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besk0e = exp(x) * (-log(0.5*x)*besi0(x)
     1  - .25 + csevl (.5*y-1., bk0cs, ntk0) )
      return
c
 20   if (x.le.8.) besk0e = (1.25 + csevl ((16./x-5.)/3., ak0cs, ntak0))
     1  / sqrt(x)
      if (x.gt.8.) besk0e = (1.25 + csevl (16./x-1., ak02cs, ntak02))
     1  / sqrt(x)
c
      return
      end
