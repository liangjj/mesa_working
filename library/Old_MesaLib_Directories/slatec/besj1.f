*deck besj1
      function besj1 (x)
c***begin prologue  besj1
c***purpose  compute the bessel function of the first kind of order one.
c***library   slatec (fnlib)
c***category  c10a1
c***type      single precision (besj1-s, dbesj1-d)
c***keywords  bessel function, first kind, fnlib, order one,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c besj1(x) calculates the bessel function of the first kind of
c order one for real argument x.
c
c series for bj1        on the interval  0.          to  1.60000d+01
c                                        with weighted error   4.48e-17
c                                         log weighted error  16.35
c                               significant figures required  15.77
c                                    decimal places required  16.89
c
c series for bm1        on the interval  0.          to  6.25000d-02
c                                        with weighted error   5.61e-17
c                                         log weighted error  16.25
c                               significant figures required  14.97
c                                    decimal places required  16.91
c
c series for bth1       on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.10e-17
c                                         log weighted error  16.39
c                               significant figures required  15.96
c                                    decimal places required  17.08
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   780601  date written
c   890210  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besj1
      dimension bj1cs(12), bm1cs(21), bth1cs(24)
      logical first
      save bj1cs, bm1cs, bth1cs, pi4, ntj1, ntm1, ntth1,
     1 xsml, xmin, xmax, first
      data bj1cs( 1) /   -.1172614151 3332787e0 /
      data bj1cs( 2) /   -.2536152183 0790640e0 /
      data bj1cs( 3) /    .0501270809 84469569e0 /
      data bj1cs( 4) /   -.0046315148 09625081e0 /
      data bj1cs( 5) /    .0002479962 29415914e0 /
      data bj1cs( 6) /   -.0000086789 48686278e0 /
      data bj1cs( 7) /    .0000002142 93917143e0 /
      data bj1cs( 8) /   -.0000000039 36093079e0 /
      data bj1cs( 9) /    .0000000000 55911823e0 /
      data bj1cs(10) /   -.0000000000 00632761e0 /
      data bj1cs(11) /    .0000000000 00005840e0 /
      data bj1cs(12) /   -.0000000000 00000044e0 /
      data bm1cs( 1) /    .1047362510 931285e0 /
      data bm1cs( 2) /    .0044244389 3702345e0 /
      data bm1cs( 3) /   -.0000566163 9504035e0 /
      data bm1cs( 4) /    .0000023134 9417339e0 /
      data bm1cs( 5) /   -.0000001737 7182007e0 /
      data bm1cs( 6) /    .0000000189 3209930e0 /
      data bm1cs( 7) /   -.0000000026 5416023e0 /
      data bm1cs( 8) /    .0000000004 4740209e0 /
      data bm1cs( 9) /   -.0000000000 8691795e0 /
      data bm1cs(10) /    .0000000000 1891492e0 /
      data bm1cs(11) /   -.0000000000 0451884e0 /
      data bm1cs(12) /    .0000000000 0116765e0 /
      data bm1cs(13) /   -.0000000000 0032265e0 /
      data bm1cs(14) /    .0000000000 0009450e0 /
      data bm1cs(15) /   -.0000000000 0002913e0 /
      data bm1cs(16) /    .0000000000 0000939e0 /
      data bm1cs(17) /   -.0000000000 0000315e0 /
      data bm1cs(18) /    .0000000000 0000109e0 /
      data bm1cs(19) /   -.0000000000 0000039e0 /
      data bm1cs(20) /    .0000000000 0000014e0 /
      data bm1cs(21) /   -.0000000000 0000005e0 /
      data bth1cs( 1) /    .7406014102 6313850e0 /
      data bth1cs( 2) /   -.0045717556 59637690e0 /
      data bth1cs( 3) /    .0001198185 10964326e0 /
      data bth1cs( 4) /   -.0000069645 61891648e0 /
      data bth1cs( 5) /    .0000006554 95621447e0 /
      data bth1cs( 6) /   -.0000000840 66228945e0 /
      data bth1cs( 7) /    .0000000133 76886564e0 /
      data bth1cs( 8) /   -.0000000024 99565654e0 /
      data bth1cs( 9) /    .0000000005 29495100e0 /
      data bth1cs(10) /   -.0000000001 24135944e0 /
      data bth1cs(11) /    .0000000000 31656485e0 /
      data bth1cs(12) /   -.0000000000 08668640e0 /
      data bth1cs(13) /    .0000000000 02523758e0 /
      data bth1cs(14) /   -.0000000000 00775085e0 /
      data bth1cs(15) /    .0000000000 00249527e0 /
      data bth1cs(16) /   -.0000000000 00083773e0 /
      data bth1cs(17) /    .0000000000 00029205e0 /
      data bth1cs(18) /   -.0000000000 00010534e0 /
      data bth1cs(19) /    .0000000000 00003919e0 /
      data bth1cs(20) /   -.0000000000 00001500e0 /
      data bth1cs(21) /    .0000000000 00000589e0 /
      data bth1cs(22) /   -.0000000000 00000237e0 /
      data bth1cs(23) /    .0000000000 00000097e0 /
      data bth1cs(24) /   -.0000000000 00000040e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data first /.true./
c***first executable statement  besj1
      if (first) then
         ntj1 = inits (bj1cs, 12, 0.1*r1mach(3))
         ntm1 = inits (bm1cs, 21, 0.1*r1mach(3))
         ntth1 = inits (bth1cs, 24, 0.1*r1mach(3))
c
         xsml = sqrt (8.0*r1mach(3))
         xmin = 2.0*r1mach(1)
         xmax = 1.0/r1mach(4)
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.4.0) go to 20
c
      besj1 = 0.
      if (y.eq.0.0) return
      if (y .le. xmin) call xermsg ('slatec', 'besj1',
     +   'abs(x) so small j1 underflows', 1, 1)
      if (y.gt.xmin) besj1 = 0.5*x
      if (y.gt.xsml) besj1 = x * (.25 + csevl(.125*y*y-1., bj1cs, ntj1))
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'besj1',
     +   'no precision because abs(x) is too big', 2, 2)
      z = 32.0/y**2 - 1.0
      ampl = (0.75 + csevl (z, bm1cs, ntm1)) / sqrt(y)
      theta = y - 3.0*pi4 + csevl (z, bth1cs, ntth1) / y
      besj1 = sign (ampl, x) * cos (theta)
c
      return
      end
