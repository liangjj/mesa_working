*deck besy1
      function besy1 (x)
c***begin prologue  besy1
c***purpose  compute the bessel function of the second kind of order
c            one.
c***library   slatec (fnlib)
c***category  c10a1
c***type      single precision (besy1-s, dbesy1-d)
c***keywords  bessel function, fnlib, order one, second kind,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c besy1(x) calculates the bessel function of the second kind of
c order one for real argument x.
c
c series for by1        on the interval  0.          to  1.60000d+01
c                                        with weighted error   1.87e-18
c                                         log weighted error  17.73
c                               significant figures required  17.83
c                                    decimal places required  18.30
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
c***routines called  besj1, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besy1
      dimension by1cs(14), bm1cs(21), bth1cs(24)
      logical first
      save by1cs, bm1cs, bth1cs, twodpi, pi4,
     1 nty1, ntm1, ntth1, xmin, xsml, xmax, first
      data by1cs( 1) /    .0320804710 0611908629e0 /
      data by1cs( 2) /   1.2627078974 33500450e0 /
      data by1cs( 3) /    .0064999618 9992317500e0 /
      data by1cs( 4) /   -.0893616452 8860504117e0 /
      data by1cs( 5) /    .0132508812 2175709545e0 /
      data by1cs( 6) /   -.0008979059 1196483523e0 /
      data by1cs( 7) /    .0000364736 1487958306e0 /
      data by1cs( 8) /   -.0000010013 7438166600e0 /
      data by1cs( 9) /    .0000000199 4539657390e0 /
      data by1cs(10) /   -.0000000003 0230656018e0 /
      data by1cs(11) /    .0000000000 0360987815e0 /
      data by1cs(12) /   -.0000000000 0003487488e0 /
      data by1cs(13) /    .0000000000 0000027838e0 /
      data by1cs(14) /   -.0000000000 0000000186e0 /
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
      data twodpi / 0.6366197723 6758134e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data first /.true./
c***first executable statement  besy1
      if (first) then
         nty1 = inits (by1cs, 14, 0.1*r1mach(3))
         ntm1 = inits (bm1cs, 21, 0.1*r1mach(3))
         ntth1 = inits (bth1cs, 24, 0.1*r1mach(3))
c
         xmin = 1.571*exp ( max(log(r1mach(1)), -log(r1mach(2)))+.01)
         xsml = sqrt (4.0*r1mach(3))
         xmax = 1.0/r1mach(4)
      endif
      first = .false.
c
      if (x .le. 0.) call xermsg ('slatec', 'besy1',
     +   'x is zero or negative', 1, 2)
      if (x.gt.4.0) go to 20
c
      if (x .lt. xmin) call xermsg ('slatec', 'besy1',
     +   'x so small y1 overflows', 3, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besy1 = twodpi*log(0.5*x)*besj1(x) +
     1  (0.5 + csevl (.125*y-1., by1cs, nty1))/x
      return
c
 20   if (x .gt. xmax) call xermsg ('slatec', 'besy1',
     +   'no precision because x is big', 2, 2)
c
      z = 32.0/x**2 - 1.0
      ampl = (0.75 + csevl (z, bm1cs, ntm1)) / sqrt(x)
      theta = x - 3.0*pi4 + csevl (z, bth1cs, ntth1) / x
      besy1 = ampl * sin (theta)
c
      return
      end
