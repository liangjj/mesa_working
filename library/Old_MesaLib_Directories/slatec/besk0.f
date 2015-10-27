*deck besk0
      function besk0 (x)
c***begin prologue  besk0
c***purpose  compute the modified (hyperbolic) bessel function of the
c            third kind of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besk0-s, dbesk0-d)
c***keywords  fnlib, hyperbolic bessel function,
c             modified bessel function, order zero, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c besk0(x) calculates the modified (hyperbolic) bessel function
c of the third kind of order zero for real argument x .gt. 0.0.
c
c series for bk0        on the interval  0.          to  4.00000d+00
c                                        with weighted error   3.57e-19
c                                         log weighted error  18.45
c                               significant figures required  17.99
c                                    decimal places required  18.97
c
c***references  (none)
c***routines called  besi0, besk0e, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besk0
      dimension bk0cs(11)
      logical first
      save bk0cs, ntk0, xsml, xmax, first
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
      data first /.true./
c***first executable statement  besk0
      if (first) then
         ntk0 = inits (bk0cs, 11, 0.1*r1mach(3))
         xsml = sqrt (4.0*r1mach(3))
         xmaxt = -log(r1mach(1))
         xmax = xmaxt - 0.5*xmaxt*log(xmaxt)/(xmaxt+0.5) - 0.01
      endif
      first = .false.
c
      if (x .le. 0.) call xermsg ('slatec', 'besk0',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besk0 = -log(0.5*x)*besi0(x) - .25 + csevl (.5*y-1., bk0cs, ntk0)
      return
c
 20   besk0 = 0.
      if (x .gt. xmax) call xermsg ('slatec', 'besk0',
     +   'x so big k0 underflows', 1, 1)
      if (x.gt.xmax) return
c
      besk0 = exp(-x) * besk0e(x)
c
      return
      end
