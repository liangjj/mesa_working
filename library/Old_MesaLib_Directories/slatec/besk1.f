*deck besk1
      function besk1 (x)
c***begin prologue  besk1
c***purpose  compute the modified (hyperbolic) bessel function of the
c            third kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besk1-s, dbesk1-d)
c***keywords  fnlib, hyperbolic bessel function,
c             modified bessel function, order one, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c besk1(x) computes the modified (hyperbolic) bessel function of third
c kind of order one for real argument x, where x .gt. 0.
c
c series for bk1        on the interval  0.          to  4.00000d+00
c                                        with weighted error   7.02e-18
c                                         log weighted error  17.15
c                               significant figures required  16.73
c                                    decimal places required  17.67
c
c***references  (none)
c***routines called  besi1, besk1e, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besk1
      dimension bk1cs(11)
      logical first
      save bk1cs, ntk1, xmin, xsml, xmax, first
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
      data first /.true./
c***first executable statement  besk1
      if (first) then
         ntk1 = inits (bk1cs, 11, 0.1*r1mach(3))
         xmin = exp (max(log(r1mach(1)), -log(r1mach(2))) + .01)
         xsml = sqrt (4.0*r1mach(3))
         xmaxt = -log(r1mach(1))
         xmax = xmaxt - 0.5*xmaxt*log(xmaxt)/(xmaxt+0.5)
      endif
      first = .false.
c
      if (x .le. 0.) call xermsg ('slatec', 'besk1',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.0) go to 20
c
      if (x .lt. xmin) call xermsg ('slatec', 'besk1',
     +   'x so small k1 overflows', 3, 2)
      y = 0.
      if (x.gt.xsml) y = x*x
      besk1 = log(0.5*x)*besi1(x) +
     1  (0.75 + csevl (.5*y-1., bk1cs, ntk1))/x
      return
c
 20   besk1 = 0.
      if (x .gt. xmax) call xermsg ('slatec', 'besk1',
     +   'x so big k1 underflows', 1, 1)
      if (x.gt.xmax) return
c
      besk1 = exp(-x) * besk1e(x)
c
      return
      end
