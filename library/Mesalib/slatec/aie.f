*deck aie
      function aie (x)
c***begin prologue  aie
c***purpose  calculate the airy function for a negative argument and an
c            exponentially scaled airy function for a non-negative
c            argument.
c***library   slatec (fnlib)
c***category  c10d
c***type      single precision (aie-s, daie-d)
c***keywords  exponentially scaled airy function, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c aie(x) computes the exponentially scaled airy function for
c non-negative x.  it evaluates ai(x) for x .le. 0.0 and
c exp(zeta)*ai(x) for x .ge. 0.0 where zeta = (2.0/3.0)*(x**1.5).
c
c series for aif        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   1.09e-19
c                                         log weighted error  18.96
c                               significant figures required  17.76
c                                    decimal places required  19.44
c
c series for aig        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   1.51e-17
c                                         log weighted error  16.82
c                               significant figures required  15.19
c                                    decimal places required  17.27
c
c series for aip        on the interval  0.          to  1.00000d+00
c                                        with weighted error   5.10e-17
c                                         log weighted error  16.29
c                               significant figures required  14.41
c                                    decimal places required  17.06
c
c***references  (none)
c***routines called  csevl, inits, r1mach, r9aimp
c***revision history  (yymmdd)
c   770701  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  aie
      dimension aifcs(9), aigcs(8), aipcs(34)
      logical first
      save aifcs, aigcs, aipcs, naif, naig,
     1 naip, x3sml, x32sml, xbig, first
      data aifcs( 1) /   -.0379713584 9666999750e0 /
      data aifcs( 2) /    .0591918885 3726363857e0 /
      data aifcs( 3) /    .0009862928 0577279975e0 /
      data aifcs( 4) /    .0000068488 4381907656e0 /
      data aifcs( 5) /    .0000000259 4202596219e0 /
      data aifcs( 6) /    .0000000000 6176612774e0 /
      data aifcs( 7) /    .0000000000 0010092454e0 /
      data aifcs( 8) /    .0000000000 0000012014e0 /
      data aifcs( 9) /    .0000000000 0000000010e0 /
      data aigcs( 1) /    .0181523655 8116127e0 /
      data aigcs( 2) /    .0215725631 6601076e0 /
      data aigcs( 3) /    .0002567835 6987483e0 /
      data aigcs( 4) /    .0000014265 2141197e0 /
      data aigcs( 5) /    .0000000045 7211492e0 /
      data aigcs( 6) /    .0000000000 0952517e0 /
      data aigcs( 7) /    .0000000000 0001392e0 /
      data aigcs( 8) /    .0000000000 0000001e0 /
      data aipcs( 1) /   -.0187519297 793868e0 /
      data aipcs( 2) /   -.0091443848 250055e0 /
      data aipcs( 3) /    .0009010457 337825e0 /
      data aipcs( 4) /   -.0001394184 127221e0 /
      data aipcs( 5) /    .0000273815 815785e0 /
      data aipcs( 6) /   -.0000062750 421119e0 /
      data aipcs( 7) /    .0000016064 844184e0 /
      data aipcs( 8) /   -.0000004476 392158e0 /
      data aipcs( 9) /    .0000001334 635874e0 /
      data aipcs(10) /   -.0000000420 735334e0 /
      data aipcs(11) /    .0000000139 021990e0 /
      data aipcs(12) /   -.0000000047 831848e0 /
      data aipcs(13) /    .0000000017 047897e0 /
      data aipcs(14) /   -.0000000006 268389e0 /
      data aipcs(15) /    .0000000002 369824e0 /
      data aipcs(16) /   -.0000000000 918641e0 /
      data aipcs(17) /    .0000000000 364278e0 /
      data aipcs(18) /   -.0000000000 147475e0 /
      data aipcs(19) /    .0000000000 060851e0 /
      data aipcs(20) /   -.0000000000 025552e0 /
      data aipcs(21) /    .0000000000 010906e0 /
      data aipcs(22) /   -.0000000000 004725e0 /
      data aipcs(23) /    .0000000000 002076e0 /
      data aipcs(24) /   -.0000000000 000924e0 /
      data aipcs(25) /    .0000000000 000417e0 /
      data aipcs(26) /   -.0000000000 000190e0 /
      data aipcs(27) /    .0000000000 000087e0 /
      data aipcs(28) /   -.0000000000 000040e0 /
      data aipcs(29) /    .0000000000 000019e0 /
      data aipcs(30) /   -.0000000000 000009e0 /
      data aipcs(31) /    .0000000000 000004e0 /
      data aipcs(32) /   -.0000000000 000002e0 /
      data aipcs(33) /    .0000000000 000001e0 /
      data aipcs(34) /   -.0000000000 000000e0 /
      data first /.true./
c***first executable statement  aie
      if (first) then
         eta = 0.1*r1mach(3)
         naif  = inits (aifcs, 9, eta)
         naig  = inits (aigcs, 8, eta)
         naip  = inits (aipcs, 34, eta)
c
         x3sml = eta**0.3333
         x32sml = 1.3104*x3sml**2
         xbig = r1mach(2)**0.6666
      endif
      first = .false.
c
      if (x.ge.(-1.0)) go to 20
      call r9aimp (x, xm, theta)
      aie = xm * cos(theta)
      return
c
 20   if (x.gt.1.0) go to 30
      z = 0.0
      if (abs(x).gt.x3sml) z = x**3
      aie = 0.375 + (csevl (z, aifcs, naif) - x*(0.25 +
     1  csevl (z, aigcs, naig)) )
      if (x.gt.x32sml) aie = aie * exp(2.0*x*sqrt(x)/3.0)
      return
c
 30   sqrtx = sqrt(x)
      z = -1.0
      if (x.lt.xbig) z = 2.0/(x*sqrtx) - 1.0
      aie = (.28125 + csevl (z, aipcs, naip))/sqrt(sqrtx)
      return
c
      end
