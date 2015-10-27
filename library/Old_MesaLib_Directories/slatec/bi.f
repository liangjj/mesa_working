*deck bi
      function bi (x)
c***begin prologue  bi
c***purpose  evaluate the bairy function (the airy function of the
c            second kind).
c***library   slatec (fnlib)
c***category  c10d
c***type      single precision (bi-s, dbi-d)
c***keywords  bairy function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c bi(x) calculates the airy function of the second kind for real
c argument x.
c
c series for bif        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   1.88e-19
c                                         log weighted error  18.72
c                               significant figures required  17.74
c                                    decimal places required  19.20
c
c series for big        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   2.61e-17
c                                         log weighted error  16.58
c                               significant figures required  15.17
c                                    decimal places required  17.03
c
c series for bif2       on the interval  1.00000d+00 to  8.00000d+00
c                                        with weighted error   1.11e-17
c                                         log weighted error  16.95
c                        approx significant figures required  16.5
c                                    decimal places required  17.45
c
c series for big2       on the interval  1.00000d+00 to  8.00000d+00
c                                        with weighted error   1.19e-18
c                                         log weighted error  17.92
c                        approx significant figures required  17.2
c                                    decimal places required  18.42
c
c***references  (none)
c***routines called  bie, csevl, inits, r1mach, r9aimp, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  bi
      dimension bifcs(9), bigcs(8), bif2cs(10), big2cs(10)
      logical first
      save bifcs, bigcs, bif2cs, big2cs, nbif, nbig, nbif2,
     1 nbig2, x3sml, xmax, first
      data bifcs( 1) /   -.0167302164 7198664948e0 /
      data bifcs( 2) /    .1025233583 424944561e0 /
      data bifcs( 3) /    .0017083092 5073815165e0 /
      data bifcs( 4) /    .0000118625 4546774468e0 /
      data bifcs( 5) /    .0000000449 3290701779e0 /
      data bifcs( 6) /    .0000000001 0698207143e0 /
      data bifcs( 7) /    .0000000000 0017480643e0 /
      data bifcs( 8) /    .0000000000 0000020810e0 /
      data bifcs( 9) /    .0000000000 0000000018e0 /
      data bigcs( 1) /    .0224662232 4857452e0 /
      data bigcs( 2) /    .0373647754 5301955e0 /
      data bigcs( 3) /    .0004447621 8957212e0 /
      data bigcs( 4) /    .0000024708 0756363e0 /
      data bigcs( 5) /    .0000000079 1913533e0 /
      data bigcs( 6) /    .0000000000 1649807e0 /
      data bigcs( 7) /    .0000000000 0002411e0 /
      data bigcs( 8) /    .0000000000 0000002e0 /
      data bif2cs( 1) /   0.0998457269 3816041e0 /
      data bif2cs( 2) /    .4786249778 63005538e0 /
      data bif2cs( 3) /    .0251552119 604330118e0 /
      data bif2cs( 4) /    .0005820693 885232645e0 /
      data bif2cs( 5) /    .0000074997 659644377e0 /
      data bif2cs( 6) /    .0000000613 460287034e0 /
      data bif2cs( 7) /    .0000000003 462753885e0 /
      data bif2cs( 8) /    .0000000000 014288910e0 /
      data bif2cs( 9) /    .0000000000 000044962e0 /
      data bif2cs(10) /    .0000000000 000000111e0 /
      data big2cs( 1) /    .0333056621 45514340e0 /
      data big2cs( 2) /    .1613092151 23197068e0 /
      data big2cs( 3) /    .0063190073 096134286e0 /
      data big2cs( 4) /    .0001187904 568162517e0 /
      data big2cs( 5) /    .0000013045 345886200e0 /
      data big2cs( 6) /    .0000000093 741259955e0 /
      data big2cs( 7) /    .0000000000 474580188e0 /
      data big2cs( 8) /    .0000000000 001783107e0 /
      data big2cs( 9) /    .0000000000 000005167e0 /
      data big2cs(10) /    .0000000000 000000011e0 /
      data first /.true./
c***first executable statement  bi
      if (first) then
         eta = 0.1*r1mach(3)
         nbif  = inits (bifcs , 9, eta)
         nbig  = inits (bigcs , 8, eta)
         nbif2 = inits (bif2cs, 10, eta)
         nbig2 = inits (big2cs, 10, eta)
c
         x3sml = eta**0.3333
         xmax = (1.5*log(r1mach(2)))**0.6666
      endif
      first = .false.
c
      if (x.ge.(-1.0)) go to 20
      call r9aimp (x, xm, theta)
      bi = xm * sin(theta)
      return
c
 20   if (x.gt.1.0) go to 30
      z = 0.0
      if (abs(x).gt.x3sml) z = x**3
      bi = 0.625 + csevl (z, bifcs, nbif) + x*(0.4375 +
     1  csevl (z, bigcs, nbig))
      return
c
 30   if (x.gt.2.0) go to 40
      z = (2.0*x**3 - 9.0) / 7.0
      bi = 1.125 + csevl (z, bif2cs, nbif2) + x*(0.625 +
     1  csevl (z, big2cs, nbig2))
      return
c
 40   if (x .gt. xmax) call xermsg ('slatec', 'bi',
     +   'x so big that bi overflows', 1, 2)
c
      bi = bie(x) * exp(2.0*x*sqrt(x)/3.0)
      return
c
      end
