*deck bie
      function bie (x)
c***begin prologue  bie
c***purpose  calculate the bairy function for a negative argument and an
c            exponentially scaled bairy function for a non-negative
c            argument.
c***library   slatec (fnlib)
c***category  c10d
c***type      single precision (bie-s, dbie-d)
c***keywords  bairy function, exponentially scaled, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate bi(x) for x .le. 0  and  bi(x)*exp(zeta)  where
c zeta = 2/3 * x**(3/2)  for x .ge. 0.0
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
c series for bip        on the interval  1.25000d-01 to  3.53553d-01
c                                        with weighted error   1.91e-17
c                                         log weighted error  16.72
c                               significant figures required  15.35
c                                    decimal places required  17.41
c
c series for bip2       on the interval  0.          to  1.25000d-01
c                                        with weighted error   1.05e-18
c                                         log weighted error  17.98
c                               significant figures required  16.74
c                                    decimal places required  18.71
c
c***references  (none)
c***routines called  csevl, inits, r1mach, r9aimp
c***revision history  (yymmdd)
c   770701  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  bie
      logical first
      dimension bifcs(9), bigcs(8), bif2cs(10), big2cs(10), bipcs(24),
     1  bip2cs(29)
      save bifcs, bigcs, bif2cs, big2cs, bipcs, bip2cs, atr, btr,
     1 nbif, nbig, nbif2, nbig2, nbip, nbip2, x3sml, x32sml, xbig, first
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
      data bipcs( 1) /   -.0832204747 7943447e0 /
      data bipcs( 2) /    .0114611892 7371174e0 /
      data bipcs( 3) /    .0004289644 0718911e0 /
      data bipcs( 4) /   -.0001490663 9379950e0 /
      data bipcs( 5) /   -.0000130765 9726787e0 /
      data bipcs( 6) /    .0000063275 9839610e0 /
      data bipcs( 7) /   -.0000004222 6696982e0 /
      data bipcs( 8) /   -.0000001914 7186298e0 /
      data bipcs( 9) /    .0000000645 3106284e0 /
      data bipcs(10) /   -.0000000078 4485467e0 /
      data bipcs(11) /   -.0000000009 6077216e0 /
      data bipcs(12) /    .0000000007 0004713e0 /
      data bipcs(13) /   -.0000000001 7731789e0 /
      data bipcs(14) /    .0000000000 2272089e0 /
      data bipcs(15) /    .0000000000 0165404e0 /
      data bipcs(16) /   -.0000000000 0185171e0 /
      data bipcs(17) /    .0000000000 0059576e0 /
      data bipcs(18) /   -.0000000000 0012194e0 /
      data bipcs(19) /    .0000000000 0001334e0 /
      data bipcs(20) /    .0000000000 0000172e0 /
      data bipcs(21) /   -.0000000000 0000145e0 /
      data bipcs(22) /    .0000000000 0000049e0 /
      data bipcs(23) /   -.0000000000 0000011e0 /
      data bipcs(24) /    .0000000000 0000001e0 /
      data bip2cs( 1) /   -.1135967375 85988679e0 /
      data bip2cs( 2) /    .0041381473 947881595e0 /
      data bip2cs( 3) /    .0001353470 622119332e0 /
      data bip2cs( 4) /    .0000104273 166530153e0 /
      data bip2cs( 5) /    .0000013474 954767849e0 /
      data bip2cs( 6) /    .0000001696 537405438e0 /
      data bip2cs( 7) /   -.0000000100 965008656e0 /
      data bip2cs( 8) /   -.0000000167 291194937e0 /
      data bip2cs( 9) /   -.0000000045 815364485e0 /
      data bip2cs(10) /    .0000000003 736681366e0 /
      data bip2cs(11) /    .0000000005 766930320e0 /
      data bip2cs(12) /    .0000000000 621812650e0 /
      data bip2cs(13) /   -.0000000000 632941202e0 /
      data bip2cs(14) /   -.0000000000 149150479e0 /
      data bip2cs(15) /    .0000000000 078896213e0 /
      data bip2cs(16) /    .0000000000 024960513e0 /
      data bip2cs(17) /   -.0000000000 012130075e0 /
      data bip2cs(18) /   -.0000000000 003740493e0 /
      data bip2cs(19) /    .0000000000 002237727e0 /
      data bip2cs(20) /    .0000000000 000474902e0 /
      data bip2cs(21) /   -.0000000000 000452616e0 /
      data bip2cs(22) /   -.0000000000 000030172e0 /
      data bip2cs(23) /    .0000000000 000091058e0 /
      data bip2cs(24) /   -.0000000000 000009814e0 /
      data bip2cs(25) /   -.0000000000 000016429e0 /
      data bip2cs(26) /    .0000000000 000005533e0 /
      data bip2cs(27) /    .0000000000 000002175e0 /
      data bip2cs(28) /   -.0000000000 000001737e0 /
      data bip2cs(29) /   -.0000000000 000000010e0 /
      data atr / 8.750690570 8484345 e0 /
      data btr / -2.093836321 356054 e0 /
      data first /.true./
c***first executable statement  bie
      if (first) then
         eta = 0.1*r1mach(3)
         nbif = inits (bifcs, 9, eta)
         nbig = inits (bigcs, 8, eta)
         nbif2 = inits (bif2cs, 10, eta)
         nbig2 = inits (big2cs, 10, eta)
         nbip  = inits (bipcs , 24, eta)
         nbip2 = inits (bip2cs, 29, eta)
c
         x3sml = eta**0.3333
         x32sml = 1.3104*x3sml**2
         xbig = r1mach(2)**0.6666
      endif
      first = .false.
c
      if (x.ge.(-1.0)) go to 20
      call r9aimp (x, xm, theta)
      bie = xm * sin(theta)
      return
c
 20   if (x.gt.1.0) go to 30
      z = 0.0
      if (abs(x).gt.x3sml) z = x**3
      bie = 0.625 + csevl (z, bifcs, nbif) + x*(0.4375 +
     1  csevl (z, bigcs, nbig))
      if (x.gt.x32sml) bie = bie * exp(-2.0*x*sqrt(x)/3.0)
      return
c
 30   if (x.gt.2.0) go to 40
      z = (2.0*x**3 - 9.0) / 7.0
      bie = exp(-2.0*x*sqrt(x)/3.0) * (1.125 + csevl (z, bif2cs, nbif2)
     1  + x*(0.625 + csevl (z, big2cs, nbig2)) )
      return
c
 40   if (x.gt.4.0) go to 50
      sqrtx = sqrt(x)
      z = atr/(x*sqrtx) + btr
      bie = (0.625 + csevl (z, bipcs, nbip)) / sqrt(sqrtx)
      return
c
 50   sqrtx = sqrt(x)
      z = -1.0
      if (x.lt.xbig) z = 16.0/(x*sqrtx) - 1.0
      bie = (0.625 + csevl (z, bip2cs, nbip2))/sqrt(sqrtx)
      return
c
      end
