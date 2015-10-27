*deck e1
      function e1 (x)
c***begin prologue  e1
c***purpose  compute the exponential integral e1(x).
c***library   slatec (fnlib)
c***category  c5
c***type      single precision (e1-s, de1-d)
c***keywords  e1 function, exponential integral, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c e1 calculates the single precision exponential integral, e1(x), for
c positive single precision argument x and the cauchy principal value
c for negative x.  if principal values are used everywhere, then, for
c all x,
c
c    e1(x) = -ei(-x)
c or
c    ei(x) = -e1(-x).
c
c
c series for ae11       on the interval -1.00000d-01 to  0.
c                                        with weighted error   1.76e-17
c                                         log weighted error  16.75
c                               significant figures required  15.70
c                                    decimal places required  17.55
c
c
c series for ae12       on the interval -2.50000d-01 to -1.00000d-01
c                                        with weighted error   5.83e-17
c                                         log weighted error  16.23
c                               significant figures required  15.76
c                                    decimal places required  16.93
c
c
c series for e11        on the interval -4.00000d+00 to -1.00000d+00
c                                        with weighted error   1.08e-18
c                                         log weighted error  17.97
c                               significant figures required  19.02
c                                    decimal places required  18.61
c
c
c series for e12        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   3.15e-18
c                                         log weighted error  17.50
c                        approx significant figures required  15.8
c                                    decimal places required  18.10
c
c
c series for ae13       on the interval  2.50000d-01 to  1.00000d+00
c                                        with weighted error   2.34e-17
c                                         log weighted error  16.63
c                               significant figures required  16.14
c                                    decimal places required  17.33
c
c
c series for ae14       on the interval  0.          to  2.50000d-01
c                                        with weighted error   5.41e-17
c                                         log weighted error  16.27
c                               significant figures required  15.38
c                                    decimal places required  16.97
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891115  modified prologue description.  (wrb)
c   891115  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  e1
      dimension ae11cs(39), ae12cs(25), e11cs(19), e12cs(16),
     1  ae13cs(25), ae14cs(26)
      logical first
      save ae11cs, ae12cs, e11cs, e12cs, ae13cs, ae14cs,
     1 ntae11, ntae12, nte11, nte12, ntae13, ntae14, xmax, first
      data ae11cs( 1) /    .1215032397 1606579e0 /
      data ae11cs( 2) /   -.0650887785 13550150e0 /
      data ae11cs( 3) /    .0048976513 57459670e0 /
      data ae11cs( 4) /   -.0006492378 43027216e0 /
      data ae11cs( 5) /    .0000938404 34587471e0 /
      data ae11cs( 6) /    .0000004202 36380882e0 /
      data ae11cs( 7) /   -.0000081133 74735904e0 /
      data ae11cs( 8) /    .0000028042 47688663e0 /
      data ae11cs( 9) /    .0000000564 87164441e0 /
      data ae11cs(10) /   -.0000003448 09174450e0 /
      data ae11cs(11) /    .0000000582 09273578e0 /
      data ae11cs(12) /    .0000000387 11426349e0 /
      data ae11cs(13) /   -.0000000124 53235014e0 /
      data ae11cs(14) /   -.0000000051 18504888e0 /
      data ae11cs(15) /    .0000000021 48771527e0 /
      data ae11cs(16) /    .0000000008 68459898e0 /
      data ae11cs(17) /   -.0000000003 43650105e0 /
      data ae11cs(18) /   -.0000000001 79796603e0 /
      data ae11cs(19) /    .0000000000 47442060e0 /
      data ae11cs(20) /    .0000000000 40423282e0 /
      data ae11cs(21) /   -.0000000000 03543928e0 /
      data ae11cs(22) /   -.0000000000 08853444e0 /
      data ae11cs(23) /   -.0000000000 00960151e0 /
      data ae11cs(24) /    .0000000000 01692921e0 /
      data ae11cs(25) /    .0000000000 00607990e0 /
      data ae11cs(26) /   -.0000000000 00224338e0 /
      data ae11cs(27) /   -.0000000000 00200327e0 /
      data ae11cs(28) /   -.0000000000 00006246e0 /
      data ae11cs(29) /    .0000000000 00045571e0 /
      data ae11cs(30) /    .0000000000 00016383e0 /
      data ae11cs(31) /   -.0000000000 00005561e0 /
      data ae11cs(32) /   -.0000000000 00006074e0 /
      data ae11cs(33) /   -.0000000000 00000862e0 /
      data ae11cs(34) /    .0000000000 00001223e0 /
      data ae11cs(35) /    .0000000000 00000716e0 /
      data ae11cs(36) /   -.0000000000 00000024e0 /
      data ae11cs(37) /   -.0000000000 00000201e0 /
      data ae11cs(38) /   -.0000000000 00000082e0 /
      data ae11cs(39) /    .0000000000 00000017e0 /
      data ae12cs( 1) /    .5824174951 3472674e0 /
      data ae12cs( 2) /   -.1583488509 0578275e0 /
      data ae12cs( 3) /   -.0067642755 90323141e0 /
      data ae12cs( 4) /    .0051258439 50185725e0 /
      data ae12cs( 5) /    .0004352324 92169391e0 /
      data ae12cs( 6) /   -.0001436133 66305483e0 /
      data ae12cs( 7) /   -.0000418013 20556301e0 /
      data ae12cs( 8) /   -.0000027133 95758640e0 /
      data ae12cs( 9) /    .0000011513 81913647e0 /
      data ae12cs(10) /    .0000004206 50022012e0 /
      data ae12cs(11) /    .0000000665 81901391e0 /
      data ae12cs(12) /    .0000000006 62143777e0 /
      data ae12cs(13) /   -.0000000028 44104870e0 /
      data ae12cs(14) /   -.0000000009 40724197e0 /
      data ae12cs(15) /   -.0000000001 77476602e0 /
      data ae12cs(16) /   -.0000000000 15830222e0 /
      data ae12cs(17) /    .0000000000 02905732e0 /
      data ae12cs(18) /    .0000000000 01769356e0 /
      data ae12cs(19) /    .0000000000 00492735e0 /
      data ae12cs(20) /    .0000000000 00093709e0 /
      data ae12cs(21) /    .0000000000 00010707e0 /
      data ae12cs(22) /   -.0000000000 00000537e0 /
      data ae12cs(23) /   -.0000000000 00000716e0 /
      data ae12cs(24) /   -.0000000000 00000244e0 /
      data ae12cs(25) /   -.0000000000 00000058e0 /
      data e11cs( 1) / -16.1134616555 71494026e0 /
      data e11cs( 2) /   7.7940727787 426802769e0 /
      data e11cs( 3) /  -1.9554058188 631419507e0 /
      data e11cs( 4) /    .3733729386 6277945612e0 /
      data e11cs( 5) /   -.0569250319 1092901938e0 /
      data e11cs( 6) /    .0072110777 6966009185e0 /
      data e11cs( 7) /   -.0007810490 1449841593e0 /
      data e11cs( 8) /    .0000738809 3356262168e0 /
      data e11cs( 9) /   -.0000062028 6187580820e0 /
      data e11cs(10) /    .0000004681 6002303176e0 /
      data e11cs(11) /   -.0000000320 9288853329e0 /
      data e11cs(12) /    .0000000020 1519974874e0 /
      data e11cs(13) /   -.0000000001 1673686816e0 /
      data e11cs(14) /    .0000000000 0627627066e0 /
      data e11cs(15) /   -.0000000000 0031481541e0 /
      data e11cs(16) /    .0000000000 0001479904e0 /
      data e11cs(17) /   -.0000000000 0000065457e0 /
      data e11cs(18) /    .0000000000 0000002733e0 /
      data e11cs(19) /   -.0000000000 0000000108e0 /
      data e12cs( 1) /  -0.0373902147 92202795e0 /
      data e12cs( 2) /   0.0427239860 62209577e0 /
      data e12cs( 3) /   -.1303182079 849700544e0 /
      data e12cs( 4) /    .0144191240 2469889073e0 /
      data e12cs( 5) /   -.0013461707 8051068022e0 /
      data e12cs( 6) /    .0001073102 9253063780e0 /
      data e12cs( 7) /   -.0000074299 9951611943e0 /
      data e12cs( 8) /    .0000004537 7325690753e0 /
      data e12cs( 9) /   -.0000000247 6417211390e0 /
      data e12cs(10) /    .0000000012 2076581374e0 /
      data e12cs(11) /   -.0000000000 5485141480e0 /
      data e12cs(12) /    .0000000000 0226362142e0 /
      data e12cs(13) /   -.0000000000 0008635897e0 /
      data e12cs(14) /    .0000000000 0000306291e0 /
      data e12cs(15) /   -.0000000000 0000010148e0 /
      data e12cs(16) /    .0000000000 0000000315e0 /
      data ae13cs( 1) /   -.6057732466 4060346e0 /
      data ae13cs( 2) /   -.1125352434 8366090e0 /
      data ae13cs( 3) /    .0134322662 47902779e0 /
      data ae13cs( 4) /   -.0019268451 87381145e0 /
      data ae13cs( 5) /    .0003091183 37720603e0 /
      data ae13cs( 6) /   -.0000535641 32129618e0 /
      data ae13cs( 7) /    .0000098278 12880247e0 /
      data ae13cs( 8) /   -.0000018853 68984916e0 /
      data ae13cs( 9) /    .0000003749 43193568e0 /
      data ae13cs(10) /   -.0000000768 23455870e0 /
      data ae13cs(11) /    .0000000161 43270567e0 /
      data ae13cs(12) /   -.0000000034 66802211e0 /
      data ae13cs(13) /    .0000000007 58754209e0 /
      data ae13cs(14) /   -.0000000001 68864333e0 /
      data ae13cs(15) /    .0000000000 38145706e0 /
      data ae13cs(16) /   -.0000000000 08733026e0 /
      data ae13cs(17) /    .0000000000 02023672e0 /
      data ae13cs(18) /   -.0000000000 00474132e0 /
      data ae13cs(19) /    .0000000000 00112211e0 /
      data ae13cs(20) /   -.0000000000 00026804e0 /
      data ae13cs(21) /    .0000000000 00006457e0 /
      data ae13cs(22) /   -.0000000000 00001568e0 /
      data ae13cs(23) /    .0000000000 00000383e0 /
      data ae13cs(24) /   -.0000000000 00000094e0 /
      data ae13cs(25) /    .0000000000 00000023e0 /
      data ae14cs( 1) /   -.1892918000 753017e0 /
      data ae14cs( 2) /   -.0864811785 5259871e0 /
      data ae14cs( 3) /    .0072241015 4374659e0 /
      data ae14cs( 4) /   -.0008097559 4575573e0 /
      data ae14cs( 5) /    .0001099913 4432661e0 /
      data ae14cs( 6) /   -.0000171733 2998937e0 /
      data ae14cs( 7) /    .0000029856 2751447e0 /
      data ae14cs( 8) /   -.0000005659 6491457e0 /
      data ae14cs( 9) /    .0000001152 6808397e0 /
      data ae14cs(10) /   -.0000000249 5030440e0 /
      data ae14cs(11) /    .0000000056 9232420e0 /
      data ae14cs(12) /   -.0000000013 5995766e0 /
      data ae14cs(13) /    .0000000003 3846628e0 /
      data ae14cs(14) /   -.0000000000 8737853e0 /
      data ae14cs(15) /    .0000000000 2331588e0 /
      data ae14cs(16) /   -.0000000000 0641148e0 /
      data ae14cs(17) /    .0000000000 0181224e0 /
      data ae14cs(18) /   -.0000000000 0052538e0 /
      data ae14cs(19) /    .0000000000 0015592e0 /
      data ae14cs(20) /   -.0000000000 0004729e0 /
      data ae14cs(21) /    .0000000000 0001463e0 /
      data ae14cs(22) /   -.0000000000 0000461e0 /
      data ae14cs(23) /    .0000000000 0000148e0 /
      data ae14cs(24) /   -.0000000000 0000048e0 /
      data ae14cs(25) /    .0000000000 0000016e0 /
      data ae14cs(26) /   -.0000000000 0000005e0 /
      data first /.true./
c***first executable statement  e1
      if (first) then
         eta = 0.1*r1mach(3)
         ntae11 = inits (ae11cs, 39, eta)
         ntae12 = inits (ae12cs, 25, eta)
         nte11 = inits (e11cs, 19, eta)
         nte12 = inits (e12cs, 16, eta)
         ntae13 = inits (ae13cs, 25, eta)
         ntae14 = inits (ae14cs, 26, eta)
c
         xmaxt = -log (r1mach(1))
         xmax = xmaxt - log(xmaxt)
      endif
      first = .false.
c
      if (x.gt.(-10.)) go to 20
c
c e1(x) = -ei(-x) for x .le. -10.
c
      e1 = exp(-x)/x * (1.+csevl (20./x+1., ae11cs, ntae11))
      return
c
 20   if (x.gt.(-4.0)) go to 30
c
c e1(x) = -ei(-x) for -10. .lt. x .le. -4.
c
      e1 = exp(-x)/x * (1.+csevl ((40./x+7.)/3., ae12cs, ntae12))
      return
c
 30   if (x.gt.(-1.0)) go to 40
c
c e1(x) = -ei(-x) for -4. .lt. x .le. -1.
c
      e1 = -log(abs(x)) + csevl ((2.*x+5.)/3., e11cs, nte11)
      return
c
 40   if (x.gt.1.) go to 50
      if (x .eq. 0.) call xermsg ('slatec', 'e1', 'x is 0', 2, 2)
c
c e1(x) = -ei(-x) for -1. .lt. x .le. 1.,  x .ne. 0.
c
      e1 = (-log(abs(x)) - 0.6875 + x) + csevl (x, e12cs, nte12)
      return
c
 50   if (x.gt.4.) go to 60
c
c e1(x) = -ei(-x) for 1. .lt. x .le. 4.
c
      e1 = exp(-x)/x * (1.+csevl ((8./x-5.)/3., ae13cs, ntae13))
      return
c
 60   if (x.gt.xmax) go to 70
c
c e1(x) = -ei(-x) for 4. .lt. x .le. xmax
c
      e1 = exp(-x)/x * (1. + csevl (8./x-1., ae14cs, ntae14))
      return
c
c e1(x) = -ei(-x) for x .gt. xmax
c
 70   call xermsg ('slatec', 'e1', 'x so big e1 underflows', 1, 1)
      e1 = 0.
      return
c
      end
