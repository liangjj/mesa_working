*deck erfc
      function erfc (x)
c***begin prologue  erfc
c***purpose  compute the complementary error function.
c***library   slatec (fnlib)
c***category  c8a, l5a1e
c***type      single precision (erfc-s, derfc-d)
c***keywords  complementary error function, erfc, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c erfc(x) calculates the single precision complementary error
c function for single precision argument x.
c
c series for erf        on the interval  0.          to  1.00000d+00
c                                        with weighted error   7.10e-18
c                                         log weighted error  17.15
c                               significant figures required  16.31
c                                    decimal places required  17.71
c
c series for erfc       on the interval  0.          to  2.50000d-01
c                                        with weighted error   4.81e-17
c                                         log weighted error  16.32
c                        approx significant figures required  15.0
c
c
c series for erc2       on the interval  2.50000d-01 to  1.00000d+00
c                                        with weighted error   5.22e-17
c                                         log weighted error  16.28
c                        approx significant figures required  15.0
c                                    decimal places required  16.96
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  erfc
      dimension erfcs(13), erfccs(24), erc2cs(23)
      logical first
      save erfcs, erc2cs, erfccs, sqrtpi, nterf, nterfc,
     1 nterc2, xsml, xmax, sqeps, first
      data erfcs( 1) /   -.0490461212 34691808e0 /
      data erfcs( 2) /   -.1422612051 0371364e0 /
      data erfcs( 3) /    .0100355821 87599796e0 /
      data erfcs( 4) /   -.0005768764 69976748e0 /
      data erfcs( 5) /    .0000274199 31252196e0 /
      data erfcs( 6) /   -.0000011043 17550734e0 /
      data erfcs( 7) /    .0000000384 88755420e0 /
      data erfcs( 8) /   -.0000000011 80858253e0 /
      data erfcs( 9) /    .0000000000 32334215e0 /
      data erfcs(10) /   -.0000000000 00799101e0 /
      data erfcs(11) /    .0000000000 00017990e0 /
      data erfcs(12) /   -.0000000000 00000371e0 /
      data erfcs(13) /    .0000000000 00000007e0 /
      data erc2cs( 1) /   -.0696013466 02309501e0 /
      data erc2cs( 2) /   -.0411013393 62620893e0 /
      data erc2cs( 3) /    .0039144958 66689626e0 /
      data erc2cs( 4) /   -.0004906395 65054897e0 /
      data erc2cs( 5) /    .0000715747 90013770e0 /
      data erc2cs( 6) /   -.0000115307 16341312e0 /
      data erc2cs( 7) /    .0000019946 70590201e0 /
      data erc2cs( 8) /   -.0000003642 66647159e0 /
      data erc2cs( 9) /    .0000000694 43726100e0 /
      data erc2cs(10) /   -.0000000137 12209021e0 /
      data erc2cs(11) /    .0000000027 88389661e0 /
      data erc2cs(12) /   -.0000000005 81416472e0 /
      data erc2cs(13) /    .0000000001 23892049e0 /
      data erc2cs(14) /   -.0000000000 26906391e0 /
      data erc2cs(15) /    .0000000000 05942614e0 /
      data erc2cs(16) /   -.0000000000 01332386e0 /
      data erc2cs(17) /    .0000000000 00302804e0 /
      data erc2cs(18) /   -.0000000000 00069666e0 /
      data erc2cs(19) /    .0000000000 00016208e0 /
      data erc2cs(20) /   -.0000000000 00003809e0 /
      data erc2cs(21) /    .0000000000 00000904e0 /
      data erc2cs(22) /   -.0000000000 00000216e0 /
      data erc2cs(23) /    .0000000000 00000052e0 /
      data erfccs( 1) /   0.0715179310 202925e0 /
      data erfccs( 2) /   -.0265324343 37606719e0 /
      data erfccs( 3) /    .0017111539 77920853e0 /
      data erfccs( 4) /   -.0001637516 63458512e0 /
      data erfccs( 5) /    .0000198712 93500549e0 /
      data erfccs( 6) /   -.0000028437 12412769e0 /
      data erfccs( 7) /    .0000004606 16130901e0 /
      data erfccs( 8) /   -.0000000822 77530261e0 /
      data erfccs( 9) /    .0000000159 21418724e0 /
      data erfccs(10) /   -.0000000032 95071356e0 /
      data erfccs(11) /    .0000000007 22343973e0 /
      data erfccs(12) /   -.0000000001 66485584e0 /
      data erfccs(13) /    .0000000000 40103931e0 /
      data erfccs(14) /   -.0000000000 10048164e0 /
      data erfccs(15) /    .0000000000 02608272e0 /
      data erfccs(16) /   -.0000000000 00699105e0 /
      data erfccs(17) /    .0000000000 00192946e0 /
      data erfccs(18) /   -.0000000000 00054704e0 /
      data erfccs(19) /    .0000000000 00015901e0 /
      data erfccs(20) /   -.0000000000 00004729e0 /
      data erfccs(21) /    .0000000000 00001432e0 /
      data erfccs(22) /   -.0000000000 00000439e0 /
      data erfccs(23) /    .0000000000 00000138e0 /
      data erfccs(24) /   -.0000000000 00000048e0 /
      data sqrtpi /1.772453850 9055160e0/
      data first /.true./
c***first executable statement  erfc
      if (first) then
         eta = 0.1*r1mach(3)
         nterf = inits (erfcs, 13, eta)
         nterfc = inits (erfccs, 24, eta)
         nterc2 = inits (erc2cs, 23, eta)
c
         xsml = -sqrt (-log(sqrtpi*r1mach(3)))
         txmax = sqrt (-log(sqrtpi*r1mach(1)))
         xmax = txmax - 0.5*log(txmax)/txmax - 0.01
         sqeps = sqrt (2.0*r1mach(3))
      endif
      first = .false.
c
      if (x.gt.xsml) go to 20
c
c erfc(x) = 1.0 - erf(x) for x .lt. xsml
c
      erfc = 2.
      return
c
 20   if (x.gt.xmax) go to 40
      y = abs(x)
      if (y.gt.1.0) go to 30
c
c erfc(x) = 1.0 - erf(x) for -1. .le. x .le. 1.
c
      if (y.lt.sqeps) erfc = 1.0 - 2.0*x/sqrtpi
      if (y.ge.sqeps) erfc = 1.0 -
     1  x*(1.0 + csevl (2.*x*x-1., erfcs, nterf) )
      return
c
c erfc(x) = 1.0 - erf(x) for 1. .lt. abs(x) .le. xmax
c
 30   y = y*y
      if (y.le.4.) erfc = exp(-y)/abs(x) * (0.5 + csevl ((8./y-5.)/3.,
     1  erc2cs, nterc2) )
      if (y.gt.4.) erfc = exp(-y)/abs(x) * (0.5 + csevl (8./y-1.,
     1  erfccs, nterfc) )
      if (x.lt.0.) erfc = 2.0 - erfc
      return
c
 40   call xermsg ('slatec', 'erfc', 'x so big erfc underflows', 1, 1)
      erfc = 0.
      return
c
      end
