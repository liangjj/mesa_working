*deck erfc      
      function erfc(x)
c***begin prologue  erfc
c***date written   770701   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c8a,l5a1e
c***keywords  complementary error function,erfc,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the complementary error function (erfc).
c***description
c
c erfc(x) calculates the single precision complementary error
c function for single precision argument x.
c
c series for erf        on the interval  0.          to  1.00000d+00
c                                        with weighted error   7.10d-18
c                                         log weighted error  17.15
c          significant figures required  16.31
c                                    decimal places required  17.71
c
c series for erfc       on the interval  0.          to  2.50000d-01
c                                        with weighted error   4.81d-17
c                                         log weighted error  16.32
c                        approx significant figures required  15.0
c
c
c series for erc2       on the interval  2.50000d-01 to  1.00000d+00
c                                        with weighted error   5.22e-17
c                                         log weighted error  16.28
c                        approx significant figures required  15.0
c                                    decimal places required  16.96
c***references  (none)
c***routines called  csevl,inits,r1mach,xerror
c***end prologue  erfc
      implicit real*8 (a-h,o-z)
      dimension erfcs(13), erfccs(24), erc2cs(23)
      data erf cs( 1) /   -.0490461212 34691808d0 /
      data erf cs( 2) /   -.1422612051 0371364d0 /
      data erf cs( 3) /    .0100355821 87599796d0 /
      data erf cs( 4) /   -.0005768764 69976748d0 /
      data erf cs( 5) /    .0000274199 31252196d0 /
      data erf cs( 6) /   -.0000011043 17550734d0 /
      data erf cs( 7) /    .0000000384 88755420d0 /
      data erf cs( 8) /   -.0000000011 80858253d0 /
      data erf cs( 9) /    .0000000000 32334215d0 /
      data erf cs(10) /   -.0000000000 00799101d0 /
      data erf cs(11) /    .0000000000 00017990d0 /
      data erf cs(12) /   -.0000000000 00000371d0 /
      data erf cs(13) /    .0000000000 00000007d0 /
      data erc2cs( 1) /   -.0696013466 02309501d0 /
      data erc2cs( 2) /   -.0411013393 62620893d0 /
      data erc2cs( 3) /    .0039144958 66689626d0 /
      data erc2cs( 4) /   -.0004906395 65054897d0 /
      data erc2cs( 5) /    .0000715747 90013770d0 /
      data erc2cs( 6) /   -.0000115307 16341312d0 /
      data erc2cs( 7) /    .0000019946 70590201d0 /
      data erc2cs( 8) /   -.0000003642 66647159d0 /
      data erc2cs( 9) /    .0000000694 43726100d0 /
      data erc2cs(10) /   -.0000000137 12209021d0 /
      data erc2cs(11) /    .0000000027 88389661d0 /
      data erc2cs(12) /   -.0000000005 81416472d0 /
      data erc2cs(13) /    .0000000001 23892049d0 /
      data erc2cs(14) /   -.0000000000 26906391d0 /
      data erc2cs(15) /    .0000000000 05942614d0 /
      data erc2cs(16) /   -.0000000000 01332386d0 /
      data erc2cs(17) /    .0000000000 00302804d0 /
      data erc2cs(18) /   -.0000000000 00069666d0 /
      data erc2cs(19) /    .0000000000 00016208d0 /
      data erc2cs(20) /   -.0000000000 00003809d0 /
      data erc2cs(21) /    .0000000000 00000904d0 /
      data erc2cs(22) /   -.0000000000 00000216d0 /
      data erc2cs(23) /    .0000000000 00000052d0 /
      data erfccs( 1) /   0.0715179310 202925d0 /
      data erfccs( 2) /   -.0265324343 37606719d0 /
      data erfccs( 3) /    .0017111539 77920853d0 /
      data erfccs( 4) /   -.0001637516 63458512d0 /
      data erfccs( 5) /    .0000198712 93500549d0 /
      data erfccs( 6) /   -.0000028437 12412769d0 /
      data erfccs( 7) /    .0000004606 16130901d0 /
      data erfccs( 8) /   -.0000000822 77530261d0 /
      data erfccs( 9) /    .0000000159 21418724d0 /
      data erfccs(10) /   -.0000000032 95071356d0 /
      data erfccs(11) /    .0000000007 22343973d0 /
      data erfccs(12) /   -.0000000001 66485584d0 /
      data erfccs(13) /    .0000000000 40103931d0 /
      data erfccs(14) /   -.0000000000 10048164d0 /
      data erfccs(15) /    .0000000000 02608272d0 /
      data erfccs(16) /   -.0000000000 00699105d0 /
      data erfccs(17) /    .0000000000 00192946d0 /
      data erfccs(18) /   -.0000000000 00054704d0 /
      data erfccs(19) /    .0000000000 00015901d0 /
      data erfccs(20) /   -.0000000000 00004729d0 /
      data erfccs(21) /    .0000000000 00001432d0 /
      data erfccs(22) /   -.0000000000 00000439d0 /
      data erfccs(23) /    .0000000000 00000138d0 /
      data erfccs(24) /   -.0000000000 00000048d0 /
      data sqrtpi /1.772453850 9055160d0/
      data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0.d0/
c***first executable statement  erfc
      if (nterf.ne.0) go to 10
      eta = 0.1d0*r1mach(3)
      nterf = inits (erfcs, 13, eta)
      nterfc = inits (erfccs, 24, eta)
      nterc2 = inits (erc2cs, 23, eta)
c
      xsml = -sqrt (-log(sqrtpi*r1mach(3)))
      xmax = sqrt (-log(sqrtpi*r1mach(1)))
      xmax = xmax - 0.5d0*log(xmax)/xmax - 0.01d0
      sqeps = sqrt (2.0d0*r1mach(3))
c
 10   if (x.gt.xsml) go to 20
c
c erfc(x) = 1.0 - erf(x) for x .lt. xsml
c
      erfc = 2.d0
      return
c
 20   if (x.gt.xmax) go to 40
      y = abs(x)
      if (y.gt.1.0d0) go to 30
c
c erfc(x) = 1.0 - erf(x) for -1. .le. x .le. 1.
c
      if (y.lt.sqeps) erfc = 1.0d0 - 2.0d0*x/sqrtpi
      if (y.ge.sqeps) erfc = 1.0d0 -
     1  x*(1.0d0 + csevl (2.d0*x*x-1.d0, erfcs, nterf) )
      return
c
c erfc(x) = 1.0 - erf(x) for 1. .lt. abs(x) .le. xmax
c
 30   y = y*y
      if (y.le.4.d0) erfc = exp(-y)/abs(x) * (0.5d0 + 
     1           csevl ((8.d0/y-5.d0)/3.d0,erc2cs, nterc2) )
      if (y.gt.4.d0) erfc = exp(-y)/abs(x) * (0.5d0 + 
     1           csevl (8.d0/y-1.d0,erfccs, nterfc) )
      if (x.lt.0.d0) erfc = 2.0d0 - erfc
      return
c
 40   call lnkerr ( 'erfc    x so big erfc underflows')
      erfc = 0.d0
      return
c
      end
