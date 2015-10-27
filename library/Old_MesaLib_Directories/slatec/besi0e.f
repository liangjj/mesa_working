*deck besi0e
      function besi0e (x)
c***begin prologue  besi0e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the first kind of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besi0e-s, dbsi0e-d)
c***keywords  exponentially scaled, first kind, fnlib,
c             hyperbolic bessel function, modified bessel function,
c             order zero, special functions
c***author  fullerton, w., (lanl)
c***description
c
c besi0e(x) calculates the exponentially scaled modified (hyperbolic)
c bessel function of the first kind of order zero for real argument x;
c i.e., exp(-abs(x))*i0(x).
c
c
c series for bi0        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.46e-18
c                                         log weighted error  17.61
c                               significant figures required  17.90
c                                    decimal places required  18.15
c
c
c series for ai0        on the interval  1.25000d-01 to  3.33333d-01
c                                        with weighted error   7.87e-17
c                                         log weighted error  16.10
c                               significant figures required  14.69
c                                    decimal places required  16.76
c
c
c series for ai02       on the interval  0.          to  1.25000d-01
c                                        with weighted error   3.79e-17
c                                         log weighted error  16.42
c                               significant figures required  14.86
c                                    decimal places required  17.09
c
c***references  (none)
c***routines called  csevl, inits, r1mach
c***revision history  (yymmdd)
c   770701  date written
c   890313  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  besi0e
      dimension bi0cs(12), ai0cs(21), ai02cs(22)
      logical first
      save bi0cs, ai0cs, ai02cs, nti0, ntai0, ntai02, xsml, first
      data bi0cs( 1) /   -.0766054725 2839144951e0 /
      data bi0cs( 2) /   1.9273379539 93808270e0 /
      data bi0cs( 3) /    .2282644586 920301339e0 /
      data bi0cs( 4) /    .0130489146 6707290428e0 /
      data bi0cs( 5) /    .0004344270 9008164874e0 /
      data bi0cs( 6) /    .0000094226 5768600193e0 /
      data bi0cs( 7) /    .0000001434 0062895106e0 /
      data bi0cs( 8) /    .0000000016 1384906966e0 /
      data bi0cs( 9) /    .0000000000 1396650044e0 /
      data bi0cs(10) /    .0000000000 0009579451e0 /
      data bi0cs(11) /    .0000000000 0000053339e0 /
      data bi0cs(12) /    .0000000000 0000000245e0 /
      data ai0cs( 1) /    .0757599449 4023796e0 /
      data ai0cs( 2) /    .0075913808 1082334e0 /
      data ai0cs( 3) /    .0004153131 3389237e0 /
      data ai0cs( 4) /    .0000107007 6463439e0 /
      data ai0cs( 5) /   -.0000079011 7997921e0 /
      data ai0cs( 6) /   -.0000007826 1435014e0 /
      data ai0cs( 7) /    .0000002783 8499429e0 /
      data ai0cs( 8) /    .0000000082 5247260e0 /
      data ai0cs( 9) /   -.0000000120 4463945e0 /
      data ai0cs(10) /    .0000000015 5964859e0 /
      data ai0cs(11) /    .0000000002 2925563e0 /
      data ai0cs(12) /   -.0000000001 1916228e0 /
      data ai0cs(13) /    .0000000000 1757854e0 /
      data ai0cs(14) /    .0000000000 0112822e0 /
      data ai0cs(15) /   -.0000000000 0114684e0 /
      data ai0cs(16) /    .0000000000 0027155e0 /
      data ai0cs(17) /   -.0000000000 0002415e0 /
      data ai0cs(18) /   -.0000000000 0000608e0 /
      data ai0cs(19) /    .0000000000 0000314e0 /
      data ai0cs(20) /   -.0000000000 0000071e0 /
      data ai0cs(21) /    .0000000000 0000007e0 /
      data ai02cs( 1) /    .0544904110 1410882e0 /
      data ai02cs( 2) /    .0033691164 7825569e0 /
      data ai02cs( 3) /    .0000688975 8346918e0 /
      data ai02cs( 4) /    .0000028913 7052082e0 /
      data ai02cs( 5) /    .0000002048 9185893e0 /
      data ai02cs( 6) /    .0000000226 6668991e0 /
      data ai02cs( 7) /    .0000000033 9623203e0 /
      data ai02cs( 8) /    .0000000004 9406022e0 /
      data ai02cs( 9) /    .0000000000 1188914e0 /
      data ai02cs(10) /   -.0000000000 3149915e0 /
      data ai02cs(11) /   -.0000000000 1321580e0 /
      data ai02cs(12) /   -.0000000000 0179419e0 /
      data ai02cs(13) /    .0000000000 0071801e0 /
      data ai02cs(14) /    .0000000000 0038529e0 /
      data ai02cs(15) /    .0000000000 0001539e0 /
      data ai02cs(16) /   -.0000000000 0004151e0 /
      data ai02cs(17) /   -.0000000000 0000954e0 /
      data ai02cs(18) /    .0000000000 0000382e0 /
      data ai02cs(19) /    .0000000000 0000176e0 /
      data ai02cs(20) /   -.0000000000 0000034e0 /
      data ai02cs(21) /   -.0000000000 0000027e0 /
      data ai02cs(22) /    .0000000000 0000003e0 /
      data first /.true./
c***first executable statement  besi0e
      if (first) then
         nti0 = inits (bi0cs, 12, 0.1*r1mach(3))
         ntai0 = inits (ai0cs, 21, 0.1*r1mach(3))
         ntai02 = inits (ai02cs, 22, 0.1*r1mach(3))
         xsml = sqrt (4.5*r1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi0e = 1.0 - x
      if (y.gt.xsml) besi0e = exp(-y) * ( 2.75 +
     1  csevl (y*y/4.5-1.0, bi0cs, nti0) )
      return
c
 20   if (y.le.8.) besi0e = (.375 + csevl ((48./y-11.)/5., ai0cs, ntai0)
     1  ) / sqrt(y)
      if (y.gt.8.) besi0e = (.375 + csevl (16./y-1., ai02cs, ntai02))
     1  / sqrt(y)
c
      return
      end
