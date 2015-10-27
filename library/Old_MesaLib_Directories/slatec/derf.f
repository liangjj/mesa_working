*deck derf
      double precision function derf (x)
c***begin prologue  derf
c***purpose  compute the error function.
c***library   slatec (fnlib)
c***category  c8a, l5a1e
c***type      double precision (erf-s, derf-d)
c***keywords  erf, error function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c derf(x) calculates the double precision error function for double
c precision argument x.
c
c series for erf        on the interval  0.          to  1.00000e+00
c                                        with weighted error   1.28e-32
c                                         log weighted error  31.89
c                               significant figures required  31.05
c                                    decimal places required  32.55
c
c***references  (none)
c***routines called  d1mach, dcsevl, derfc, initds
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900727  added external statement.  (wrb)
c   920618  removed space from variable name.  (rwc, wrb)
c***end prologue  derf
      double precision x, erfcs(21), sqeps, sqrtpi, xbig, y, d1mach,
     1  dcsevl, derfc
      logical first
      external derfc
      save erfcs, sqrtpi, nterf, xbig, sqeps, first
      data erfcs(  1) / -.4904612123 4691808039 9845440333 76 d-1     /
      data erfcs(  2) / -.1422612051 0371364237 8247418996 31 d+0     /
      data erfcs(  3) / +.1003558218 7599795575 7546767129 33 d-1     /
      data erfcs(  4) / -.5768764699 7674847650 8270255091 67 d-3     /
      data erfcs(  5) / +.2741993125 2196061034 4221607914 71 d-4     /
      data erfcs(  6) / -.1104317550 7344507604 1353812959 05 d-5     /
      data erfcs(  7) / +.3848875542 0345036949 9613114981 74 d-7     /
      data erfcs(  8) / -.1180858253 3875466969 6317518015 81 d-8     /
      data erfcs(  9) / +.3233421582 6050909646 4029309533 54 d-10    /
      data erfcs( 10) / -.7991015947 0045487581 6073747085 95 d-12    /
      data erfcs( 11) / +.1799072511 3961455611 9672454866 34 d-13    /
      data erfcs( 12) / -.3718635487 8186926382 3168282094 93 d-15    /
      data erfcs( 13) / +.7103599003 7142529711 6899083946 66 d-17    /
      data erfcs( 14) / -.1261245511 9155225832 4954248533 33 d-18    /
      data erfcs( 15) / +.2091640694 1769294369 1705002666 66 d-20    /
      data erfcs( 16) / -.3253973102 9314072982 3641600000 00 d-22    /
      data erfcs( 17) / +.4766867209 7976748332 3733333333 33 d-24    /
      data erfcs( 18) / -.6598012078 2851343155 1999999999 99 d-26    /
      data erfcs( 19) / +.8655011469 9637626197 3333333333 33 d-28    /
      data erfcs( 20) / -.1078892517 7498064213 3333333333 33 d-29    /
      data erfcs( 21) / +.1281188399 3017002666 6666666666 66 d-31    /
      data sqrtpi / 1.772453850 9055160272 9816748334 115d0 /
      data first /.true./
c***first executable statement  derf
      if (first) then
         nterf = initds (erfcs, 21, 0.1*real(d1mach(3)))
         xbig = sqrt(-log(sqrtpi*d1mach(3)))
         sqeps = sqrt(2.0d0*d1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.d0) go to 20
c
c erf(x) = 1.0 - erfc(x)  for  -1.0 .le. x .le. 1.0
c
      if (y.le.sqeps) derf = 2.0d0*x*x/sqrtpi
      if (y.gt.sqeps) derf = x*(1.0d0 + dcsevl (2.d0*x*x-1.d0,
     1  erfcs, nterf))
      return
c
c erf(x) = 1.0 - erfc(x) for abs(x) .gt. 1.0
c
 20   if (y.le.xbig) derf = sign (1.0d0-derfc(y), x)
      if (y.gt.xbig) derf = sign (1.0d0, x)
c
      return
      end
