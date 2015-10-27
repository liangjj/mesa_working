*deck dgamma
      double precision function dgamma (x)
c***begin prologue  dgamma
c***purpose  compute the complete gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      double precision (gamma-s, dgamma-d, cgamma-c)
c***keywords  complete gamma function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dgamma(x) calculates the double precision complete gamma function
c for double precision argument x.
c
c series for gam        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.00
c                                    decimal places required  32.05
c
c***references  (none)
c***routines called  d1mach, d9lgmc, dcsevl, dgamlm, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable name.  (rwc, wrb)
c***end prologue  dgamma
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1  xmin, y, d9lgmc, dcsevl, d1mach
      logical first
c
      save gamcs, pi, sq2pil, ngam, xmin, xmax, dxrel, first
      data gamcs(  1) / +.8571195590 9893314219 2006239994 2 d-2      /
      data gamcs(  2) / +.4415381324 8410067571 9131577165 2 d-2      /
      data gamcs(  3) / +.5685043681 5993633786 3266458878 9 d-1      /
      data gamcs(  4) / -.4219835396 4185605010 1250018662 4 d-2      /
      data gamcs(  5) / +.1326808181 2124602205 8400679635 2 d-2      /
      data gamcs(  6) / -.1893024529 7988804325 2394702388 6 d-3      /
      data gamcs(  7) / +.3606925327 4412452565 7808221722 5 d-4      /
      data gamcs(  8) / -.6056761904 4608642184 8554829036 5 d-5      /
      data gamcs(  9) / +.1055829546 3022833447 3182350909 3 d-5      /
      data gamcs( 10) / -.1811967365 5423840482 9185589116 6 d-6      /
      data gamcs( 11) / +.3117724964 7153222777 9025459316 9 d-7      /
      data gamcs( 12) / -.5354219639 0196871408 7408102434 7 d-8      /
      data gamcs( 13) / +.9193275519 8595889468 8778682594 0 d-9      /
      data gamcs( 14) / -.1577941280 2883397617 6742327395 3 d-9      /
      data gamcs( 15) / +.2707980622 9349545432 6654043308 9 d-10     /
      data gamcs( 16) / -.4646818653 8257301440 8166105893 3 d-11     /
      data gamcs( 17) / +.7973350192 0074196564 6076717535 9 d-12     /
      data gamcs( 18) / -.1368078209 8309160257 9949917230 9 d-12     /
      data gamcs( 19) / +.2347319486 5638006572 3347177168 8 d-13     /
      data gamcs( 20) / -.4027432614 9490669327 6657053469 9 d-14     /
      data gamcs( 21) / +.6910051747 3721009121 3833697525 7 d-15     /
      data gamcs( 22) / -.1185584500 2219929070 5238712619 2 d-15     /
      data gamcs( 23) / +.2034148542 4963739552 0102605193 2 d-16     /
      data gamcs( 24) / -.3490054341 7174058492 7401294910 8 d-17     /
      data gamcs( 25) / +.5987993856 4853055671 3505106602 6 d-18     /
      data gamcs( 26) / -.1027378057 8722280744 9006977843 1 d-18     /
      data gamcs( 27) / +.1762702816 0605298249 4275966074 8 d-19     /
      data gamcs( 28) / -.3024320653 7353062609 5877211204 2 d-20     /
      data gamcs( 29) / +.5188914660 2183978397 1783355050 6 d-21     /
      data gamcs( 30) / -.8902770842 4565766924 4925160106 6 d-22     /
      data gamcs( 31) / +.1527474068 4933426022 7459689130 6 d-22     /
      data gamcs( 32) / -.2620731256 1873629002 5732833279 9 d-23     /
      data gamcs( 33) / +.4496464047 8305386703 3104657066 6 d-24     /
      data gamcs( 34) / -.7714712731 3368779117 0390152533 3 d-25     /
      data gamcs( 35) / +.1323635453 1260440364 8657271466 6 d-25     /
      data gamcs( 36) / -.2270999412 9429288167 0231381333 3 d-26     /
      data gamcs( 37) / +.3896418998 0039914493 2081663999 9 d-27     /
      data gamcs( 38) / -.6685198115 1259533277 9212799999 9 d-28     /
      data gamcs( 39) / +.1146998663 1400243843 4761386666 6 d-28     /
      data gamcs( 40) / -.1967938586 3451346772 9510399999 9 d-29     /
      data gamcs( 41) / +.3376448816 5853380903 3489066666 6 d-30     /
      data gamcs( 42) / -.5793070335 7821357846 2549333333 3 d-31     /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data first /.true./
c***first executable statement  dgamma
      if (first) then
         ngam = initds (gamcs, 42, 0.1*real(d1mach(3)) )
c
         call dgamlm (xmin, xmax)
         dxrel = sqrt(d1mach(4))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.10.d0) go to 50
c
c compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
c gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - n
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.0
c
      n = -n
      if (x .eq. 0.d0) call xermsg ('slatec', 'dgamma', 'x is 0', 4, 2)
      if (x .lt. 0.0 .and. x+n-2 .eq. 0.d0) call xermsg ('slatec',
     +   'dgamma', 'x is a negative integer', 4, 2)
      if (x .lt. (-0.5d0) .and. abs((x-aint(x-0.5d0))/x) .lt. dxrel)
     +   call xermsg ('slatec', 'dgamma',
     +   'answer lt half precision because x too near negative integer',
     +   1, 1)
c
      do 20 i=1,n
        dgamma = dgamma/(x+i-1 )
 20   continue
      return
c
c gamma(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dgamma = (y+i) * dgamma
 40   continue
      return
c
c gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).
c
 50   if (x .gt. xmax) call xermsg ('slatec', 'dgamma',
     +   'x so big gamma overflows', 3, 2)
c
      dgamma = 0.d0
      if (x .lt. xmin) call xermsg ('slatec', 'dgamma',
     +   'x so small gamma underflows', 2, 1)
      if (x.lt.xmin) return
c
      dgamma = exp ((y-0.5d0)*log(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
c
      if (abs((x-aint(x-0.5d0))/x) .lt. dxrel) call xermsg ('slatec',
     +   'dgamma',
     +   'answer lt half precision, x too near negative integer', 1, 1)
c
      sinpiy = sin (pi*y)
      if (sinpiy .eq. 0.d0) call xermsg ('slatec', 'dgamma',
     +   'x is a negative integer', 4, 2)
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end
