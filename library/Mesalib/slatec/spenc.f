*deck spenc
      function spenc (x)
c***begin prologue  spenc
c***purpose  compute a form of spence's integral due to k. mitchell.
c***library   slatec (fnlib)
c***category  c5
c***type      single precision (spenc-s, dspenc-d)
c***keywords  fnlib, special functions, spence's integral
c***author  fullerton, w., (lanl)
c***description
c
c evaluate a form of spence's function defined by
c        integral from 0 to x of  -log(1-y)/y  dy.
c for abs(x) .le. 1, the uniformly convergent expansion
c        spenc = sum k=1,infinity  x**k / k**2     is valid.
c
c spence's function can be used to evaluate much more general integral
c forms.  for example,
c        integral from 0 to z of  log(a*x+b)/(c*x+d)  dx  =
c             log(abs(b-a*d/c))*log(abs(a*(c*x+d)/(a*d-b*c)))/c
c             - spenc (a*(c*z+d)/(a*d-b*c)) / c.
c
c ref -- k. mitchell, philosophical magazine, 40, p. 351 (1949).
c        stegun and abromowitz, ams 55, p. 1004.
c
c
c series for spen       on the interval  0.          to  5.00000d-01
c                                        with weighted error   6.82e-17
c                                         log weighted error  16.17
c                               significant figures required  15.22
c                                    decimal places required  16.81
c
c***references  (none)
c***routines called  csevl, inits, r1mach
c***revision history  (yymmdd)
c   780201  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  spenc
      dimension spencs(19)
      logical first
      save spencs, pi26, nspenc, xbig, first
      data spencs( 1) /    .1527365598 892406e0 /
      data spencs( 2) /    .0816965805 8051014e0 /
      data spencs( 3) /    .0058141571 4077873e0 /
      data spencs( 4) /    .0005371619 8145415e0 /
      data spencs( 5) /    .0000572470 4675185e0 /
      data spencs( 6) /    .0000066745 4612164e0 /
      data spencs( 7) /    .0000008276 4673397e0 /
      data spencs( 8) /    .0000001073 3156730e0 /
      data spencs( 9) /    .0000000144 0077294e0 /
      data spencs(10) /    .0000000019 8444202e0 /
      data spencs(11) /    .0000000002 7940058e0 /
      data spencs(12) /    .0000000000 4003991e0 /
      data spencs(13) /    .0000000000 0582346e0 /
      data spencs(14) /    .0000000000 0085767e0 /
      data spencs(15) /    .0000000000 0012768e0 /
      data spencs(16) /    .0000000000 0001918e0 /
      data spencs(17) /    .0000000000 0000290e0 /
      data spencs(18) /    .0000000000 0000044e0 /
      data spencs(19) /    .0000000000 0000006e0 /
      data pi26 / 1.644934066 848226e0 /
      data first /.true./
c***first executable statement  spenc
      if (first) then
         nspenc = inits (spencs, 19, 0.1*r1mach(3))
         xbig = 1.0/r1mach(3)
      endif
      first = .false.
c
      if (x.gt.2.0) go to 60
      if (x.gt.1.0) go to 50
      if (x.gt.0.5) go to 40
      if (x.ge.0.0) go to 30
      if (x.gt.(-1.)) go to 20
c
c here if x .le. -1.0
c
      aln = log(1.0-x)
      spenc = -pi26 - 0.5*aln*(2.0*log(-x)-aln)
      if (x.gt.(-xbig)) spenc = spenc
     1  + (1.0 + csevl (4.0/(1.0-x)-1.0, spencs, nspenc)) / (1.0-x)
      return
c
c -1.0 .lt. x .lt. 0.0
c
 20   spenc = -0.5*log(1.0-x)**2
     1  - x*(1.0 + csevl (4.0*x/(x-1.0)-1.0, spencs, nspenc)) / (x-1.0)
      return
c
c 0.0 .le. x .le. 0.5
c
 30   spenc = x*(1.0 + csevl (4.0*x-1.0, spencs, nspenc))
      return
c
c 0.5 .lt. x .le. 1.0
c
 40   spenc = pi26
      if (x.ne.1.0) spenc = pi26 - log(x)*log(1.0-x)
     1  - (1.0-x)*(1.0 + csevl (4.0*(1.0-x)-1.0, spencs, nspenc))
      return
c
c 1.0 .lt. x .le. 2.0
c
 50   spenc = pi26 - 0.5*log(x)*log((x-1.0)**2/x)
     1  + (x-1.)*(1.0 + csevl (4.0*(x-1.)/x-1.0, spencs, nspenc))/x
      return
c
c x .gt. 2.0
c
 60   spenc = 2.0*pi26 - 0.5*log(x)**2
      if (x.lt.xbig) spenc = spenc
     1  - (1.0 + csevl (4.0/x-1.0, spencs, nspenc))/x
      return
c
      end
