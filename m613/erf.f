*deck @(#)erf.f	5.1 11/6/94
      function erf(x)
c***begin prologue  erf.f
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c8a,l5a1e
c***keywords  erf,error function,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the error function, erf.
c***description
c
c erf(x) calculates the single precision error function for
c single precision argument x.
c
c series for erf        on the interval  0.          to  1.00000d+00
c                                        with weighted error   7.10d-18
c                                         log weighted error  17.15
c          significant figures required  16.31
c                                    decimal places required  17.71
c***references  (none)
c***routines called  csevl, erfc, inits, r1mach
c***end prologue  erf
      implicit real*8 (a-h,o-z)
      dimension erfcs(13)
      external erfc
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
      data sqrtpi /1.772453850 9055160d0/
      data nterf, xbig, sqeps / 0, 0.d0, 0.d0/
c***first executable statement  erf
      if (nterf.ne.0) go to 10
      nterf = inits (erfcs, 13, 0.1d0*r1mach(3))
      xbig = sqrt(-log(sqrtpi*r1mach(3)))
      sqeps = sqrt(2.0d0*r1mach(3))
c
 10   y = abs(x)
      if (y.gt.1.d0) go to 20
c
c erf(x) = 1. - erfc(x) for -1. .le. x .le. 1.
c
      if (y.le.sqeps) erf = 2.0d0*x/sqrtpi
      if (y.gt.sqeps) erf = x*(1.0d0 + 
     1                    csevl(2.d0*x**2-1.d0, erfcs, nterf))
      return
c
c erf(x) = 1. - erfc(x) for  abs(x) .gt. 1.
c
 20   if (y.le.xbig) erf = sign (1.d0-erfc(y), x)
      if (y.gt.xbig) erf = sign (1.d0, x)
c
      return
      end
