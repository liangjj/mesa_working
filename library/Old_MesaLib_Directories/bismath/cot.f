      function cot(x)
c***begin prologue  cot
c***date written   770601   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c4a
c***keywords  cotangent,elementary function
c***author  fullerton, w., (lanl)
c***purpose  computes the cotangent.
c***description
c
c cot(x) calculates the cotangent of the real argument x.  x is in
c units of radians.
c
c series for cot        on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.76e-17
c                                         log weighted error  16.42
c          significant figures required  15.51
c                                    decimal places required  16.88
c***references  (none)
c***routines called  csevl,inits,r1mach,xerror
c***end prologue  cot
      implicit real*8(a-h,o-z)
      dimension cotcs(8)
      data cot cs( 1) /    .2402591609 8295630d0 /
      data cot cs( 2) /   -.0165330316 01500228d0 /
      data cot cs( 3) /   -.0000429983 91931724d0 /
      data cot cs( 4) /   -.0000001592 83223327d0 /
      data cot cs( 5) /   -.0000000006 19109313d0 /
      data cot cs( 6) /   -.0000000000 02430197d0 /
      data cot cs( 7) /   -.0000000000 00009560d0 /
      data cot cs( 8) /   -.0000000000 00000037d0 /
      data pi2rec / .01161977236 75813430 d0 /
      data nterms, xmax, xsml, xmin, sqeps / 0, 4*0.0d0 /
c***first executable statement  cot
      if (nterms.ne.0) go to 10
      nterms = inits (cotcs, 8, 0.1d0*r1mach(3))
      xmax = 1.0d0/r1mach(4)
      xsml = sqrt (3.0d0*r1mach(3))
      xmin = exp ( max(log(r1mach(1)), -log(r1mach(2))) + 0.01d0)
      sqeps = sqrt (r1mach(4))
c
 10   y = abs(x)
      if (abs(x).lt.xmin) call lnkerr ( 'cot     abs(x) is zero or so sm
     1all cot overflows')
      if (y.gt.xmax) call lnkerr (  'cot     no precision because abs(x)
     1 is big')
c
c carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625d0 + pi2rec)
c = aint(.625d0*y) + rem(.625d0*y) + y*pi2rec  =  aint(.625d0*y) + z
c = aint(.625d0*y) + aint(z) + rem(z)
c
      ainty = aint (y)
      yrem = y - ainty
      prodbg = 0.625d0*ainty
      ainty = aint (prodbg)
      y = (prodbg-ainty) + 0.625d0*yrem + y*pi2rec
      ainty2 = aint (y)
      ainty = ainty + ainty2
      y = y - ainty2
c
      ifn = mod (ainty, 2.d0)
      if (ifn.eq.1) y = 1.0d0 - y
c
      if (abs(x).gt.0.5d0 .and. y.lt.abs(x)*sqeps)
     1    call lnkerr (      'cot answer lt half precision, '//
     2                 'abs(x) too big or x near n*pi (n.ne2.0)')
c
      if (y.gt.0.25d0) go to 20
      cot = 1.0d0/x
      if (y.gt.xsml) cot = (0.5d0 +
     1                      csevl (32.0d0*y*y-1.d0, cotcs, nterms)) /y
      go to 40
c
 20   if (y.gt.0.5d0) go to 30
      cot = (0.5d0 + csevl (8.0d0*y*y-1.d0, cotcs, nterms)) / (0.5d0*y)
      cot = (cot**2 - 1.0d0) * 0.5d0 / cot
      go to 40
c
 30   cot = (0.5d0 + csevl (2.0d0*y*y-1.d0, cotcs, nterms)) / (0.25d0*y)
      cot = (cot**2 - 1.0d0) * 0.5d0 / cot
      cot = (cot**2 - 1.0d0) * 0.5d0 / cot
c
 40   if (x.ne.0.d0) cot = sign (cot, x)
      if (ifn.eq.1) cot = -cot
c
      return
      end
