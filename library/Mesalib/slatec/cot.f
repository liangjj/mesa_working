*deck cot
      function cot (x)
c***begin prologue  cot
c***purpose  compute the cotangent.
c***library   slatec (fnlib)
c***category  c4a
c***type      single precision (cot-s, dcot-d, ccot-c)
c***keywords  cotangent, elementary functions, fnlib, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c cot(x) calculates the cotangent of the real argument x.  x is in
c units of radians.
c
c series for cot        on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.76e-17
c                                         log weighted error  16.42
c                               significant figures required  15.51
c                                    decimal places required  16.88
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  cot
      dimension cotcs(8)
      logical first
      save cotcs, pi2rec, nterms, xmax, xsml, xmin, sqeps, first
      data cotcs( 1) /    .2402591609 8295630e0 /
      data cotcs( 2) /   -.0165330316 01500228e0 /
      data cotcs( 3) /   -.0000429983 91931724e0 /
      data cotcs( 4) /   -.0000001592 83223327e0 /
      data cotcs( 5) /   -.0000000006 19109313e0 /
      data cotcs( 6) /   -.0000000000 02430197e0 /
      data cotcs( 7) /   -.0000000000 00009560e0 /
      data cotcs( 8) /   -.0000000000 00000037e0 /
      data pi2rec / .01161977236 75813430 e0 /
      data first /.true./
c***first executable statement  cot
      if (first) then
         nterms = inits (cotcs, 8, 0.1*r1mach(3))
         xmax = 1.0/r1mach(4)
         xsml = sqrt (3.0*r1mach(3))
         xmin = exp ( max(log(r1mach(1)), -log(r1mach(2))) + 0.01)
         sqeps = sqrt (r1mach(4))
      endif
      first = .false.
c
      y = abs(x)
      if (abs(x) .lt. xmin) call xermsg ('slatec', 'cot',
     +   'abs(x) is zero or so small cot overflows', 2, 2)
      if (y .gt. xmax) call xermsg ('slatec', 'cot',
     +   'no precision because abs(x) is too big', 3, 2)
c
c carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
c = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
c = aint(.625*y) + aint(z) + rem(z)
c
      ainty = aint (y)
      yrem = y - ainty
      prodbg = 0.625*ainty
      ainty = aint (prodbg)
      y = (prodbg-ainty) + 0.625*yrem + y*pi2rec
      ainty2 = aint (y)
      ainty = ainty + ainty2
      y = y - ainty2
c
      ifn = mod (ainty, 2.)
      if (ifn.eq.1) y = 1.0 - y
c
      if (abs(x) .gt. 0.5 .and. y .lt. abs(x)*sqeps) call xermsg
     +   ('slatec', 'cot',
     +   'answer lt half precision, abs(x) too big or x near n*pi ' //
     +   '(n.ne.0)' , 1, 1)
c
      if (y.gt.0.25) go to 20
      cot = 1.0/x
      if (y.gt.xsml) cot = (0.5 + csevl (32.0*y*y-1., cotcs, nterms)) /y
      go to 40
c
 20   if (y.gt.0.5) go to 30
      cot = (0.5 + csevl (8.0*y*y-1., cotcs, nterms)) / (0.5*y)
      cot = (cot**2 - 1.0) * 0.5 / cot
      go to 40
c
 30   cot = (0.5 + csevl (2.0*y*y-1., cotcs, nterms)) / (0.25*y)
      cot = (cot**2 - 1.0) * 0.5 / cot
      cot = (cot**2 - 1.0) * 0.5 / cot
c
 40   if (x.ne.0.) cot = sign (cot, x)
      if (ifn.eq.1) cot = -cot
c
      return
      end
