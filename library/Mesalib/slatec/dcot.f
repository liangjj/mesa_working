*deck dcot
      double precision function dcot (x)
c***begin prologue  dcot
c***purpose  compute the cotangent.
c***library   slatec (fnlib)
c***category  c4a
c***type      double precision (cot-s, dcot-d, ccot-c)
c***keywords  cotangent, elementary functions, fnlib, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c dcot(x) calculates the double precision trigonometric cotangent
c for double precision argument x.  x is in units of radians.
c
c series for cot        on the interval  0.          to  6.25000e-02
c                                        with weighted error   5.52e-34
c                                         log weighted error  33.26
c                               significant figures required  32.34
c                                    decimal places required  33.85
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  dcot
      double precision x, cotcs(15), ainty, ainty2, pi2rec, sqeps,
     1  xmax, xmin, xsml, y, yrem, prodbg, dcsevl, d1mach
      logical first
      save cotcs, pi2rec, nterms, xmax, xsml, xmin, sqeps, first
      data cotcs(  1) / +.2402591609 8295630250 9553617744 970 d+0    /
      data cotcs(  2) / -.1653303160 1500227845 4746025255 758 d-1    /
      data cotcs(  3) / -.4299839193 1724018935 6476228239 895 d-4    /
      data cotcs(  4) / -.1592832233 2754104602 3490851122 445 d-6    /
      data cotcs(  5) / -.6191093135 1293487258 8620579343 187 d-9    /
      data cotcs(  6) / -.2430197415 0726460433 1702590579 575 d-11   /
      data cotcs(  7) / -.9560936758 8000809842 7062083100 000 d-14   /
      data cotcs(  8) / -.3763537981 9458058041 6291539706 666 d-16   /
      data cotcs(  9) / -.1481665746 4674657885 2176794666 666 d-18   /
      data cotcs( 10) / -.5833356589 0366657947 7984000000 000 d-21   /
      data cotcs( 11) / -.2296626469 6464577392 8533333333 333 d-23   /
      data cotcs( 12) / -.9041970573 0748332671 9999999999 999 d-26   /
      data cotcs( 13) / -.3559885519 2060006400 0000000000 000 d-28   /
      data cotcs( 14) / -.1401551398 2429866666 6666666666 666 d-30   /
      data cotcs( 15) / -.5518004368 7253333333 3333333333 333 d-33   /
      data pi2rec / .01161977236 7581343075 5350534900 57 d0 /
      data first /.true./
c***first executable statement  dcot
      if (first) then
         nterms = initds (cotcs, 15, 0.1*real(d1mach(3)) )
         xmax = 1.0d0/d1mach(4)
         xsml = sqrt(3.0d0*d1mach(3))
         xmin = exp (max(log(d1mach(1)), -log(d1mach(2))) + 0.01d0)
         sqeps = sqrt(d1mach(4))
      endif
      first = .false.
c
      y = abs(x)
      if (y .lt. xmin) call xermsg ('slatec', 'dcot',
     +   'abs(x) is zero or so small dcot overflows', 2, 2)
      if (y .gt. xmax) call xermsg ('slatec', 'dcot',
     +   'no precision because abs(x) is too big', 3, 2)
c
c carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
c = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
c = aint(.625*y) + aint(z) + rem(z)
c
      ainty = aint (y)
      yrem = y - ainty
      prodbg = 0.625d0*ainty
      ainty = aint (prodbg)
      y = (prodbg-ainty) + 0.625d0*yrem + pi2rec*y
      ainty2 = aint (y)
      ainty = ainty + ainty2
      y = y - ainty2
c
      ifn = mod (ainty, 2.0d0)
      if (ifn.eq.1) y = 1.0d0 - y
c
      if (abs(x) .gt. 0.5d0 .and. y .lt. abs(x)*sqeps) call xermsg
     +   ('slatec', 'dcot',
     +   'answer lt half precision, abs(x) too big or x near n*pi ' //
     +   '(n.ne.0)', 1, 1)
c
      if (y.gt.0.25d0) go to 20
      dcot = 1.0d0/x
      if (y.gt.xsml) dcot = (0.5d0 + dcsevl (32.0d0*y*y-1.d0, cotcs,
     1  nterms)) / y
      go to 40
c
 20   if (y.gt.0.5d0) go to 30
      dcot = (0.5d0 + dcsevl (8.d0*y*y-1.d0, cotcs, nterms))/(0.5d0*y)
      dcot = (dcot*dcot-1.d0)*0.5d0/dcot
      go to 40
c
 30   dcot = (0.5d0 + dcsevl (2.d0*y*y-1.d0, cotcs, nterms))/(.25d0*y)
      dcot = (dcot*dcot-1.d0)*0.5d0/dcot
      dcot = (dcot*dcot-1.d0)*0.5d0/dcot
c
 40   if (x.ne.0.d0) dcot = sign (dcot, x)
      if (ifn.eq.1) dcot = -dcot
c
      return
      end
