*deck dbesj1
      double precision function dbesj1 (x)
c***begin prologue  dbesj1
c***purpose  compute the bessel function of the first kind of order one.
c***library   slatec (fnlib)
c***category  c10a1
c***type      double precision (besj1-s, dbesj1-d)
c***keywords  bessel function, first kind, fnlib, order one,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbesj1(x) calculates the double precision bessel function of the
c first kind of order one for double precision argument x.
c
c series for bj1        on the interval  0.          to  1.60000e+01
c                                        with weighted error   1.16e-33
c                                         log weighted error  32.93
c                               significant figures required  32.36
c                                    decimal places required  33.57
c
c***references  (none)
c***routines called  d1mach, d9b1mp, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   780601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   910401  corrected error in code which caused values to have the
c           wrong sign for arguments less than 4.0.  (wrb)
c***end prologue  dbesj1
      double precision x, bj1cs(19), ampl, theta, xsml, xmin, y,
     1  d1mach, dcsevl
      logical first
      save bj1cs, ntj1, xsml, xmin, first
      data bj1cs(  1) / -.1172614151 3332786560 6240574524 003 d+0    /
      data bj1cs(  2) / -.2536152183 0790639562 3030884554 698 d+0    /
      data bj1cs(  3) / +.5012708098 4469568505 3656363203 743 d-1    /
      data bj1cs(  4) / -.4631514809 6250819184 2619728789 772 d-2    /
      data bj1cs(  5) / +.2479962294 1591402453 9124064592 364 d-3    /
      data bj1cs(  6) / -.8678948686 2788258452 1246435176 416 d-5    /
      data bj1cs(  7) / +.2142939171 4379369150 2766250991 292 d-6    /
      data bj1cs(  8) / -.3936093079 1831797922 9322764073 061 d-8    /
      data bj1cs(  9) / +.5591182317 9468800401 8248059864 032 d-10   /
      data bj1cs( 10) / -.6327616404 6613930247 7695274014 880 d-12   /
      data bj1cs( 11) / +.5840991610 8572470032 6945563268 266 d-14   /
      data bj1cs( 12) / -.4482533818 7012581903 9135059199 999 d-16   /
      data bj1cs( 13) / +.2905384492 6250246630 6018688000 000 d-18   /
      data bj1cs( 14) / -.1611732197 8414416541 2118186666 666 d-20   /
      data bj1cs( 15) / +.7739478819 3927463729 8346666666 666 d-23   /
      data bj1cs( 16) / -.3248693782 1119984114 3466666666 666 d-25   /
      data bj1cs( 17) / +.1202237677 2274102272 0000000000 000 d-27   /
      data bj1cs( 18) / -.3952012212 6513493333 3333333333 333 d-30   /
      data bj1cs( 19) / +.1161678082 2664533333 3333333333 333 d-32   /
      data first /.true./
c***first executable statement  dbesj1
      if (first) then
         ntj1 = initds (bj1cs, 19, 0.1*real(d1mach(3)))
c
         xsml = sqrt(8.0d0*d1mach(3))
         xmin = 2.0d0*d1mach(1)
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.4.0d0) go to 20
c
      dbesj1 = 0.0d0
      if (y.eq.0.0d0) return
      if (y .le. xmin) call xermsg ('slatec', 'dbesj1',
     +   'abs(x) so small j1 underflows', 1, 1)
      if (y.gt.xmin) dbesj1 = 0.5d0*x
      if (y.gt.xsml) dbesj1 = x*(.25d0 + dcsevl (.125d0*y*y-1.d0,
     1  bj1cs, ntj1) )
      return
c
 20   call d9b1mp (y, ampl, theta)
      dbesj1 = sign (ampl, x) * cos(theta)
c
      return
      end
