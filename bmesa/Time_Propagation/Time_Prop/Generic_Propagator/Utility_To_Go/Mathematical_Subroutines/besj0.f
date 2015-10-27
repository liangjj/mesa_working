*deck besj0
      function besj0 (x)
c***begin prologue  besj0
c***purpose  compute the bessel function of the first kind of order
c            zero.
c***library   slatec (fnlib)
c***category  c10a1
c***type      double precision 
c***keywords  bessel function, first kind, fnlib, order zero,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c besj0(x) calculates the double precision bessel function of
c the first kind of order zero for double precision argument x.
c
c series for bj0        on the interval  0.          to  1.60000e+01
c                                        with weighted error   4.39e-32
c                                         log weighted error  31.36
c                               significant figures required  31.21
c                                    decimal places required  32.00
c
c***references  (none)
c***routines called  r1mach, r9b0mp, csevl, initds
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  besj0
      real*8 x, bj0cs(19), ampl, theta, xsml, y, r1mach,
     1  csevl
      logical first
      save bj0cs, ntj0, xsml, first
      data bj0cs(  1) / +.1002541619 6893913701 0731272640 74 d+0     /
      data bj0cs(  2) / -.6652230077 6440513177 6787578311 24 d+0     /
      data bj0cs(  3) / +.2489837034 9828131370 4604687266 80 d+0     /
      data bj0cs(  4) / -.3325272317 0035769653 8843415038 54 d-1     /
      data bj0cs(  5) / +.2311417930 4694015462 9049241177 29 d-2     /
      data bj0cs(  6) / -.9911277419 9508092339 0485193365 49 d-4     /
      data bj0cs(  7) / +.2891670864 3998808884 7339037470 78 d-5     /
      data bj0cs(  8) / -.6121085866 3032635057 8184074815 16 d-7     /
      data bj0cs(  9) / +.9838650793 8567841324 7687486364 15 d-9     /
      data bj0cs( 10) / -.1242355159 7301765145 5158970068 36 d-10    /
      data bj0cs( 11) / +.1265433630 2559045797 9158272103 63 d-12    /
      data bj0cs( 12) / -.1061945649 5287244546 9148175129 59 d-14    /
      data bj0cs( 13) / +.7470621075 8024567437 0989155840 00 d-17    /
      data bj0cs( 14) / -.4469703227 4412780547 6270079999 99 d-19    /
      data bj0cs( 15) / +.2302428158 4337436200 5230933333 33 d-21    /
      data bj0cs( 16) / -.1031914479 4166698148 5226666666 66 d-23    /
      data bj0cs( 17) / +.4060817827 4873322700 8000000000 00 d-26    /
      data bj0cs( 18) / -.1414383600 5240913919 9999999999 99 d-28    /
      data bj0cs( 19) / +.4391090549 6698880000 0000000000 00 d-31    /
      data first /.true./
c***first executable statement  besj0
      if (first) then
         ntj0 = initds (bj0cs, 19, 0.1*real(r1mach(3)))
         xsml = sqrt(8.0d0*r1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.4.0d0) go to 20
c
      besj0 = 1.0d0
      if (y.gt.xsml) besj0 = csevl (.125d0*y*y-1.d0, bj0cs, ntj0)
      return
c
 20   call d9b0mp (y, ampl, theta)
      besj0 = ampl * cos(theta)
c
      return
      end
