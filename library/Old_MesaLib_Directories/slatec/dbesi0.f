*deck dbesi0
      double precision function dbesi0 (x)
c***begin prologue  dbesi0
c***purpose  compute the hyperbolic bessel function of the first kind
c            of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besi0-s, dbesi0-d)
c***keywords  first kind, fnlib, hyperbolic bessel function,
c             modified bessel function, order zero, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbesi0(x) calculates the double precision modified (hyperbolic)
c bessel function of the first kind of order zero and double
c precision argument x.
c
c series for bi0        on the interval  0.          to  9.00000e+00
c                                        with weighted error   9.51e-34
c                                         log weighted error  33.02
c                               significant figures required  33.31
c                                    decimal places required  33.65
c
c***references  (none)
c***routines called  d1mach, dbsi0e, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesi0
      double precision x, bi0cs(18), xmax, xsml, y, d1mach,
     1  dcsevl, dbsi0e
      logical first
      save bi0cs, nti0, xsml, xmax, first
      data bi0cs(  1) / -.7660547252 8391449510 8189497624 3285 d-1   /
      data bi0cs(  2) / +.1927337953 9938082699 5240875088 1196 d+1   /
      data bi0cs(  3) / +.2282644586 9203013389 3702929233 0415 d+0   /
      data bi0cs(  4) / +.1304891466 7072904280 7933421069 1888 d-1   /
      data bi0cs(  5) / +.4344270900 8164874513 7868268102 6107 d-3   /
      data bi0cs(  6) / +.9422657686 0019346639 2317174411 8766 d-5   /
      data bi0cs(  7) / +.1434006289 5106910799 6209187817 9957 d-6   /
      data bi0cs(  8) / +.1613849069 6617490699 1541971999 4611 d-8   /
      data bi0cs(  9) / +.1396650044 5356696994 9509270814 2522 d-10  /
      data bi0cs( 10) / +.9579451725 5054453446 2752317189 3333 d-13  /
      data bi0cs( 11) / +.5333981859 8625021310 1510774400 0000 d-15  /
      data bi0cs( 12) / +.2458716088 4374707746 9678591999 9999 d-17  /
      data bi0cs( 13) / +.9535680890 2487700269 4434133333 3333 d-20  /
      data bi0cs( 14) / +.3154382039 7214273367 8933333333 3333 d-22  /
      data bi0cs( 15) / +.9004564101 0946374314 6666666666 6666 d-25  /
      data bi0cs( 16) / +.2240647369 1236700160 0000000000 0000 d-27  /
      data bi0cs( 17) / +.4903034603 2428373333 3333333333 3333 d-30  /
      data bi0cs( 18) / +.9508172606 1226666666 6666666666 6666 d-33  /
      data first /.true./
c***first executable statement  dbesi0
      if (first) then
         nti0 = initds (bi0cs, 18, 0.1*real(d1mach(3)))
         xsml = sqrt(4.5d0*d1mach(3))
         xmax = log (d1mach(2))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0d0) go to 20
c
      dbesi0 = 1.0d0
      if (y.gt.xsml) dbesi0 = 2.75d0 + dcsevl (y*y/4.5d0-1.d0, bi0cs,
     1  nti0)
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'dbesi0',
     +   'abs(x) so big i0 overflows', 2, 2)
c
      dbesi0 = exp(y) * dbsi0e(x)
c
      return
      end
