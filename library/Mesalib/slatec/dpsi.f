*deck dpsi
      double precision function dpsi (x)
c***begin prologue  dpsi
c***purpose  compute the psi (or digamma) function.
c***library   slatec (fnlib)
c***category  c7c
c***type      double precision (psi-s, dpsi-d, cpsi-c)
c***keywords  digamma function, fnlib, psi function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dpsi calculates the double precision psi (or digamma) function for
c double precision argument x.  psi(x) is the logarithmic derivative
c of the gamma function of x.
c
c series for psi        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.93
c                                    decimal places required  32.05
c
c
c series for apsi       on the interval  0.          to  1.00000e-02
c                                        with weighted error   7.75e-33
c                                         log weighted error  32.11
c                               significant figures required  28.88
c                                    decimal places required  32.71
c
c***references  (none)
c***routines called  d1mach, dcot, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c   920618  removed space from variable name.  (rwc, wrb)
c***end prologue  dpsi
      double precision x, psics(42), apsics(16), aux, dxrel, pi, xbig,
     1  y, dcot, dcsevl, d1mach
      logical first
      external dcot
      save psics, apsics, pi, ntpsi, ntapsi, xbig, dxrel, first
      data psics(  1) / -.3805708083 5217921520 4376776670 39 d-1     /
      data psics(  2) / +.4914153930 2938712748 2046996542 77 d+0     /
      data psics(  3) / -.5681574782 1244730242 8920647340 81 d-1     /
      data psics(  4) / +.8357821225 9143131362 7756507478 62 d-2     /
      data psics(  5) / -.1333232857 9943425998 0792741723 93 d-2     /
      data psics(  6) / +.2203132870 6930824892 8723979795 21 d-3     /
      data psics(  7) / -.3704023817 8456883592 8890869492 29 d-4     /
      data psics(  8) / +.6283793654 8549898933 6514187176 90 d-5     /
      data psics(  9) / -.1071263908 5061849855 2835417470 74 d-5     /
      data psics( 10) / +.1831283946 5484165805 7315898103 78 d-6     /
      data psics( 11) / -.3135350936 1808509869 0057797968 85 d-7     /
      data psics( 12) / +.5372808776 2007766260 4719191436 15 d-8     /
      data psics( 13) / -.9211681415 9784275717 8806326247 30 d-9     /
      data psics( 14) / +.1579812652 1481822782 2528840328 23 d-9     /
      data psics( 15) / -.2709864613 2380443065 4405894097 07 d-10    /
      data psics( 16) / +.4648722859 9096834872 9473195295 49 d-11    /
      data psics( 17) / -.7975272563 8303689726 5047977727 37 d-12    /
      data psics( 18) / +.1368272385 7476992249 2510538928 38 d-12    /
      data psics( 19) / -.2347515606 0658972717 3206779807 19 d-13    /
      data psics( 20) / +.4027630715 5603541107 9079250062 81 d-14    /
      data psics( 21) / -.6910251853 1179037846 5474229747 71 d-15    /
      data psics( 22) / +.1185604713 8863349552 9291395257 68 d-15    /
      data psics( 23) / -.2034168961 6261559308 1542104842 23 d-16    /
      data psics( 24) / +.3490074968 6463043850 3742329323 51 d-17    /
      data psics( 25) / -.5988014693 4976711003 0110813934 93 d-18    /
      data psics( 26) / +.1027380162 8080588258 3980057122 13 d-18    /
      data psics( 27) / -.1762704942 4561071368 3592601053 86 d-19    /
      data psics( 28) / +.3024322801 8156920457 4540354901 33 d-20    /
      data psics( 29) / -.5188916830 2092313774 2860888746 66 d-21    /
      data psics( 30) / +.8902773034 5845713905 0058874879 99 d-22    /
      data psics( 31) / -.1527474289 9426728392 8949719040 00 d-22    /
      data psics( 32) / +.2620731479 8962083136 3583180799 99 d-23    /
      data psics( 33) / -.4496464273 8220696772 5983880533 33 d-24    /
      data psics( 34) / +.7714712959 6345107028 9193642666 66 d-25    /
      data psics( 35) / -.1323635476 1887702968 1026389333 33 d-25    /
      data psics( 36) / +.2270999436 2408300091 2773119999 99 d-26    /
      data psics( 37) / -.3896419021 5374115954 4913919999 99 d-27    /
      data psics( 38) / +.6685198138 8855302310 6798933333 33 d-28    /
      data psics( 39) / -.1146998665 4920864872 5299199999 99 d-28    /
      data psics( 40) / +.1967938588 6541405920 5154133333 33 d-29    /
      data psics( 41) / -.3376448818 9750979801 9072000000 00 d-30    /
      data psics( 42) / +.5793070319 3214159246 6773333333 33 d-31    /
      data apsics(  1) / -.8327107910 6929076017 4456932269 d-3        /
      data apsics(  2) / -.4162518421 9273935282 1627121990 d-3        /
      data apsics(  3) / +.1034315609 7874129117 4463193961 d-6        /
      data apsics(  4) / -.1214681841 3590415298 7299556365 d-9        /
      data apsics(  5) / +.3113694319 9835615552 1240278178 d-12       /
      data apsics(  6) / -.1364613371 9317704177 6516100945 d-14       /
      data apsics(  7) / +.9020517513 1541656513 0837974000 d-17       /
      data apsics(  8) / -.8315429974 2159146482 9933635466 d-19       /
      data apsics(  9) / +.1012242570 7390725418 8479482666 d-20       /
      data apsics( 10) / -.1562702494 3562250762 0478933333 d-22       /
      data apsics( 11) / +.2965427168 0890389613 3226666666 d-24       /
      data apsics( 12) / -.6746868867 6570216374 1866666666 d-26       /
      data apsics( 13) / +.1803453116 9718990421 3333333333 d-27       /
      data apsics( 14) / -.5569016182 4598360746 6666666666 d-29       /
      data apsics( 15) / +.1958679226 0773625173 3333333333 d-30       /
      data apsics( 16) / -.7751958925 2333568000 0000000000 d-32       /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data first /.true./
c***first executable statement  dpsi
      if (first) then
         ntpsi = initds (psics, 42, 0.1*real(d1mach(3)) )
         ntapsi = initds (apsics, 16, 0.1*real(d1mach(3)) )
c
         xbig = 1.0d0/sqrt(d1mach(3))
         dxrel = sqrt(d1mach(4))
      endif
      first = .false.
c
      y = abs(x)
c
      if (y.gt.10.0d0) go to 50
c
c dpsi(x) for abs(x) .le. 2
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - n
      n = n - 1
      dpsi = dcsevl (2.d0*y-1.d0, psics, ntpsi)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
      n = -n
      if (x .eq. 0.d0) call xermsg ('slatec', 'dpsi', 'x is 0', 2, 2)
      if (x .lt. 0.d0 .and. x+n-2 .eq. 0.d0) call xermsg ('slatec',
     +   'dpsi', 'x is a negative integer', 3, 2)
      if (x .lt. (-0.5d0) .and. abs((x-aint(x-0.5d0))/x) .lt. dxrel)
     +   call xermsg ('slatec', 'dpsi',
     +   'answer lt half precision because x too near negative integer',
     +   1, 1)
c
      do 20 i=1,n
        dpsi = dpsi - 1.d0/(x+i-1)
 20   continue
      return
c
c dpsi(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dpsi = dpsi + 1.0d0/(y+i)
 40   continue
      return
c
c dpsi(x) for abs(x) .gt. 10.0
c
 50   aux = 0.d0
      if (y.lt.xbig) aux = dcsevl (2.d0*(10.d0/y)**2-1.d0, apsics,
     1  ntapsi)
c
      if (x.lt.0.d0) dpsi = log(abs(x)) - 0.5d0/x + aux
     1  - pi*dcot(pi*x)
      if (x.gt.0.d0) dpsi = log(x) - 0.5d0/x + aux
      return
c
      end
