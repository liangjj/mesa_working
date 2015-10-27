*deck psi
      function psi (x)
c***begin prologue  psi
c***purpose  compute the psi (or digamma) function.
c***library   slatec (fnlib)
c***category  c7c
c***type      single precision (psi-s, dpsi-d, cpsi-c)
c***keywords  digamma function, fnlib, psi function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c psi(x) calculates the psi (or digamma) function for real argument x.
c psi(x) is the logarithmic derivative of the gamma function of x.
c
c series for psi        on the interval  0.          to  1.00000d+00
c                                        with weighted error   2.03e-17
c                                         log weighted error  16.69
c                               significant figures required  16.39
c                                    decimal places required  17.37
c
c series for apsi       on the interval  0.          to  2.50000d-01
c                                        with weighted error   5.54e-17
c                                         log weighted error  16.26
c                               significant figures required  14.42
c                                    decimal places required  16.86
c
c***references  (none)
c***routines called  cot, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  psi
      dimension psics(23), apsics(16)
      logical first
      external cot
      save psics, apsics, pi, ntpsi, ntapsi, xbig, dxrel, first
      data psics( 1) /   -.0380570808 35217922e0 /
      data psics( 2) /    .4914153930 2938713e0 /
      data psics( 3) /   -.0568157478 21244730e0 /
      data psics( 4) /    .0083578212 25914313e0 /
      data psics( 5) /   -.0013332328 57994342e0 /
      data psics( 6) /    .0002203132 87069308e0 /
      data psics( 7) /   -.0000370402 38178456e0 /
      data psics( 8) /    .0000062837 93654854e0 /
      data psics( 9) /   -.0000010712 63908506e0 /
      data psics(10) /    .0000001831 28394654e0 /
      data psics(11) /   -.0000000313 53509361e0 /
      data psics(12) /    .0000000053 72808776e0 /
      data psics(13) /   -.0000000009 21168141e0 /
      data psics(14) /    .0000000001 57981265e0 /
      data psics(15) /   -.0000000000 27098646e0 /
      data psics(16) /    .0000000000 04648722e0 /
      data psics(17) /   -.0000000000 00797527e0 /
      data psics(18) /    .0000000000 00136827e0 /
      data psics(19) /   -.0000000000 00023475e0 /
      data psics(20) /    .0000000000 00004027e0 /
      data psics(21) /   -.0000000000 00000691e0 /
      data psics(22) /    .0000000000 00000118e0 /
      data psics(23) /   -.0000000000 00000020e0 /
      data apsics( 1) /   -.0204749044 678185e0 /
      data apsics( 2) /   -.0101801271 534859e0 /
      data apsics( 3) /    .0000559718 725387e0 /
      data apsics( 4) /   -.0000012917 176570e0 /
      data apsics( 5) /    .0000000572 858606e0 /
      data apsics( 6) /   -.0000000038 213539e0 /
      data apsics( 7) /    .0000000003 397434e0 /
      data apsics( 8) /   -.0000000000 374838e0 /
      data apsics( 9) /    .0000000000 048990e0 /
      data apsics(10) /   -.0000000000 007344e0 /
      data apsics(11) /    .0000000000 001233e0 /
      data apsics(12) /   -.0000000000 000228e0 /
      data apsics(13) /    .0000000000 000045e0 /
      data apsics(14) /   -.0000000000 000009e0 /
      data apsics(15) /    .0000000000 000002e0 /
      data apsics(16) /   -.0000000000 000000e0 /
      data pi     / 3.1415926535 8979324e0/
      data first /.true./
c***first executable statement  psi
      if (first) then
         ntpsi = inits (psics, 23, 0.1*r1mach(3))
         ntapsi = inits (apsics, 16, 0.1*r1mach(3))
c
         xbig = 1.0/sqrt(r1mach(3))
         dxrel = sqrt (r1mach(4))
      endif
      first = .false.
c
      y = abs(x)
      if (y.ge.2.0) go to 30
c
c psi(x) for -2. .lt. x .lt. 2.
c
      n = x
      if (x.lt.0.) n = n - 1
      y = x - n
      n = n - 1
      psi = csevl (2.*y-1., psics, ntpsi)
      if (n.eq.0) return
c
      n = -n
      if (x .eq. 0.) call xermsg ('slatec', 'psi', 'x is 0', 2, 2)
      if (x .lt. 0. .and. x+n-2 .eq. 0.) call xermsg ('slatec', 'psi',
     +   'x is a negative integer', 3, 2)
      if (x .lt. (-0.5) .and. abs((x-aint(x-0.5))/x) .lt. dxrel)
     +   call xermsg ('slatec', 'psi',
     +   'answer lt half precision because x too near negative integer',
     +   1, 1)
c
      do 20 i=1,n
        psi = psi - 1.0/(x+i-1)
 20   continue
      return
c
c psi(x) for abs(x) .ge. 2.
c
 30   aux = 0.
      if (y.lt.xbig) aux = csevl (8./y**2-1., apsics, ntapsi)
      if (x.lt.0.) psi = log(abs(x)) - 0.5/x + aux - pi*cot(pi*x)
      if (x.gt.0.) psi = log(x) - 0.5/x + aux
      return
c
      end
