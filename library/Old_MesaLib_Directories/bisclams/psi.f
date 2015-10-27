      function psi(x)
c***begin prologue  psi
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7c
c***keywords  digamma function,psi function,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the psi (or digamma) function.
c***description
c
c psi(x) calculates the psi (or digamma) function for real argument x.
c psi(x) is the logarithmic derivative of the gamma function of x.
c
c series for psi        on the interval  0.          to  1.00000d+00
c                                        with weighted error   2.03e-17
c                                         log weighted error  16.69
c          significant figures required  16.39
c                                    decimal places required  17.37
c
c series for apsi       on the interval  0.          to  2.50000d-01
c                                        with weighted error   5.54e-17
c                                         log weighted error  16.26
c          significant figures required  14.42
c                                    decimal places required  16.86
c***references  (none)
c***routines called  csevl,inits,r1mach,xerror
c***end prologue  psi
      implicit real*8(a-h,o-z)
      dimension psics(23), apsics(16)
      data psi cs( 1) /   -.0380570808 35217922d0 /
      data psi cs( 2) /    .4914153930 2938713d0 /
      data psi cs( 3) /   -.0568157478 21244730d0 /
      data psi cs( 4) /    .0083578212 25914313d0 /
      data psi cs( 5) /   -.0013332328 57994342d0 /
      data psi cs( 6) /    .0002203132 87069308d0 /
      data psi cs( 7) /   -.0000370402 38178456d0 /
      data psi cs( 8) /    .0000062837 93654854d0 /
      data psi cs( 9) /   -.0000010712 63908506d0 /
      data psi cs(10) /    .0000001831 28394654d0 /
      data psi cs(11) /   -.0000000313 53509361d0 /
      data psi cs(12) /    .0000000053 72808776d0 /
      data psi cs(13) /   -.0000000009 21168141d0 /
      data psi cs(14) /    .0000000001 57981265d0 /
      data psi cs(15) /   -.0000000000 27098646d0 /
      data psi cs(16) /    .0000000000 04648722d0 /
      data psi cs(17) /   -.0000000000 00797527d0 /
      data psi cs(18) /    .0000000000 00136827d0 /
      data psi cs(19) /   -.0000000000 00023475d0 /
      data psi cs(20) /    .0000000000 00004027d0 /
      data psi cs(21) /   -.0000000000 00000691d0 /
      data psi cs(22) /    .0000000000 00000118d0 /
      data psi cs(23) /   -.0000000000 00000020d0 /
      data apsics( 1) /   -.0204749044 678185d0 /
      data apsics( 2) /   -.0101801271 534859d0 /
      data apsics( 3) /    .0000559718 725387d0 /
      data apsics( 4) /   -.0000012917 176570d0 /
      data apsics( 5) /    .0000000572 858606d0 /
      data apsics( 6) /   -.0000000038 213539d0 /
      data apsics( 7) /    .0000000003 397434d0 /
      data apsics( 8) /   -.0000000000 374838d0 /
      data apsics( 9) /    .0000000000 048990d0 /
      data apsics(10) /   -.0000000000 007344d0 /
      data apsics(11) /    .0000000000 001233d0 /
      data apsics(12) /   -.0000000000 000228d0 /
      data apsics(13) /    .0000000000 000045d0 /
      data apsics(14) /   -.0000000000 000009d0 /
      data apsics(15) /    .0000000000 000002d0 /
      data apsics(16) /   -.0000000000 000000d0 /
      data pi     / 3.1415926535 8979324d0/
      data ntpsi, ntapsi, xbig, dxrel /0, 0, 0.d0, 0.d0/
      common/io/inp,iout
c***first executable statement  psi
      if (ntpsi.ne.0) go to 10
      ntpsi = inits (psics, 23, 0.1d0*r1mach(3))
      ntapsi = inits (apsics, 16, 0.1d0*r1mach(3))
c
      xbig = 1.0d0/sqrt(r1mach(3))
      dxrel = sqrt (r1mach(4))
c
 10   y = abs(x)
      if (y.ge.2.0d0) go to 30
c
c psi(x) for -2. .lt. x .lt. 2.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - float(n)
      n = n - 1
      psi = csevl (2.d0*y-1.d0, psics, ntpsi)
      if (n.eq.0) return
c
      n = -n
      if (x.eq.0.d0) call lnkerr ( 'psi     x is 0')
      if (x.lt.0.d0 .and. x+float(n-2).eq.0.d0) 
     1    call lnkerr (  'psi     x is a negative integer')
      if (x.lt.(-0.5d0) .and. abs((x-aint(x-0.5d0))/x).lt.dxrel) 
     1    call lnkerr(  'psi     answer lt half precision because '//
     1                  ' x too near negative integer')
c
      do 20 i=1,n
        psi = psi - 1.0d0/(x+float(i-1))
 20   continue
      return
c
c psi(x) for abs(x) .ge. 2.
c
 30   aux = 0.d0
      if (y.lt.xbig) aux = csevl (8.d0/y**2-1.d0, apsics, ntapsi)
      if (x.lt.0.d0) psi = log(abs(x)) - 0.5d0/x + aux - pi*cot(pi*x)
      if (x.gt.0.d0) psi = log(x) - 0.5d0/x + aux
      return
c
      end
