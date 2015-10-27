      function cpsi(zin)
c***begin prologue  cpsi
c***date written   780501   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7c
c***keywords  complex,digamma function,psi function,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the psi function of complex argument.
c***description
c
c psi(x) calculates the psi (or digamma) function of x.  psi(x)
c is the logarithmic derivative of the gamma function of x.
c***references  (none)
c***routines called  ccot,r1mach,xerror
c***end prologue  cpsi
      implicit real*8 (a-h,o-z)
      complex*16 cpsi, zin, z, z2inv, corr, ccot
      dimension bern(13)
      data bern( 1) /   .8333333333 3333333 d-1 /
      data bern( 2) /  -.8333333333 3333333 d-2 /
      data bern( 3) /   .3968253968 2539683 d-2 /
      data bern( 4) /  -.4166666666 6666667 d-2 /
      data bern( 5) /   .7575757575 7575758 d-2 /
      data bern( 6) /  -.2109279609 2796093 d-1 /
      data bern( 7) /   .8333333333 3333333 d-1 /
      data bern( 8) /  -.4432598039 2156863 d0 /
      data bern( 9) /   .3053954330 2701197 d1 /
      data bern(10) /  -.2645621212 1212121 d2 /
      data bern(11) /   .2814601449 2753623 d3 /
      data bern(12) /  -.3454885393 7728938 d4 /
      data bern(13) /   .5482758333 3333333 d5 /
      data pi / 3.141592653 589793 d0 /
      data nterm, bound, dxrel, rmin, rbig / 0, 4*0.d0 /
c***first executable statement  cpsi
      if (nterm.ne.0) go to 10
      nterm = -0.30d0*log(r1mach(3))
c maybe bound = n*(0.1d0*eps)**(-1/(2*n-1)) / (pi*exp(1))
      bound = 0.1171d0*float(nterm) *
     1  (0.1d0*r1mach(3))**(-1.d0/(2.d0*float(nterm)-1.d0))
      dxrel = sqrt(r1mach(4))
      rmin = exp (max (log(r1mach(1)), -log(r1mach(2))) + 0.011d0 )
      rbig = 1.0d0/r1mach(3)
c
 10   z = zin
      x = real(z)
      y = imag(z)
      if (y.lt.0.d0) z = conjg(z)
c
      corr = dcmplx(0.d0, 0.d0)
      cabsz = abs(z)
      if (x.ge.0.d0 .and. cabsz.gt.bound) go to 50
      if (x.lt.0.d0 .and. abs(y).gt.bound) go to 50
c
      if (cabsz.lt.bound) go to 20
c
c use the reflection formula for real(z) negative, cabs(z) large, and
c abs(aimag(y)) small.
c
      corr = -pi*ccot(pi*z)
      z = 1.0d0 - z
      go to 50
c
c use the recursion relation for cabs(z) small.
c
 20   if (cabsz.lt.rmin) then
          call lnkerr ( 'cpsi called with z so near 0 that cpsi '// 
     1                  'overflows')
      endif
c
      if (x.ge.(-0.5d0) .or. abs(y).gt.dxrel) go to 30
      if (abs((z-aint(x-0.5d0))/x).lt.dxrel) call lnkerr ( 'cpsi  answe
     1r lt half precision because z too near negative integer')
      if (y.eq.0.d0 .and. x.eq.aint(x)) call lnkerr ( 'cpsi   z is a neg
     1ative integer')
c
 30   n = sqrt(bound**2-y**2) - x + 1.0d0
      do 40 i=1,n
        corr = corr - 1.0d0/z
        z = z + 1.0d0
 40   continue
c
c now evaluate the asymptotic series for suitably large z.
c
 50   if (cabsz.gt.rbig) cpsi = log(z) + corr
      if (cabsz.gt.rbig) go to 70
c
      cpsi = dcmplx(0.d0, 0.d0)
      z2inv = 1.d0/z**2
      do 60 i=1,nterm
        ndx = nterm + 1 - i
        cpsi = bern(ndx) + z2inv*cpsi
 60   continue
      cpsi = log(z) - 0.5d0/z - cpsi*z2inv + corr
c
 70   if (y.lt.0.d0) cpsi = conjg(cpsi)
c
      return
      end

