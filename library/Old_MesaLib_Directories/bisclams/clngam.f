      function clngam(zin)
c***begin prologue  clngam
c***date written   780401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7a
c***keywords  absolute value,complete gamma function,complex,
c             gamma function,logarithm,special function
c***author  fullerton, w., (lanl)
c***purpose  clngam computes the natural log of the complex valued gamma
c            function at zin, where zin is a complex number.
c***description
c
c clngam computes the natural log of the complex valued gamma function
c at zin, where zin is a complex number.  this is a preliminary version,
c which is not accurate.
c***references  (none)
c***routines called  c9lgmc,carg,clnrel,r1mach,xerror
c***end prologue  clngam
      implicit real*8(a-h,o-z)
      complex*16 zin, z, corr, clngam, clnrel, c9lgmc
      data pi / 3.1415926535 8979324d0 /
      data sq2pil / 0.9189385332 0467274d0 /
      data bound, dxrel / 2*0.0d0 /
c***first executable statement  clngam
      if (bound.ne.0.d0) go to 10
      n = -0.30d0*log(r1mach(3))
c bound = n*(0.1d0*eps)**(-1/(2*n-1))/(pi*exp(1))
      bound = 0.1171d0*float(n)*(0.1d0*r1mach(3))**
     1                          (-1.d0/(2.d0*float(n)-1.d0))
      dxrel = sqrt (r1mach(4))
c
 10   z = zin
      x = real(zin)
      y = imag(zin)
c
      corr = (0.0d0, 0.0d0)
      cabsz = abs(z)
      if (x.ge.0.0d0 .and. cabsz.gt.bound) go to 50
      if (x.lt.0.0d0 .and. abs(y).gt.bound) go to 50
c
      if (cabsz.lt.bound) go to 20
c
c use the reflection formula for real(z) negative, cabs(z) large, and
c abs(aimag(y)) small.
c
      if (y.gt.0.0d0) z = conjg (z)
      corr = exp (-dcmplx(0.0d0,2.0d0*pi)*z)
      if (real(corr).eq.1.0d0 .and. imag(corr).eq.0.0d0)
     1    call lnkerr ('clngam  z is a negative integer')
c
      clngam = sq2pil + 1.0d0 - dcmplx(0.0d0,pi)*(z-0.5d0) - 
     1         clnrel(-corr) + (z-0.5d0)*log(1.0d0-z) - 
     2         z - c9lgmc(1.0d0-z)
      if (y.gt.0.0d0) clngam = conjg (clngam)
      return
c
c use the recursion relation for cabs(z) small.
c
 20   if (x.ge.(-0.5d0) .or. abs(y).gt.dxrel) go to 30
      if (abs((z-int(x-0.5d0))/x).lt.dxrel) call lnkerr ( 'clngam  answe
     1r lt half precision because z too near negative integer')
c
 30   n = sqrt (bound**2 - y**2) - x + 1.0d0
      argsum = 0.0d0
      corr = (1.0d0, 0.0d0)
      do 40 i=1,n
        argsum = argsum + carg(z)
        corr = z*corr
        z = 1.0d0 + z
 40   continue
c
      if (real(corr).eq.0.0d0 .and. imag(corr).eq.0.0d0) call lnkerr('cl
     1ngam  z is a negative integer')
      corr = -dcmplx (log(abs(corr)), argsum)
c
c use stirling-s approximation for large z.
c
 50   clngam = sq2pil + (z-0.5d0)*log(z) - z + corr + c9lgmc(z)
      return
c
      end
