      function c9lgmc(zin)
c***begin prologue  c9lgmc
c***date written   780401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7a
c***keywords  complete gamma function,complex,correction term,
c             gamma function,logarithm,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the log gamma correction term for most z so that
c            clog(cgamma(z)) = 0.5*alog(2.*pi) + (z-0.5)*clog(z) - z
c            + c9lgmc(z)
c***description
c
c compute the log gamma correction term for large cabs(z) when real(z)
c .ge. 0.0 and for large abs(aimag(y)) when real(z) .lt. 0.0.  we find
c c9lgmc so that
c   clog((z)) = 0.5*alog(2.*pi) + (z-0.5)*clog(z) - z + c9lgmc(z)
c***references  (none)
c***routines called  r1mach,xerror
c***end prologue  c9lgmc
      implicit real*8(a-h,o-z)
      complex*16 c9lgmc, zin, z, z2inv
      dimension bern(11)
      data bern( 1) /    .08333333333 3333333d0   /
      data bern( 2) /   -.002777777777 7777778d0  /
      data bern( 3) /    .0007936507936 5079365d0 /
      data bern( 4) /   -.0005952380952 3809524d0 /
      data bern( 5) /    .0008417508417 5084175d0 /
      data bern( 6) /   -.001917526917 5269175d0  /
      data bern( 7) /    .006410256410 2564103d0  /
      data bern( 8) /   -.02955065359 4771242d0   /
      data bern( 9) /    .1796443723 6883057d0    /
      data bern(10) /  -1.392432216 9059011d0     /
      data bern(11) /  13.40286404 4168392d0      /
      data nterm, bound, xbig, xmax / 0, 3*0.0d0 /
c***first executable statement  c9lgmc
      if (nterm.ne.0) go to 10
c
      nterm = -0.30d0*log(r1mach(3))
      bound = 0.1170d0*float(nterm)*
     1  (0.1d0*r1mach(3))**(-1.d0/(2.d0*float(nterm)-1.d0))
      xbig = 1.0d0/sqrt(r1mach(3))
      xmax = exp (min(log(r1mach(2)/12.0d0), -log(12.d0*r1mach(1))))
c
 10   z = zin
      x = real (z)
      y = imag(z)
      cabsz = abs(z)
c
      if (x.lt.0.0d0 .and. abs(y).lt.bound) call lnkerr ('c9lgmc  c9lgmc
     1 not valid for negative real(z) and small abs(imag(z))')
      if (cabsz.lt.bound) call lnkerr ( 'c9lgmc  c9lgmc not valid for sm
     1all abs(z)')
c
      if (cabsz.ge.xmax) go to 50
c
      if (cabsz.ge.xbig) c9lgmc = 1.0d0/(12.0d0*z)
      if (cabsz.ge.xbig) return
c
      z2inv = 1.0d0/z**2
      c9lgmc = (0.0d0, 0.0d0)
      do 40 i=1,nterm
        ndx = nterm + 1 - i
        c9lgmc = bern(ndx) + c9lgmc*z2inv
 40   continue
c
      c9lgmc = c9lgmc/z
      return
c
 50   c9lgmc = (0.0d0, 0.0d0)
      call lnkerr ( 'c9lgmc  z so big c9lgmc underflows')
      return
c
      end








