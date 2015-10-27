*deck casyi
      subroutine casyi (z, fnu, kode, n, y, nz, rl, tol, elim, alim)
c***begin prologue  casyi
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (casyi-a, zasyi-a)
c***author  amos, d. e., (snl)
c***description
c
c     casyi computes the i bessel function for real(z).ge.0.0 by
c     means of the asymptotic expansion for large abs(z) in the
c     region abs(z).gt.max(rl,fnu*fnu/2). nz=0 is a normal return.
c     nz.lt.0 indicates an overflow on kode=1.
c
c***see also  cbesi, cbesk
c***routines called  r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  casyi
      complex ak1, ck, cone, cs1, cs2, cz, czero, dk, ez, p1, rz, s2,
     * y, z
      real aa, acz, aez, ak, alim, arg, arm, atol, az, bb, bk, dfnu,
     * dnu2, elim, fdn, fnu, pi, rl, rtpi, rtr1, s, sgn, sqk, tol, x,
     * yy, r1mach
      integer i, ib, il, inu, j, jl, k, kode, koded, m, n, nn, nz
      dimension y(n)
      data pi, rtpi  /3.14159265358979324e0 , 0.159154943091895336e0 /
      data czero, cone / (0.0e0,0.0e0), (1.0e0,0.0e0) /
c***first executable statement  casyi
      nz = 0
      az = abs(z)
      x = real(z)
      arm = 1.0e+3*r1mach(1)
      rtr1 = sqrt(arm)
      il = min(2,n)
      dfnu = fnu + (n-il)
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      ak1 = cmplx(rtpi,0.0e0)/z
      ak1 = csqrt(ak1)
      cz = z
      if (kode.eq.2) cz = z - cmplx(x,0.0e0)
      acz = real(cz)
      if (abs(acz).gt.elim) go to 80
      dnu2 = dfnu + dfnu
      koded = 1
      if ((abs(acz).gt.alim) .and. (n.gt.2)) go to 10
      koded = 0
      ak1 = ak1*cexp(cz)
   10 continue
      fdn = 0.0e0
      if (dnu2.gt.rtr1) fdn = dnu2*dnu2
      ez = z*cmplx(8.0e0,0.0e0)
c-----------------------------------------------------------------------
c     when z is imaginary, the error test must be made relative to the
c     first reciprocal power since this is the leading term of the
c     expansion for the imaginary part.
c-----------------------------------------------------------------------
      aez = 8.0e0*az
      s = tol/aez
      jl = rl+rl + 2
      yy = aimag(z)
      p1 = czero
      if (yy.eq.0.0e0) go to 20
c-----------------------------------------------------------------------
c     calculate exp(pi*(0.5+fnu+n-il)*i) to minimize losses of
c     significance when fnu or n is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*pi
      inu = inu + n - il
      ak = -sin(arg)
      bk = cos(arg)
      if (yy.lt.0.0e0) bk = -bk
      p1 = cmplx(ak,bk)
      if (mod(inu,2).eq.1) p1 = -p1
   20 continue
      do 50 k=1,il
        sqk = fdn - 1.0e0
        atol = s*abs(sqk)
        sgn = 1.0e0
        cs1 = cone
        cs2 = cone
        ck = cone
        ak = 0.0e0
        aa = 1.0e0
        bb = aez
        dk = ez
        do 30 j=1,jl
          ck = ck*cmplx(sqk,0.0e0)/dk
          cs2 = cs2 + ck
          sgn = -sgn
          cs1 = cs1 + ck*cmplx(sgn,0.0e0)
          dk = dk + ez
          aa = aa*abs(sqk)/bb
          bb = bb + aez
          ak = ak + 8.0e0
          sqk = sqk - ak
          if (aa.le.atol) go to 40
   30   continue
        go to 90
   40   continue
        s2 = cs1
        if (x+x.lt.elim) s2 = s2 + p1*cs2*cexp(-z-z)
        fdn = fdn + 8.0e0*dfnu + 4.0e0
        p1 = -p1
        m = n - il + k
        y(m) = s2*ak1
   50 continue
      if (n.le.2) return
      nn = n
      k = nn - 2
      ak = k
      rz = (cone+cone)/z
      ib = 3
      do 60 i=ib,nn
        y(k) = cmplx(ak+fnu,0.0e0)*rz*y(k+1) + y(k+2)
        ak = ak - 1.0e0
        k = k - 1
   60 continue
      if (koded.eq.0) return
      ck = cexp(cz)
      do 70 i=1,nn
        y(i) = y(i)*ck
   70 continue
      return
   80 continue
      nz = -1
      return
   90 continue
      nz=-2
      return
      end
