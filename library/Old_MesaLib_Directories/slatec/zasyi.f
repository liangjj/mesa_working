*deck zasyi
      subroutine zasyi (zr, zi, fnu, kode, n, yr, yi, nz, rl, tol, elim,
     +   alim)
c***begin prologue  zasyi
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (casyi-a, zasyi-a)
c***author  amos, d. e., (snl)
c***description
c
c     zasyi computes the i bessel function for real(z).ge.0.0 by
c     means of the asymptotic expansion for large abs(z) in the
c     region abs(z).gt.max(rl,fnu*fnu/2). nz=0 is a normal return.
c     nz.lt.0 indicates an overflow on kode=1.
c
c***see also  zbesi, zbesk
c***routines called  d1mach, zabs, zdiv, zexp, zmlt, zsqrt
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zexp and zsqrt to external statement.  (rwc)
c***end prologue  zasyi
c     complex ak1,ck,cone,cs1,cs2,cz,czero,dk,ez,p1,rz,s2,y,z
      double precision aa, aez, ak, ak1i, ak1r, alim, arg, arm, atol,
     * az, bb, bk, cki, ckr, conei, coner, cs1i, cs1r, cs2i, cs2r, czi,
     * czr, dfnu, dki, dkr, dnu2, elim, ezi, ezr, fdn, fnu, pi, p1i,
     * p1r, raz, rl, rtpi, rtr1, rzi, rzr, s, sgn, sqk, sti, str, s2i,
     * s2r, tol, tzi, tzr, yi, yr, zeroi, zeror, zi, zr, d1mach, zabs
      integer i, ib, il, inu, j, jl, k, kode, koded, m, n, nn, nz
      dimension yr(n), yi(n)
      external zabs, zexp, zsqrt
      data pi, rtpi  /3.14159265358979324d0 , 0.159154943091895336d0 /
      data zeror,zeroi,coner,conei / 0.0d0, 0.0d0, 1.0d0, 0.0d0 /
c***first executable statement  zasyi
      nz = 0
      az = zabs(zr,zi)
      arm = 1.0d+3*d1mach(1)
      rtr1 = sqrt(arm)
      il = min(2,n)
      dfnu = fnu + (n-il)
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      raz = 1.0d0/az
      str = zr*raz
      sti = -zi*raz
      ak1r = rtpi*str*raz
      ak1i = rtpi*sti*raz
      call zsqrt(ak1r, ak1i, ak1r, ak1i)
      czr = zr
      czi = zi
      if (kode.ne.2) go to 10
      czr = zeror
      czi = zi
   10 continue
      if (abs(czr).gt.elim) go to 100
      dnu2 = dfnu + dfnu
      koded = 1
      if ((abs(czr).gt.alim) .and. (n.gt.2)) go to 20
      koded = 0
      call zexp(czr, czi, str, sti)
      call zmlt(ak1r, ak1i, str, sti, ak1r, ak1i)
   20 continue
      fdn = 0.0d0
      if (dnu2.gt.rtr1) fdn = dnu2*dnu2
      ezr = zr*8.0d0
      ezi = zi*8.0d0
c-----------------------------------------------------------------------
c     when z is imaginary, the error test must be made relative to the
c     first reciprocal power since this is the leading term of the
c     expansion for the imaginary part.
c-----------------------------------------------------------------------
      aez = 8.0d0*az
      s = tol/aez
      jl = rl+rl + 2
      p1r = zeror
      p1i = zeroi
      if (zi.eq.0.0d0) go to 30
c-----------------------------------------------------------------------
c     calculate exp(pi*(0.5+fnu+n-il)*i) to minimize losses of
c     significance when fnu or n is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*pi
      inu = inu + n - il
      ak = -sin(arg)
      bk = cos(arg)
      if (zi.lt.0.0d0) bk = -bk
      p1r = ak
      p1i = bk
      if (mod(inu,2).eq.0) go to 30
      p1r = -p1r
      p1i = -p1i
   30 continue
      do 70 k=1,il
        sqk = fdn - 1.0d0
        atol = s*abs(sqk)
        sgn = 1.0d0
        cs1r = coner
        cs1i = conei
        cs2r = coner
        cs2i = conei
        ckr = coner
        cki = conei
        ak = 0.0d0
        aa = 1.0d0
        bb = aez
        dkr = ezr
        dki = ezi
        do 40 j=1,jl
          call zdiv(ckr, cki, dkr, dki, str, sti)
          ckr = str*sqk
          cki = sti*sqk
          cs2r = cs2r + ckr
          cs2i = cs2i + cki
          sgn = -sgn
          cs1r = cs1r + ckr*sgn
          cs1i = cs1i + cki*sgn
          dkr = dkr + ezr
          dki = dki + ezi
          aa = aa*abs(sqk)/bb
          bb = bb + aez
          ak = ak + 8.0d0
          sqk = sqk - ak
          if (aa.le.atol) go to 50
   40   continue
        go to 110
   50   continue
        s2r = cs1r
        s2i = cs1i
        if (zr+zr.ge.elim) go to 60
        tzr = zr + zr
        tzi = zi + zi
        call zexp(-tzr, -tzi, str, sti)
        call zmlt(str, sti, p1r, p1i, str, sti)
        call zmlt(str, sti, cs2r, cs2i, str, sti)
        s2r = s2r + str
        s2i = s2i + sti
   60   continue
        fdn = fdn + 8.0d0*dfnu + 4.0d0
        p1r = -p1r
        p1i = -p1i
        m = n - il + k
        yr(m) = s2r*ak1r - s2i*ak1i
        yi(m) = s2r*ak1i + s2i*ak1r
   70 continue
      if (n.le.2) return
      nn = n
      k = nn - 2
      ak = k
      str = zr*raz
      sti = -zi*raz
      rzr = (str+str)*raz
      rzi = (sti+sti)*raz
      ib = 3
      do 80 i=ib,nn
        yr(k) = (ak+fnu)*(rzr*yr(k+1)-rzi*yi(k+1)) + yr(k+2)
        yi(k) = (ak+fnu)*(rzr*yi(k+1)+rzi*yr(k+1)) + yi(k+2)
        ak = ak - 1.0d0
        k = k - 1
   80 continue
      if (koded.eq.0) return
      call zexp(czr, czi, ckr, cki)
      do 90 i=1,nn
        str = yr(i)*ckr - yi(i)*cki
        yi(i) = yr(i)*cki + yi(i)*ckr
        yr(i) = str
   90 continue
      return
  100 continue
      nz = -1
      return
  110 continue
      nz=-2
      return
      end
