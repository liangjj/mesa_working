*deck zseri
      subroutine zseri (zr, zi, fnu, kode, n, yr, yi, nz, tol, elim,
     +   alim)
c***begin prologue  zseri
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cseri-a, zseri-a)
c***author  amos, d. e., (snl)
c***description
c
c     zseri computes the i bessel function for real(z).ge.0.0 by
c     means of the power series for large abs(z) in the
c     region abs(z).le.2*sqrt(fnu+1). nz=0 is a normal return.
c     nz.gt.0 means that the last nz components were set to zero
c     due to underflow. nz.lt.0 means underflow occurred, but the
c     condition abs(z).le.2*sqrt(fnu+1) was violated and the
c     computation must be completed in another routine with n=n-abs(nz).
c
c***see also  zbesi, zbesk
c***routines called  d1mach, dgamln, zabs, zdiv, zlog, zmlt, zuchk
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zlog to external statement.  (rwc)
c***end prologue  zseri
c     complex ak1,ck,coef,cone,crsc,cscl,cz,czero,hz,rz,s1,s2,y,z
      double precision aa, acz, ak, ak1i, ak1r, alim, arm, ascle, atol,
     * az, cki, ckr, coefi, coefr, conei, coner, crscr, czi, czr, dfnu,
     * elim, fnu, fnup, hzi, hzr, raz, rs, rtr1, rzi, rzr, s, ss, sti,
     * str, s1i, s1r, s2i, s2r, tol, yi, yr, wi, wr, zeroi, zeror, zi,
     * zr, dgamln, d1mach, zabs
      integer i, ib, idum, iflag, il, k, kode, l, m, n, nn, nz, nw
      dimension yr(n), yi(n), wr(2), wi(2)
      external zabs, zlog
      data zeror,zeroi,coner,conei / 0.0d0, 0.0d0, 1.0d0, 0.0d0 /
c***first executable statement  zseri
      nz = 0
      az = zabs(zr,zi)
      if (az.eq.0.0d0) go to 160
      arm = 1.0d+3*d1mach(1)
      rtr1 = sqrt(arm)
      crscr = 1.0d0
      iflag = 0
      if (az.lt.arm) go to 150
      hzr = 0.5d0*zr
      hzi = 0.5d0*zi
      czr = zeror
      czi = zeroi
      if (az.le.rtr1) go to 10
      call zmlt(hzr, hzi, hzr, hzi, czr, czi)
   10 continue
      acz = zabs(czr,czi)
      nn = n
      call zlog(hzr, hzi, ckr, cki, idum)
   20 continue
      dfnu = fnu + (nn-1)
      fnup = dfnu + 1.0d0
c-----------------------------------------------------------------------
c     underflow test
c-----------------------------------------------------------------------
      ak1r = ckr*dfnu
      ak1i = cki*dfnu
      ak = dgamln(fnup,idum)
      ak1r = ak1r - ak
      if (kode.eq.2) ak1r = ak1r - zr
      if (ak1r.gt.(-elim)) go to 40
   30 continue
      nz = nz + 1
      yr(nn) = zeror
      yi(nn) = zeroi
      if (acz.gt.dfnu) go to 190
      nn = nn - 1
      if (nn.eq.0) return
      go to 20
   40 continue
      if (ak1r.gt.(-alim)) go to 50
      iflag = 1
      ss = 1.0d0/tol
      crscr = tol
      ascle = arm*ss
   50 continue
      aa = exp(ak1r)
      if (iflag.eq.1) aa = aa*ss
      coefr = aa*cos(ak1i)
      coefi = aa*sin(ak1i)
      atol = tol*acz/fnup
      il = min(2,nn)
      do 90 i=1,il
        dfnu = fnu + (nn-i)
        fnup = dfnu + 1.0d0
        s1r = coner
        s1i = conei
        if (acz.lt.tol*fnup) go to 70
        ak1r = coner
        ak1i = conei
        ak = fnup + 2.0d0
        s = fnup
        aa = 2.0d0
   60   continue
        rs = 1.0d0/s
        str = ak1r*czr - ak1i*czi
        sti = ak1r*czi + ak1i*czr
        ak1r = str*rs
        ak1i = sti*rs
        s1r = s1r + ak1r
        s1i = s1i + ak1i
        s = s + ak
        ak = ak + 2.0d0
        aa = aa*acz*rs
        if (aa.gt.atol) go to 60
   70   continue
        s2r = s1r*coefr - s1i*coefi
        s2i = s1r*coefi + s1i*coefr
        wr(i) = s2r
        wi(i) = s2i
        if (iflag.eq.0) go to 80
        call zuchk(s2r, s2i, nw, ascle, tol)
        if (nw.ne.0) go to 30
   80   continue
        m = nn - i + 1
        yr(m) = s2r*crscr
        yi(m) = s2i*crscr
        if (i.eq.il) go to 90
        call zdiv(coefr, coefi, hzr, hzi, str, sti)
        coefr = str*dfnu
        coefi = sti*dfnu
   90 continue
      if (nn.le.2) return
      k = nn - 2
      ak = k
      raz = 1.0d0/az
      str = zr*raz
      sti = -zi*raz
      rzr = (str+str)*raz
      rzi = (sti+sti)*raz
      if (iflag.eq.1) go to 120
      ib = 3
  100 continue
      do 110 i=ib,nn
        yr(k) = (ak+fnu)*(rzr*yr(k+1)-rzi*yi(k+1)) + yr(k+2)
        yi(k) = (ak+fnu)*(rzr*yi(k+1)+rzi*yr(k+1)) + yi(k+2)
        ak = ak - 1.0d0
        k = k - 1
  110 continue
      return
c-----------------------------------------------------------------------
c     recur backward with scaled values
c-----------------------------------------------------------------------
  120 continue
c-----------------------------------------------------------------------
c     exp(-alim)=exp(-elim)/tol=approx. one precision above the
c     underflow limit = ascle = d1mach(1)*ss*1.0d+3
c-----------------------------------------------------------------------
      s1r = wr(1)
      s1i = wi(1)
      s2r = wr(2)
      s2i = wi(2)
      do 130 l=3,nn
        ckr = s2r
        cki = s2i
        s2r = s1r + (ak+fnu)*(rzr*ckr-rzi*cki)
        s2i = s1i + (ak+fnu)*(rzr*cki+rzi*ckr)
        s1r = ckr
        s1i = cki
        ckr = s2r*crscr
        cki = s2i*crscr
        yr(k) = ckr
        yi(k) = cki
        ak = ak - 1.0d0
        k = k - 1
        if (zabs(ckr,cki).gt.ascle) go to 140
  130 continue
      return
  140 continue
      ib = l + 1
      if (ib.gt.nn) return
      go to 100
  150 continue
      nz = n
      if (fnu.eq.0.0d0) nz = nz - 1
  160 continue
      yr(1) = zeror
      yi(1) = zeroi
      if (fnu.ne.0.0d0) go to 170
      yr(1) = coner
      yi(1) = conei
  170 continue
      if (n.eq.1) return
      do 180 i=2,n
        yr(i) = zeror
        yi(i) = zeroi
  180 continue
      return
c-----------------------------------------------------------------------
c     return with nz.lt.0 if abs(z*z/4).gt.fnu+n-nz-1 complete
c     the calculation in cbinu with n=n-abs(nz)
c-----------------------------------------------------------------------
  190 continue
      nz = -nz
      return
      end
