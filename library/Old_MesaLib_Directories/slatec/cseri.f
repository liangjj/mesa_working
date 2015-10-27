*deck cseri
      subroutine cseri (z, fnu, kode, n, y, nz, tol, elim, alim)
c***begin prologue  cseri
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (cseri-a, zseri-a)
c***author  amos, d. e., (snl)
c***description
c
c     cseri computes the i bessel function for real(z).ge.0.0 by
c     means of the power series for large abs(z) in the
c     region abs(z).le.2*sqrt(fnu+1). nz=0 is a normal return.
c     nz.gt.0 means that the last nz components were set to zero
c     due to underflow. nz.lt.0 means underflow occurred, but the
c     condition abs(z).le.2*sqrt(fnu+1) was violated and the
c     computation must be completed in another routine with n=n-abs(nz).
c
c***see also  cbesi, cbesk
c***routines called  cuchk, gamln, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cseri
      complex ak1, ck, coef, cone, crsc, cz, czero, hz, rz, s1, s2, w,
     * y, z
      real aa, acz, ak, alim, arm, ascle, atol, az, dfnu, elim, fnu,
     * fnup, rak1, rs, rtr1, s, ss, tol, x, gamln, r1mach
      integer i, ib, idum, iflag, il, k, kode, l, m, n, nn, nw, nz
      dimension y(n), w(2)
      data czero, cone / (0.0e0,0.0e0), (1.0e0,0.0e0) /
c***first executable statement  cseri
      nz = 0
      az = abs(z)
      if (az.eq.0.0e0) go to 150
      x = real(z)
      arm = 1.0e+3*r1mach(1)
      rtr1 = sqrt(arm)
      crsc = cmplx(1.0e0,0.0e0)
      iflag = 0
      if (az.lt.arm) go to 140
      hz = z*cmplx(0.5e0,0.0e0)
      cz = czero
      if (az.gt.rtr1) cz = hz*hz
      acz = abs(cz)
      nn = n
      ck = clog(hz)
   10 continue
      dfnu = fnu + (nn-1)
      fnup = dfnu + 1.0e0
c-----------------------------------------------------------------------
c     underflow test
c-----------------------------------------------------------------------
      ak1 = ck*cmplx(dfnu,0.0e0)
      ak = gamln(fnup,idum)
      ak1 = ak1 - cmplx(ak,0.0e0)
      if (kode.eq.2) ak1 = ak1 - cmplx(x,0.0e0)
      rak1 = real(ak1)
      if (rak1.gt.(-elim)) go to 30
   20 continue
      nz = nz + 1
      y(nn) = czero
      if (acz.gt.dfnu) go to 170
      nn = nn - 1
      if (nn.eq.0) return
      go to 10
   30 continue
      if (rak1.gt.(-alim)) go to 40
      iflag = 1
      ss = 1.0e0/tol
      crsc = cmplx(tol,0.0e0)
      ascle = arm*ss
   40 continue
      ak = aimag(ak1)
      aa = exp(rak1)
      if (iflag.eq.1) aa = aa*ss
      coef = cmplx(aa,0.0e0)*cmplx(cos(ak),sin(ak))
      atol = tol*acz/fnup
      il = min(2,nn)
      do 80 i=1,il
        dfnu = fnu + (nn-i)
        fnup = dfnu + 1.0e0
        s1 = cone
        if (acz.lt.tol*fnup) go to 60
        ak1 = cone
        ak = fnup + 2.0e0
        s = fnup
        aa = 2.0e0
   50   continue
        rs = 1.0e0/s
        ak1 = ak1*cz*cmplx(rs,0.0e0)
        s1 = s1 + ak1
        s = s + ak
        ak = ak + 2.0e0
        aa = aa*acz*rs
        if (aa.gt.atol) go to 50
   60   continue
        m = nn - i + 1
        s2 = s1*coef
        w(i) = s2
        if (iflag.eq.0) go to 70
        call cuchk(s2, nw, ascle, tol)
        if (nw.ne.0) go to 20
   70   continue
        y(m) = s2*crsc
        if (i.ne.il) coef = coef*cmplx(dfnu,0.0e0)/hz
   80 continue
      if (nn.le.2) return
      k = nn - 2
      ak = k
      rz = (cone+cone)/z
      if (iflag.eq.1) go to 110
      ib = 3
   90 continue
      do 100 i=ib,nn
        y(k) = cmplx(ak+fnu,0.0e0)*rz*y(k+1) + y(k+2)
        ak = ak - 1.0e0
        k = k - 1
  100 continue
      return
c-----------------------------------------------------------------------
c     recur backward with scaled values
c-----------------------------------------------------------------------
  110 continue
c-----------------------------------------------------------------------
c     exp(-alim)=exp(-elim)/tol=approx. one precision above the
c     underflow limit = ascle = r1mach(1)*cscl*1.0e+3
c-----------------------------------------------------------------------
      s1 = w(1)
      s2 = w(2)
      do 120 l=3,nn
        ck = s2
        s2 = s1 + cmplx(ak+fnu,0.0e0)*rz*s2
        s1 = ck
        ck = s2*crsc
        y(k) = ck
        ak = ak - 1.0e0
        k = k - 1
        if (abs(ck).gt.ascle) go to 130
  120 continue
      return
  130 continue
      ib = l + 1
      if (ib.gt.nn) return
      go to 90
  140 continue
      nz = n
      if (fnu.eq.0.0e0) nz = nz - 1
  150 continue
      y(1) = czero
      if (fnu.eq.0.0e0) y(1) = cone
      if (n.eq.1) return
      do 160 i=2,n
        y(i) = czero
  160 continue
      return
c-----------------------------------------------------------------------
c     return with nz.lt.0 if abs(z*z/4).gt.fnu+n-nz-1 complete
c     the calculation in cbinu with n=n-abs(nz)
c-----------------------------------------------------------------------
  170 continue
      nz = -nz
      return
      end
