*deck cuni2
      subroutine cuni2 (z, fnu, kode, n, y, nz, nlast, fnul, tol, elim,
     +   alim)
c***begin prologue  cuni2
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (cuni2-a, zuni2-a)
c***author  amos, d. e., (snl)
c***description
c
c     cuni2 computes i(fnu,z) in the right half plane by means of
c     uniform asymptotic expansion for j(fnu,zn) where zn is z*i
c     or -z*i and zn is in the right half plane also.
c
c     fnul is the smallest order permitted for the asymptotic
c     expansion. nlast=0 means all of the y values were set.
c     nlast.ne.0 is the number left to be computed by another
c     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1.lt.fnul.
c     y(i)=czero for i=nlast+1,n
c
c***see also  cbesi, cbesk
c***routines called  cairy, cuchk, cunhj, cuoik, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cuni2
      complex ai, arg, asum, bsum, cfn, ci, cid, cip, cone, crsc, cscl,
     * csr, css, cy, czero, c1, c2, dai, phi, rz, s1, s2, y, z, zb,
     * zeta1, zeta2, zn, zar
      real aarg, aic, alim, ang, aphi, ascle, ay, bry, car, c2i, c2m,
     * c2r, elim, fn, fnu, fnul, hpi, rs1, sar, tol, yy, r1mach
      integer i, iflag, in, inu, j, k, kode, n, nai, nd, ndai, nlast,
     * nn, nuf, nw, nz, idum
      dimension bry(3), y(n), cip(4), css(3), csr(3), cy(2)
      data czero,cone,ci/(0.0e0,0.0e0),(1.0e0,0.0e0),(0.0e0,1.0e0)/
      data cip(1),cip(2),cip(3),cip(4)/
     1 (1.0e0,0.0e0), (0.0e0,1.0e0), (-1.0e0,0.0e0), (0.0e0,-1.0e0)/
      data hpi, aic  /
     1      1.57079632679489662e+00,     1.265512123484645396e+00/
c***first executable statement  cuni2
      nz = 0
      nd = n
      nlast = 0
c-----------------------------------------------------------------------
c     computed values with exponents between alim and elim in mag-
c     nitude are scaled to keep intermediate arithmetic on scale,
c     exp(alim)=exp(elim)*tol
c-----------------------------------------------------------------------
      cscl = cmplx(1.0e0/tol,0.0e0)
      crsc = cmplx(tol,0.0e0)
      css(1) = cscl
      css(2) = cone
      css(3) = crsc
      csr(1) = crsc
      csr(2) = cone
      csr(3) = cscl
      bry(1) = 1.0e+3*r1mach(1)/tol
      yy = aimag(z)
c-----------------------------------------------------------------------
c     zn is in the right half plane after rotation by ci or -ci
c-----------------------------------------------------------------------
      zn = -z*ci
      zb = z
      cid = -ci
      inu = fnu
      ang = hpi*(fnu-inu)
      car = cos(ang)
      sar = sin(ang)
      c2 = cmplx(car,sar)
      zar = c2
      in = inu + n - 1
      in = mod(in,4)
      c2 = c2*cip(in+1)
      if (yy.gt.0.0e0) go to 10
      zn = conjg(-zn)
      zb = conjg(zb)
      cid = -cid
      c2 = conjg(c2)
   10 continue
c-----------------------------------------------------------------------
c     check for underflow and overflow on first member
c-----------------------------------------------------------------------
      fn = max(fnu,1.0e0)
      call cunhj(zn, fn, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
      if (kode.eq.1) go to 20
      cfn = cmplx(fnu,0.0e0)
      s1 = -zeta1 + cfn*(cfn/(zb+zeta2))
      go to 30
   20 continue
      s1 = -zeta1 + zeta2
   30 continue
      rs1 = real(s1)
      if (abs(rs1).gt.elim) go to 150
   40 continue
      nn = min(2,nd)
      do 90 i=1,nn
        fn = fnu + (nd-i)
        call cunhj(zn, fn, 0, tol, phi, arg, zeta1, zeta2, asum, bsum)
        if (kode.eq.1) go to 50
        cfn = cmplx(fn,0.0e0)
        ay = abs(yy)
        s1 = -zeta1 + cfn*(cfn/(zb+zeta2)) + cmplx(0.0e0,ay)
        go to 60
   50   continue
        s1 = -zeta1 + zeta2
   60   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = real(s1)
        if (abs(rs1).gt.elim) go to 120
        if (i.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 70
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
        aphi = abs(phi)
        aarg = abs(arg)
        rs1 = rs1 + alog(aphi) - 0.25e0*alog(aarg) - aic
        if (abs(rs1).gt.elim) go to 120
        if (i.eq.1) iflag = 1
        if (rs1.lt.0.0e0) go to 70
        if (i.eq.1) iflag = 3
   70   continue
c-----------------------------------------------------------------------
c     scale s1 to keep intermediate arithmetic on scale near
c     exponent extremes
c-----------------------------------------------------------------------
        call cairy(arg, 0, 2, ai, nai, idum)
        call cairy(arg, 1, 2, dai, ndai, idum)
        s2 = phi*(ai*asum+dai*bsum)
        c2r = real(s1)
        c2i = aimag(s1)
        c2m = exp(c2r)*real(css(iflag))
        s1 = cmplx(c2m,0.0e0)*cmplx(cos(c2i),sin(c2i))
        s2 = s2*s1
        if (iflag.ne.1) go to 80
        call cuchk(s2, nw, bry(1), tol)
        if (nw.ne.0) go to 120
   80   continue
        if (yy.le.0.0e0) s2 = conjg(s2)
        j = nd - i + 1
        s2 = s2*c2
        cy(i) = s2
        y(j) = s2*csr(iflag)
        c2 = c2*cid
   90 continue
      if (nd.le.2) go to 110
      rz = cmplx(2.0e0,0.0e0)/z
      bry(2) = 1.0e0/bry(1)
      bry(3) = r1mach(2)
      s1 = cy(1)
      s2 = cy(2)
      c1 = csr(iflag)
      ascle = bry(iflag)
      k = nd - 2
      fn = k
      do 100 i=3,nd
        c2 = s2
        s2 = s1 + cmplx(fnu+fn,0.0e0)*rz*s2
        s1 = c2
        c2 = s2*c1
        y(k) = c2
        k = k - 1
        fn = fn - 1.0e0
        if (iflag.ge.3) go to 100
        c2r = real(c2)
        c2i = aimag(c2)
        c2r = abs(c2r)
        c2i = abs(c2i)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 100
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1*c1
        s2 = c2
        s1 = s1*css(iflag)
        s2 = s2*css(iflag)
        c1 = csr(iflag)
  100 continue
  110 continue
      return
  120 continue
      if (rs1.gt.0.0e0) go to 140
c-----------------------------------------------------------------------
c     set underflow and update parameters
c-----------------------------------------------------------------------
      y(nd) = czero
      nz = nz + 1
      nd = nd - 1
      if (nd.eq.0) go to 110
      call cuoik(z, fnu, kode, 1, nd, y, nuf, tol, elim, alim)
      if (nuf.lt.0) go to 140
      nd = nd - nuf
      nz = nz + nuf
      if (nd.eq.0) go to 110
      fn = fnu + (nd-1)
      if (fn.lt.fnul) go to 130
c      fn = aimag(cid)
c      j = nuf + 1
c      k = mod(j,4) + 1
c      s1 = cip(k)
c      if (fn.lt.0.0e0) s1 = conjg(s1)
c      c2 = c2*s1
      in = inu + nd - 1
      in = mod(in,4) + 1
      c2 = zar*cip(in)
      if (yy.le.0.0e0)c2=conjg(c2)
      go to 40
  130 continue
      nlast = nd
      return
  140 continue
      nz = -1
      return
  150 continue
      if (rs1.gt.0.0e0) go to 140
      nz = n
      do 160 i=1,n
        y(i) = czero
  160 continue
      return
      end
