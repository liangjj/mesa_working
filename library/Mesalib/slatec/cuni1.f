*deck cuni1
      subroutine cuni1 (z, fnu, kode, n, y, nz, nlast, fnul, tol, elim,
     +   alim)
c***begin prologue  cuni1
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (cuni1-a, zuni1-a)
c***author  amos, d. e., (snl)
c***description
c
c     cuni1 computes i(fnu,z)  by means of the uniform asymptotic
c     expansion for i(fnu,z) in -pi/3.le.arg z.le.pi/3.
c
c     fnul is the smallest order permitted for the asymptotic
c     expansion. nlast=0 means all of the y values were set.
c     nlast.ne.0 is the number left to be computed by another
c     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1.lt.fnul.
c     y(i)=czero for i=nlast+1,n
c
c***see also  cbesi, cbesk
c***routines called  cuchk, cunik, cuoik, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cuni1
      complex cfn, cone, crsc, cscl, csr, css, cwrk, czero, c1, c2,
     * phi, rz, sum, s1, s2, y, z, zeta1, zeta2, cy
      real alim, aphi, ascle, bry, c2i, c2m, c2r, elim, fn, fnu, fnul,
     * rs1, tol, yy, r1mach
      integer i, iflag, init, k, kode, m, n, nd, nlast, nn, nuf, nw, nz
      dimension bry(3), y(n), cwrk(16), css(3), csr(3), cy(2)
      data czero, cone / (0.0e0,0.0e0), (1.0e0,0.0e0) /
c***first executable statement  cuni1
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
c-----------------------------------------------------------------------
c     check for underflow and overflow on first member
c-----------------------------------------------------------------------
      fn = max(fnu,1.0e0)
      init = 0
      call cunik(z, fn, 1, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
      if (kode.eq.1) go to 10
      cfn = cmplx(fn,0.0e0)
      s1 = -zeta1 + cfn*(cfn/(z+zeta2))
      go to 20
   10 continue
      s1 = -zeta1 + zeta2
   20 continue
      rs1 = real(s1)
      if (abs(rs1).gt.elim) go to 130
   30 continue
      nn = min(2,nd)
      do 80 i=1,nn
        fn = fnu + (nd-i)
        init = 0
        call cunik(z, fn, 1, 0, tol, init, phi, zeta1, zeta2, sum, cwrk)
        if (kode.eq.1) go to 40
        cfn = cmplx(fn,0.0e0)
        yy = aimag(z)
        s1 = -zeta1 + cfn*(cfn/(z+zeta2)) + cmplx(0.0e0,yy)
        go to 50
   40   continue
        s1 = -zeta1 + zeta2
   50   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = real(s1)
        if (abs(rs1).gt.elim) go to 110
        if (i.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 60
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = abs(phi)
        rs1 = rs1 + alog(aphi)
        if (abs(rs1).gt.elim) go to 110
        if (i.eq.1) iflag = 1
        if (rs1.lt.0.0e0) go to 60
        if (i.eq.1) iflag = 3
   60   continue
c-----------------------------------------------------------------------
c     scale s1 if abs(s1).lt.ascle
c-----------------------------------------------------------------------
        s2 = phi*sum
        c2r = real(s1)
        c2i = aimag(s1)
        c2m = exp(c2r)*real(css(iflag))
        s1 = cmplx(c2m,0.0e0)*cmplx(cos(c2i),sin(c2i))
        s2 = s2*s1
        if (iflag.ne.1) go to 70
        call cuchk(s2, nw, bry(1), tol)
        if (nw.ne.0) go to 110
   70   continue
        m = nd - i + 1
        cy(i) = s2
        y(m) = s2*csr(iflag)
   80 continue
      if (nd.le.2) go to 100
      rz = cmplx(2.0e0,0.0e0)/z
      bry(2) = 1.0e0/bry(1)
      bry(3) = r1mach(2)
      s1 = cy(1)
      s2 = cy(2)
      c1 = csr(iflag)
      ascle = bry(iflag)
      k = nd - 2
      fn = k
      do 90 i=3,nd
        c2 = s2
        s2 = s1 + cmplx(fnu+fn,0.0e0)*rz*s2
        s1 = c2
        c2 = s2*c1
        y(k) = c2
        k = k - 1
        fn = fn - 1.0e0
        if (iflag.ge.3) go to 90
        c2r = real(c2)
        c2i = aimag(c2)
        c2r = abs(c2r)
        c2i = abs(c2i)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 90
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1*c1
        s2 = c2
        s1 = s1*css(iflag)
        s2 = s2*css(iflag)
        c1 = csr(iflag)
   90 continue
  100 continue
      return
c-----------------------------------------------------------------------
c     set underflow and update parameters
c-----------------------------------------------------------------------
  110 continue
      if (rs1.gt.0.0e0) go to 120
      y(nd) = czero
      nz = nz + 1
      nd = nd - 1
      if (nd.eq.0) go to 100
      call cuoik(z, fnu, kode, 1, nd, y, nuf, tol, elim, alim)
      if (nuf.lt.0) go to 120
      nd = nd - nuf
      nz = nz + nuf
      if (nd.eq.0) go to 100
      fn = fnu + (nd-1)
      if (fn.ge.fnul) go to 30
      nlast = nd
      return
  120 continue
      nz = -1
      return
  130 continue
      if (rs1.gt.0.0e0) go to 120
      nz = n
      do 140 i=1,n
        y(i) = czero
  140 continue
      return
      end
