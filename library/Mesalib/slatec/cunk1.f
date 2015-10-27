*deck cunk1
      subroutine cunk1 (z, fnu, kode, mr, n, y, nz, tol, elim, alim)
c***begin prologue  cunk1
c***subsidiary
c***purpose  subsidiary to cbesk
c***library   slatec
c***type      all (cunk1-a, zunk1-a)
c***author  amos, d. e., (snl)
c***description
c
c     cunk1 computes k(fnu,z) and its analytic continuation from the
c     right half plane to the left half plane by means of the
c     uniform asymptotic expansion.
c     mr indicates the direction of rotation for analytic continuation.
c     nz=-1 means an overflow will occur
c
c***see also  cbesk
c***routines called  cs1s2, cuchk, cunik, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cunk1
      complex cfn, ck, cone, crsc, cs, cscl, csgn, cspn, csr, css,
     * cwrk, cy, czero, c1, c2, phi,  rz, sum,  s1, s2, y, z,
     * zeta1,  zeta2,  zr, phid, zeta1d, zeta2d, sumd
      real alim, ang, aphi, asc, ascle, bry, cpn, c2i, c2m, c2r, elim,
     * fmr, fn, fnf, fnu, pi, rs1, sgn, spn, tol, x, r1mach
      integer i, ib, iflag, ifn, il, init, inu, iuf, k, kdflg, kflag,
     * kk, kode, mr, n, nw, nz, j, ipard, initd, ic, m
      dimension bry(3), init(2), y(n), sum(2), phi(2), zeta1(2),
     * zeta2(2), cy(2), cwrk(16,3), css(3), csr(3)
      data czero, cone / (0.0e0,0.0e0) , (1.0e0,0.0e0) /
      data pi / 3.14159265358979324e0 /
c***first executable statement  cunk1
      kdflg = 1
      nz = 0
c-----------------------------------------------------------------------
c     exp(-alim)=exp(-elim)/tol=approx. one precision greater than
c     the underflow limit
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
      bry(2) = 1.0e0/bry(1)
      bry(3) = r1mach(2)
      x = real(z)
      zr = z
      if (x.lt.0.0e0) zr = -z
      j=2
      do 70 i=1,n
c-----------------------------------------------------------------------
c     j flip flops between 1 and 2 in j = 3 - j
c-----------------------------------------------------------------------
        j = 3 - j
        fn = fnu + (i-1)
        init(j) = 0
        call cunik(zr, fn, 2, 0, tol, init(j), phi(j), zeta1(j),
     *   zeta2(j), sum(j), cwrk(1,j))
        if (kode.eq.1) go to 20
        cfn = cmplx(fn,0.0e0)
        s1 = zeta1(j) - cfn*(cfn/(zr+zeta2(j)))
        go to 30
   20   continue
        s1 = zeta1(j) - zeta2(j)
   30   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = real(s1)
        if (abs(rs1).gt.elim) go to 60
        if (kdflg.eq.1) kflag = 2
        if (abs(rs1).lt.alim) go to 40
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = abs(phi(j))
        rs1 = rs1 + alog(aphi)
        if (abs(rs1).gt.elim) go to 60
        if (kdflg.eq.1) kflag = 1
        if (rs1.lt.0.0e0) go to 40
        if (kdflg.eq.1) kflag = 3
   40   continue
c-----------------------------------------------------------------------
c     scale s1 to keep intermediate arithmetic on scale near
c     exponent extremes
c-----------------------------------------------------------------------
        s2 = phi(j)*sum(j)
        c2r = real(s1)
        c2i = aimag(s1)
        c2m = exp(c2r)*real(css(kflag))
        s1 = cmplx(c2m,0.0e0)*cmplx(cos(c2i),sin(c2i))
        s2 = s2*s1
        if (kflag.ne.1) go to 50
        call cuchk(s2, nw, bry(1), tol)
        if (nw.ne.0) go to 60
   50   continue
        cy(kdflg) = s2
        y(i) = s2*csr(kflag)
        if (kdflg.eq.2) go to 75
        kdflg = 2
        go to 70
   60   continue
        if (rs1.gt.0.0e0) go to 290
c-----------------------------------------------------------------------
c     for x.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
        if (x.lt.0.0e0) go to 290
        kdflg = 1
        y(i) = czero
        nz=nz+1
        if (i.eq.1) go to 70
        if (y(i-1).eq.czero) go to 70
        y(i-1) = czero
        nz=nz+1
   70 continue
      i=n
   75 continue
      rz = cmplx(2.0e0,0.0e0)/zr
      ck = cmplx(fn,0.0e0)*rz
      ib = i+1
      if (n.lt.ib) go to 160
c-----------------------------------------------------------------------
c     test last member for underflow and overflow, set sequence to zero
c     on underflow
c-----------------------------------------------------------------------
      fn = fnu+(n-1)
      ipard = 1
      if (mr.ne.0) ipard = 0
      initd = 0
      call cunik(zr,fn,2,ipard,tol,initd,phid,zeta1d,zeta2d,sumd,
     *cwrk(1,3))
      if (kode.eq.1) go to 80
      cfn=cmplx(fn,0.0e0)
      s1=zeta1d-cfn*(cfn/(zr+zeta2d))
      go to 90
   80 continue
      s1=zeta1d-zeta2d
   90 continue
      rs1=real(s1)
      if (abs(rs1).gt.elim) go to 95
      if (abs(rs1).lt.alim) go to 100
c-----------------------------------------------------------------------
c     refine estimate and test
c-----------------------------------------------------------------------
      aphi=abs(phid)
      rs1=rs1+alog(aphi)
      if (abs(rs1).lt.elim) go to 100
   95 continue
      if (rs1.gt.0.0e0) go to 290
c-----------------------------------------------------------------------
c     for x.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
      if (x.lt.0.0e0) go to 290
      nz=n
      do 96 i=1,n
        y(i) = czero
   96 continue
      return
  100 continue
c-----------------------------------------------------------------------
c     recur forward for remainder of the sequence
c-----------------------------------------------------------------------
      s1 = cy(1)
      s2 = cy(2)
      c1 = csr(kflag)
      ascle = bry(kflag)
      do 120 i=ib,n
        c2 = s2
        s2 = ck*s2 + s1
        s1 = c2
        ck = ck + rz
        c2 = s2*c1
        y(i) = c2
        if (kflag.ge.3) go to 120
        c2r = real(c2)
        c2i = aimag(c2)
        c2r = abs(c2r)
        c2i = abs(c2i)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 120
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*c1
        s2 = c2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        c1 = csr(kflag)
  120 continue
  160 continue
      if (mr.eq.0) return
c-----------------------------------------------------------------------
c     analytic continuation for re(z).lt.0.0e0
c-----------------------------------------------------------------------
      nz = 0
      fmr = mr
      sgn = -sign(pi,fmr)
c-----------------------------------------------------------------------
c     cspn and csgn are coeff of k and i functions resp.
c-----------------------------------------------------------------------
      csgn = cmplx(0.0e0,sgn)
      inu = fnu
      fnf = fnu - inu
      ifn = inu + n - 1
      ang = fnf*sgn
      cpn = cos(ang)
      spn = sin(ang)
      cspn = cmplx(cpn,spn)
      if (mod(ifn,2).eq.1) cspn = -cspn
      asc = bry(1)
      kk = n
      iuf = 0
      kdflg = 1
      ib = ib-1
      ic = ib-1
      do 260 k=1,n
        fn = fnu + (kk-1)
c-----------------------------------------------------------------------
c     logic to sort out cases whose parameters were set for the k
c     function above
c-----------------------------------------------------------------------
        m=3
        if (n.gt.2) go to 175
  170   continue
        initd = init(j)
        phid = phi(j)
        zeta1d = zeta1(j)
        zeta2d = zeta2(j)
        sumd = sum(j)
        m = j
        j = 3 - j
        go to 180
  175   continue
        if ((kk.eq.n).and.(ib.lt.n)) go to 180
        if ((kk.eq.ib).or.(kk.eq.ic)) go to 170
        initd = 0
  180   continue
        call cunik(zr, fn, 1, 0, tol, initd, phid, zeta1d,
     *   zeta2d, sumd, cwrk(1,m))
        if (kode.eq.1) go to 190
        cfn = cmplx(fn,0.0e0)
        s1 = -zeta1d + cfn*(cfn/(zr+zeta2d))
        go to 200
  190   continue
        s1 = -zeta1d + zeta2d
  200   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = real(s1)
        if (abs(rs1).gt.elim) go to 250
        if (kdflg.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 210
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = abs(phid)
        rs1 = rs1 + alog(aphi)
        if (abs(rs1).gt.elim) go to 250
        if (kdflg.eq.1) iflag = 1
        if (rs1.lt.0.0e0) go to 210
        if (kdflg.eq.1) iflag = 3
  210   continue
        s2 = csgn*phid*sumd
        c2r = real(s1)
        c2i = aimag(s1)
        c2m = exp(c2r)*real(css(iflag))
        s1 = cmplx(c2m,0.0e0)*cmplx(cos(c2i),sin(c2i))
        s2 = s2*s1
        if (iflag.ne.1) go to 220
        call cuchk(s2, nw, bry(1), tol)
        if (nw.ne.0) s2 = cmplx(0.0e0,0.0e0)
  220   continue
        cy(kdflg) = s2
        c2 = s2
        s2 = s2*csr(iflag)
c-----------------------------------------------------------------------
c     add i and k functions, k sequence in y(i), i=1,n
c-----------------------------------------------------------------------
        s1 = y(kk)
        if (kode.eq.1) go to 240
        call cs1s2(zr, s1, s2, nw, asc, alim, iuf)
        nz = nz + nw
  240   continue
        y(kk) = s1*cspn + s2
        kk = kk - 1
        cspn = -cspn
        if (c2.ne.czero) go to 245
        kdflg = 1
        go to 260
  245   continue
        if (kdflg.eq.2) go to 265
        kdflg = 2
        go to 260
  250   continue
        if (rs1.gt.0.0e0) go to 290
        s2 = czero
        go to 220
  260 continue
      k = n
  265 continue
      il = n - k
      if (il.eq.0) return
c-----------------------------------------------------------------------
c     recur backward for remainder of i sequence and add in the
c     k functions, scaling the i sequence during recurrence to keep
c     intermediate arithmetic on scale near exponent extremes.
c-----------------------------------------------------------------------
      s1 = cy(1)
      s2 = cy(2)
      cs = csr(iflag)
      ascle = bry(iflag)
      fn = (inu+il)
      do 280 i=1,il
        c2 = s2
        s2 = s1 + cmplx(fn+fnf,0.0e0)*rz*s2
        s1 = c2
        fn = fn - 1.0e0
        c2 = s2*cs
        ck = c2
        c1 = y(kk)
        if (kode.eq.1) go to 270
        call cs1s2(zr, c1, c2, nw, asc, alim, iuf)
        nz = nz + nw
  270   continue
        y(kk) = c1*cspn + c2
        kk = kk - 1
        cspn = -cspn
        if (iflag.ge.3) go to 280
        c2r = real(ck)
        c2i = aimag(ck)
        c2r = abs(c2r)
        c2i = abs(c2i)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 280
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1*cs
        s2 = ck
        s1 = s1*css(iflag)
        s2 = s2*css(iflag)
        cs = csr(iflag)
  280 continue
      return
  290 continue
      nz = -1
      return
      end
