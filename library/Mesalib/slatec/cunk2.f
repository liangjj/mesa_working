*deck cunk2
      subroutine cunk2 (z, fnu, kode, mr, n, y, nz, tol, elim, alim)
c***begin prologue  cunk2
c***subsidiary
c***purpose  subsidiary to cbesk
c***library   slatec
c***type      all (cunk2-a, zunk2-a)
c***author  amos, d. e., (snl)
c***description
c
c     cunk2 computes k(fnu,z) and its analytic continuation from the
c     right half plane to the left half plane by means of the
c     uniform asymptotic expansions for h(kind,fnu,zn) and j(fnu,zn)
c     where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or
c     -1. here zn=zr*i or -zr*i where zr=z if z is in the right
c     half plane or zr=-z if z is in the left half plane. mr indic-
c     ates the direction of rotation for analytic continuation.
c     nz=-1 means an overflow will occur
c
c***see also  cbesk
c***routines called  cairy, cs1s2, cuchk, cunhj, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cunk2
      complex ai, arg, asum, bsum, cfn, ci, cip,
     * ck, cone, crsc, cr1, cr2, cs, cscl, csgn, cspn, csr, css, cy,
     * czero, c1, c2, dai, phi,  rz, s1, s2, y, z, zb, zeta1,
     * zeta2, zn, zr, phid, argd, zeta1d, zeta2d, asumd, bsumd
      real aarg, aic, alim, ang, aphi, asc, ascle, bry, car, cpn, c2i,
     * c2m, c2r, elim, fmr, fn, fnf, fnu, hpi, pi, rs1, sar, sgn, spn,
     * tol, x, yy, r1mach
      integer i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, kk,
     * kode, mr, n, nai, ndai, nw, nz, idum, j, ipard, ic
      dimension bry(3), y(n), asum(2), bsum(2), phi(2), arg(2),
     * zeta1(2), zeta2(2), cy(2), cip(4), css(3), csr(3)
      data czero, cone, ci, cr1, cr2 /
     1         (0.0e0,0.0e0),(1.0e0,0.0e0),(0.0e0,1.0e0),
     1(1.0e0,1.73205080756887729e0),(-0.5e0,-8.66025403784438647e-01)/
      data hpi, pi, aic /
     1     1.57079632679489662e+00,     3.14159265358979324e+00,
     1     1.26551212348464539e+00/
      data cip(1),cip(2),cip(3),cip(4)/
     1 (1.0e0,0.0e0), (0.0e0,-1.0e0), (-1.0e0,0.0e0), (0.0e0,1.0e0)/
c***first executable statement  cunk2
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
      yy = aimag(zr)
      zn = -zr*ci
      zb = zr
      inu = fnu
      fnf = fnu - inu
      ang = -hpi*fnf
      car = cos(ang)
      sar = sin(ang)
      cpn = -hpi*car
      spn = -hpi*sar
      c2 = cmplx(-spn,cpn)
      kk = mod(inu,4) + 1
      cs = cr1*c2*cip(kk)
      if (yy.gt.0.0e0) go to 10
      zn = conjg(-zn)
      zb = conjg(zb)
   10 continue
c-----------------------------------------------------------------------
c     k(fnu,z) is computed from h(2,fnu,-i*z) where z is in the first
c     quadrant. fourth quadrant values (yy.le.0.0e0) are computed by
c     conjugation since the k function is real on the positive real axis
c-----------------------------------------------------------------------
      j = 2
      do 70 i=1,n
c-----------------------------------------------------------------------
c     j flip flops between 1 and 2 in j = 3 - j
c-----------------------------------------------------------------------
        j = 3 - j
        fn = fnu + (i-1)
        call cunhj(zn, fn, 0, tol, phi(j), arg(j), zeta1(j), zeta2(j),
     *   asum(j), bsum(j))
        if (kode.eq.1) go to 20
        cfn = cmplx(fn,0.0e0)
        s1 = zeta1(j) - cfn*(cfn/(zb+zeta2(j)))
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
        aarg = abs(arg(j))
        rs1 = rs1 + alog(aphi) - 0.25e0*alog(aarg) - aic
        if (abs(rs1).gt.elim) go to 60
        if (kdflg.eq.1) kflag = 1
        if (rs1.lt.0.0e0) go to 40
        if (kdflg.eq.1) kflag = 3
   40   continue
c-----------------------------------------------------------------------
c     scale s1 to keep intermediate arithmetic on scale near
c     exponent extremes
c-----------------------------------------------------------------------
        c2 = arg(j)*cr2
        call cairy(c2, 0, 2, ai, nai, idum)
        call cairy(c2, 1, 2, dai, ndai, idum)
        s2 = cs*phi(j)*(ai*asum(j)+cr2*dai*bsum(j))
        c2r = real(s1)
        c2i = aimag(s1)
        c2m = exp(c2r)*real(css(kflag))
        s1 = cmplx(c2m,0.0e0)*cmplx(cos(c2i),sin(c2i))
        s2 = s2*s1
        if (kflag.ne.1) go to 50
        call cuchk(s2, nw, bry(1), tol)
        if (nw.ne.0) go to 60
   50   continue
        if (yy.le.0.0e0) s2 = conjg(s2)
        cy(kdflg) = s2
        y(i) = s2*csr(kflag)
        cs = -ci*cs
        if (kdflg.eq.2) go to 75
        kdflg = 2
        go to 70
   60   continue
        if (rs1.gt.0.0e0) go to 300
c-----------------------------------------------------------------------
c     for x.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
        if (x.lt.0.0e0) go to 300
        kdflg = 1
        y(i) = czero
        cs = -ci*cs
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
      ib = i + 1
      if (n.lt.ib) go to 170
c-----------------------------------------------------------------------
c     test last member for underflow and overflow, set sequence to zero
c     on underflow
c-----------------------------------------------------------------------
      fn = fnu+(n-1)
      ipard = 1
      if (mr.ne.0) ipard = 0
      call cunhj(zn,fn,ipard,tol,phid,argd,zeta1d,zeta2d,asumd,bsumd)
      if (kode.eq.1) go to 80
      cfn=cmplx(fn,0.0e0)
      s1=zeta1d-cfn*(cfn/(zb+zeta2d))
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
      aarg = abs(argd)
      rs1=rs1+alog(aphi)-0.25e0*alog(aarg)-aic
      if (abs(rs1).lt.elim) go to 100
   95 continue
      if (rs1.gt.0.0e0) go to 300
c-----------------------------------------------------------------------
c     for x.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
      if (x.lt.0.0e0) go to 300
      nz=n
      do 96 i=1,n
        y(i) = czero
   96 continue
      return
  100 continue
c-----------------------------------------------------------------------
c     scaled forward recurrence for remainder of the sequence
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
  170 continue
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
      if (yy.le.0.0e0) csgn = conjg(csgn)
      ifn = inu + n - 1
      ang = fnf*sgn
      cpn = cos(ang)
      spn = sin(ang)
      cspn = cmplx(cpn,spn)
      if (mod(ifn,2).eq.1) cspn = -cspn
c-----------------------------------------------------------------------
c     cs=coeff of the j function to get the i function. i(fnu,z) is
c     computed from exp(i*fnu*hpi)*j(fnu,-i*z) where z is in the first
c     quadrant. fourth quadrant values (yy.le.0.0e0) are computed by
c     conjugation since the i function is real on the positive real axis
c-----------------------------------------------------------------------
      cs = cmplx(car,-sar)*csgn
      in = mod(ifn,4) + 1
      c2 = cip(in)
      cs = cs*conjg(c2)
      asc = bry(1)
      kk = n
      kdflg = 1
      ib = ib-1
      ic = ib-1
      iuf = 0
      do 270 k=1,n
c-----------------------------------------------------------------------
c     logic to sort out cases whose parameters were set for the k
c     function above
c-----------------------------------------------------------------------
        fn = fnu+(kk-1)
        if (n.gt.2) go to 180
  175   continue
        phid = phi(j)
        argd = arg(j)
        zeta1d = zeta1(j)
        zeta2d = zeta2(j)
        asumd = asum(j)
        bsumd = bsum(j)
        j = 3 - j
        go to 190
  180   continue
        if ((kk.eq.n).and.(ib.lt.n)) go to 190
        if ((kk.eq.ib).or.(kk.eq.ic)) go to 175
        call cunhj(zn, fn, 0, tol, phid, argd, zeta1d, zeta2d,
     *   asumd, bsumd)
  190   continue
        if (kode.eq.1) go to 200
        cfn = cmplx(fn,0.0e0)
        s1 = -zeta1d + cfn*(cfn/(zb+zeta2d))
        go to 210
  200   continue
        s1 = -zeta1d + zeta2d
  210   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = real(s1)
        if (abs(rs1).gt.elim) go to 260
        if (kdflg.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 220
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = abs(phid)
        aarg = abs(argd)
        rs1 = rs1 + alog(aphi) - 0.25e0*alog(aarg) - aic
        if (abs(rs1).gt.elim) go to 260
        if (kdflg.eq.1) iflag = 1
        if (rs1.lt.0.0e0) go to 220
        if (kdflg.eq.1) iflag = 3
  220   continue
        call cairy(argd, 0, 2, ai, nai, idum)
        call cairy(argd, 1, 2, dai, ndai, idum)
        s2 = cs*phid*(ai*asumd+dai*bsumd)
        c2r = real(s1)
        c2i = aimag(s1)
        c2m = exp(c2r)*real(css(iflag))
        s1 = cmplx(c2m,0.0e0)*cmplx(cos(c2i),sin(c2i))
        s2 = s2*s1
        if (iflag.ne.1) go to 230
        call cuchk(s2, nw, bry(1), tol)
        if (nw.ne.0) s2 = cmplx(0.0e0,0.0e0)
  230   continue
        if (yy.le.0.0e0) s2 = conjg(s2)
        cy(kdflg) = s2
        c2 = s2
        s2 = s2*csr(iflag)
c-----------------------------------------------------------------------
c     add i and k functions, k sequence in y(i), i=1,n
c-----------------------------------------------------------------------
        s1 = y(kk)
        if (kode.eq.1) go to 250
        call cs1s2(zr, s1, s2, nw, asc, alim, iuf)
        nz = nz + nw
  250   continue
        y(kk) = s1*cspn + s2
        kk = kk - 1
        cspn = -cspn
        cs = -cs*ci
        if (c2.ne.czero) go to 255
        kdflg = 1
        go to 270
  255   continue
        if (kdflg.eq.2) go to 275
        kdflg = 2
        go to 270
  260   continue
        if (rs1.gt.0.0e0) go to 300
        s2 = czero
        go to 230
  270 continue
      k = n
  275 continue
      il = n-k
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
      fn = inu+il
      do 290 i=1,il
        c2 = s2
        s2 = s1 + cmplx(fn+fnf,0.0e0)*rz*s2
        s1 = c2
        fn = fn - 1.0e0
        c2 = s2*cs
        ck = c2
        c1 = y(kk)
        if (kode.eq.1) go to 280
        call cs1s2(zr, c1, c2, nw, asc, alim, iuf)
        nz = nz + nw
  280   continue
        y(kk) = c1*cspn + c2
        kk = kk - 1
        cspn = -cspn
        if (iflag.ge.3) go to 290
        c2r = real(ck)
        c2i = aimag(ck)
        c2r = abs(c2r)
        c2i = abs(c2i)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 290
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1*cs
        s2 = ck
        s1 = s1*css(iflag)
        s2 = s2*css(iflag)
        cs = csr(iflag)
  290 continue
      return
  300 continue
      nz = -1
      return
      end
