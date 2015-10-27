*deck zuni1
      subroutine zuni1 (zr, zi, fnu, kode, n, yr, yi, nz, nlast, fnul,
     +   tol, elim, alim)
c***begin prologue  zuni1
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cuni1-a, zuni1-a)
c***author  amos, d. e., (snl)
c***description
c
c     zuni1 computes i(fnu,z)  by means of the uniform asymptotic
c     expansion for i(fnu,z) in -pi/3.le.arg z.le.pi/3.
c
c     fnul is the smallest order permitted for the asymptotic
c     expansion. nlast=0 means all of the y values were set.
c     nlast.ne.0 is the number left to be computed by another
c     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1.lt.fnul.
c     y(i)=czero for i=nlast+1,n
c
c***see also  zbesi, zbesk
c***routines called  d1mach, zabs, zuchk, zunik, zuoik
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zuni1
c     complex cfn,cone,crsc,cscl,csr,css,cwrk,czero,c1,c2,phi,rz,sum,s1,
c    *s2,y,z,zeta1,zeta2
      double precision alim, aphi, ascle, bry, coner, crsc,
     * cscl, csrr, cssr, cwrki, cwrkr, c1r, c2i, c2m, c2r, elim, fn,
     * fnu, fnul, phii, phir, rast, rs1, rzi, rzr, sti, str, sumi,
     * sumr, s1i, s1r, s2i, s2r, tol, yi, yr, zeroi, zeror, zeta1i,
     * zeta1r, zeta2i, zeta2r, zi, zr, cyr, cyi, d1mach, zabs
      integer i, iflag, init, k, kode, m, n, nd, nlast, nn, nuf, nw, nz
      dimension bry(3), yr(n), yi(n), cwrkr(16), cwrki(16), cssr(3),
     * csrr(3), cyr(2), cyi(2)
      external zabs
      data zeror,zeroi,coner / 0.0d0, 0.0d0, 1.0d0 /
c***first executable statement  zuni1
      nz = 0
      nd = n
      nlast = 0
c-----------------------------------------------------------------------
c     computed values with exponents between alim and elim in mag-
c     nitude are scaled to keep intermediate arithmetic on scale,
c     exp(alim)=exp(elim)*tol
c-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      crsc = tol
      cssr(1) = cscl
      cssr(2) = coner
      cssr(3) = crsc
      csrr(1) = crsc
      csrr(2) = coner
      csrr(3) = cscl
      bry(1) = 1.0d+3*d1mach(1)/tol
c-----------------------------------------------------------------------
c     check for underflow and overflow on first member
c-----------------------------------------------------------------------
      fn = max(fnu,1.0d0)
      init = 0
      call zunik(zr, zi, fn, 1, 1, tol, init, phir, phii, zeta1r,
     * zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki)
      if (kode.eq.1) go to 10
      str = zr + zeta2r
      sti = zi + zeta2i
      rast = fn/zabs(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = -zeta1r + str
      s1i = -zeta1i + sti
      go to 20
   10 continue
      s1r = -zeta1r + zeta2r
      s1i = -zeta1i + zeta2i
   20 continue
      rs1 = s1r
      if (abs(rs1).gt.elim) go to 130
   30 continue
      nn = min(2,nd)
      do 80 i=1,nn
        fn = fnu + (nd-i)
        init = 0
        call zunik(zr, zi, fn, 1, 0, tol, init, phir, phii, zeta1r,
     *   zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki)
        if (kode.eq.1) go to 40
        str = zr + zeta2r
        sti = zi + zeta2i
        rast = fn/zabs(str,sti)
        str = str*rast*rast
        sti = -sti*rast*rast
        s1r = -zeta1r + str
        s1i = -zeta1i + sti + zi
        go to 50
   40   continue
        s1r = -zeta1r + zeta2r
        s1i = -zeta1i + zeta2i
   50   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = s1r
        if (abs(rs1).gt.elim) go to 110
        if (i.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 60
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = zabs(phir,phii)
        rs1 = rs1 + log(aphi)
        if (abs(rs1).gt.elim) go to 110
        if (i.eq.1) iflag = 1
        if (rs1.lt.0.0d0) go to 60
        if (i.eq.1) iflag = 3
   60   continue
c-----------------------------------------------------------------------
c     scale s1 if abs(s1).lt.ascle
c-----------------------------------------------------------------------
        s2r = phir*sumr - phii*sumi
        s2i = phir*sumi + phii*sumr
        str = exp(s1r)*cssr(iflag)
        s1r = str*cos(s1i)
        s1i = str*sin(s1i)
        str = s2r*s1r - s2i*s1i
        s2i = s2r*s1i + s2i*s1r
        s2r = str
        if (iflag.ne.1) go to 70
        call zuchk(s2r, s2i, nw, bry(1), tol)
        if (nw.ne.0) go to 110
   70   continue
        cyr(i) = s2r
        cyi(i) = s2i
        m = nd - i + 1
        yr(m) = s2r*csrr(iflag)
        yi(m) = s2i*csrr(iflag)
   80 continue
      if (nd.le.2) go to 100
      rast = 1.0d0/zabs(zr,zi)
      str = zr*rast
      sti = -zi*rast
      rzr = (str+str)*rast
      rzi = (sti+sti)*rast
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach(2)
      s1r = cyr(1)
      s1i = cyi(1)
      s2r = cyr(2)
      s2i = cyi(2)
      c1r = csrr(iflag)
      ascle = bry(iflag)
      k = nd - 2
      fn = k
      do 90 i=3,nd
        c2r = s2r
        c2i = s2i
        s2r = s1r + (fnu+fn)*(rzr*c2r-rzi*c2i)
        s2i = s1i + (fnu+fn)*(rzr*c2i+rzi*c2r)
        s1r = c2r
        s1i = c2i
        c2r = s2r*c1r
        c2i = s2i*c1r
        yr(k) = c2r
        yi(k) = c2i
        k = k - 1
        fn = fn - 1.0d0
        if (iflag.ge.3) go to 90
        str = abs(c2r)
        sti = abs(c2i)
        c2m = max(str,sti)
        if (c2m.le.ascle) go to 90
        iflag = iflag + 1
        ascle = bry(iflag)
        s1r = s1r*c1r
        s1i = s1i*c1r
        s2r = c2r
        s2i = c2i
        s1r = s1r*cssr(iflag)
        s1i = s1i*cssr(iflag)
        s2r = s2r*cssr(iflag)
        s2i = s2i*cssr(iflag)
        c1r = csrr(iflag)
   90 continue
  100 continue
      return
c-----------------------------------------------------------------------
c     set underflow and update parameters
c-----------------------------------------------------------------------
  110 continue
      if (rs1.gt.0.0d0) go to 120
      yr(nd) = zeror
      yi(nd) = zeroi
      nz = nz + 1
      nd = nd - 1
      if (nd.eq.0) go to 100
      call zuoik(zr, zi, fnu, kode, 1, nd, yr, yi, nuf, tol, elim, alim)
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
      if (rs1.gt.0.0d0) go to 120
      nz = n
      do 140 i=1,n
        yr(i) = zeror
        yi(i) = zeroi
  140 continue
      return
      end
