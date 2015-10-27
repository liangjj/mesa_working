*deck zuni2
      subroutine zuni2 (zr, zi, fnu, kode, n, yr, yi, nz, nlast, fnul,
     +   tol, elim, alim)
c***begin prologue  zuni2
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cuni2-a, zuni2-a)
c***author  amos, d. e., (snl)
c***description
c
c     zuni2 computes i(fnu,z) in the right half plane by means of
c     uniform asymptotic expansion for j(fnu,zn) where zn is z*i
c     or -z*i and zn is in the right half plane also.
c
c     fnul is the smallest order permitted for the asymptotic
c     expansion. nlast=0 means all of the y values were set.
c     nlast.ne.0 is the number left to be computed by another
c     formula for orders fnu to fnu+nlast-1 because fnu+nlast-1.lt.fnul.
c     y(i)=czero for i=nlast+1,n
c
c***see also  zbesi, zbesk
c***routines called  d1mach, zabs, zairy, zuchk, zunhj, zuoik
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zuni2
c     complex ai,arg,asum,bsum,cfn,ci,cid,cip,cone,crsc,cscl,csr,css,
c    *czero,c1,c2,dai,phi,rz,s1,s2,y,z,zb,zeta1,zeta2,zn
      double precision aarg, aic, aii, air, alim, ang, aphi, argi,
     * argr, ascle, asumi, asumr, bry, bsumi, bsumr, cidi, cipi, cipr,
     * coner, crsc, cscl, csrr, cssr, c1r, c2i, c2m, c2r, daii,
     * dair, elim, fn, fnu, fnul, hpi, phii, phir, rast, raz, rs1, rzi,
     * rzr, sti, str, s1i, s1r, s2i, s2r, tol, yi, yr, zbi, zbr, zeroi,
     * zeror, zeta1i, zeta1r, zeta2i, zeta2r, zi, zni, znr, zr, cyr,
     * cyi, d1mach, zabs, car, sar
      integer i, iflag, in, inu, j, k, kode, n, nai, nd, ndai, nlast,
     * nn, nuf, nw, nz, idum
      dimension bry(3), yr(n), yi(n), cipr(4), cipi(4), cssr(3),
     * csrr(3), cyr(2), cyi(2)
      external zabs
      data zeror,zeroi,coner / 0.0d0, 0.0d0, 1.0d0 /
      data cipr(1),cipi(1),cipr(2),cipi(2),cipr(3),cipi(3),cipr(4),
     * cipi(4)/ 1.0d0,0.0d0, 0.0d0,1.0d0, -1.0d0,0.0d0, 0.0d0,-1.0d0/
      data hpi, aic  /
     1      1.57079632679489662d+00,     1.265512123484645396d+00/
c***first executable statement  zuni2
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
c     zn is in the right half plane after rotation by ci or -ci
c-----------------------------------------------------------------------
      znr = zi
      zni = -zr
      zbr = zr
      zbi = zi
      cidi = -coner
      inu = fnu
      ang = hpi*(fnu-inu)
      c2r = cos(ang)
      c2i = sin(ang)
      car = c2r
      sar = c2i
      in = inu + n - 1
      in = mod(in,4) + 1
      str = c2r*cipr(in) - c2i*cipi(in)
      c2i = c2r*cipi(in) + c2i*cipr(in)
      c2r = str
      if (zi.gt.0.0d0) go to 10
      znr = -znr
      zbi = -zbi
      cidi = -cidi
      c2i = -c2i
   10 continue
c-----------------------------------------------------------------------
c     check for underflow and overflow on first member
c-----------------------------------------------------------------------
      fn = max(fnu,1.0d0)
      call zunhj(znr, zni, fn, 1, tol, phir, phii, argr, argi, zeta1r,
     * zeta1i, zeta2r, zeta2i, asumr, asumi, bsumr, bsumi)
      if (kode.eq.1) go to 20
      str = zbr + zeta2r
      sti = zbi + zeta2i
      rast = fn/zabs(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = -zeta1r + str
      s1i = -zeta1i + sti
      go to 30
   20 continue
      s1r = -zeta1r + zeta2r
      s1i = -zeta1i + zeta2i
   30 continue
      rs1 = s1r
      if (abs(rs1).gt.elim) go to 150
   40 continue
      nn = min(2,nd)
      do 90 i=1,nn
        fn = fnu + (nd-i)
        call zunhj(znr, zni, fn, 0, tol, phir, phii, argr, argi,
     *   zeta1r, zeta1i, zeta2r, zeta2i, asumr, asumi, bsumr, bsumi)
        if (kode.eq.1) go to 50
        str = zbr + zeta2r
        sti = zbi + zeta2i
        rast = fn/zabs(str,sti)
        str = str*rast*rast
        sti = -sti*rast*rast
        s1r = -zeta1r + str
        s1i = -zeta1i + sti + abs(zi)
        go to 60
   50   continue
        s1r = -zeta1r + zeta2r
        s1i = -zeta1i + zeta2i
   60   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = s1r
        if (abs(rs1).gt.elim) go to 120
        if (i.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 70
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
        aphi = zabs(phir,phii)
        aarg = zabs(argr,argi)
        rs1 = rs1 + log(aphi) - 0.25d0*log(aarg) - aic
        if (abs(rs1).gt.elim) go to 120
        if (i.eq.1) iflag = 1
        if (rs1.lt.0.0d0) go to 70
        if (i.eq.1) iflag = 3
   70   continue
c-----------------------------------------------------------------------
c     scale s1 to keep intermediate arithmetic on scale near
c     exponent extremes
c-----------------------------------------------------------------------
        call zairy(argr, argi, 0, 2, air, aii, nai, idum)
        call zairy(argr, argi, 1, 2, dair, daii, ndai, idum)
        str = dair*bsumr - daii*bsumi
        sti = dair*bsumi + daii*bsumr
        str = str + (air*asumr-aii*asumi)
        sti = sti + (air*asumi+aii*asumr)
        s2r = phir*str - phii*sti
        s2i = phir*sti + phii*str
        str = exp(s1r)*cssr(iflag)
        s1r = str*cos(s1i)
        s1i = str*sin(s1i)
        str = s2r*s1r - s2i*s1i
        s2i = s2r*s1i + s2i*s1r
        s2r = str
        if (iflag.ne.1) go to 80
        call zuchk(s2r, s2i, nw, bry(1), tol)
        if (nw.ne.0) go to 120
   80   continue
        if (zi.le.0.0d0) s2i = -s2i
        str = s2r*c2r - s2i*c2i
        s2i = s2r*c2i + s2i*c2r
        s2r = str
        cyr(i) = s2r
        cyi(i) = s2i
        j = nd - i + 1
        yr(j) = s2r*csrr(iflag)
        yi(j) = s2i*csrr(iflag)
        str = -c2i*cidi
        c2i = c2r*cidi
        c2r = str
   90 continue
      if (nd.le.2) go to 110
      raz = 1.0d0/zabs(zr,zi)
      str = zr*raz
      sti = -zi*raz
      rzr = (str+str)*raz
      rzi = (sti+sti)*raz
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
      do 100 i=3,nd
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
        if (iflag.ge.3) go to 100
        str = abs(c2r)
        sti = abs(c2i)
        c2m = max(str,sti)
        if (c2m.le.ascle) go to 100
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
  100 continue
  110 continue
      return
  120 continue
      if (rs1.gt.0.0d0) go to 140
c-----------------------------------------------------------------------
c     set underflow and update parameters
c-----------------------------------------------------------------------
      yr(nd) = zeror
      yi(nd) = zeroi
      nz = nz + 1
      nd = nd - 1
      if (nd.eq.0) go to 110
      call zuoik(zr, zi, fnu, kode, 1, nd, yr, yi, nuf, tol, elim, alim)
      if (nuf.lt.0) go to 140
      nd = nd - nuf
      nz = nz + nuf
      if (nd.eq.0) go to 110
      fn = fnu + (nd-1)
      if (fn.lt.fnul) go to 130
c      fn = cidi
c      j = nuf + 1
c      k = mod(j,4) + 1
c      s1r = cipr(k)
c      s1i = cipi(k)
c      if (fn.lt.0.0d0) s1i = -s1i
c      str = c2r*s1r - c2i*s1i
c      c2i = c2r*s1i + c2i*s1r
c      c2r = str
      in = inu + nd - 1
      in = mod(in,4) + 1
      c2r = car*cipr(in) - sar*cipi(in)
      c2i = car*cipi(in) + sar*cipr(in)
      if (zi.le.0.0d0) c2i = -c2i
      go to 40
  130 continue
      nlast = nd
      return
  140 continue
      nz = -1
      return
  150 continue
      if (rs1.gt.0.0d0) go to 140
      nz = n
      do 160 i=1,n
        yr(i) = zeror
        yi(i) = zeroi
  160 continue
      return
      end
