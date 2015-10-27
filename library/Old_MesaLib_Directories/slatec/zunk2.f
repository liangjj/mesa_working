*deck zunk2
      subroutine zunk2 (zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim,
     +   alim)
c***begin prologue  zunk2
c***subsidiary
c***purpose  subsidiary to zbesk
c***library   slatec
c***type      all (cunk2-a, zunk2-a)
c***author  amos, d. e., (snl)
c***description
c
c     zunk2 computes k(fnu,z) and its analytic continuation from the
c     right half plane to the left half plane by means of the
c     uniform asymptotic expansions for h(kind,fnu,zn) and j(fnu,zn)
c     where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or
c     -1. here zn=zr*i or -zr*i where zr=z if z is in the right
c     half plane or zr=-z if z is in the left half plane. mr indic-
c     ates the direction of rotation for analytic continuation.
c     nz=-1 means an overflow will occur
c
c***see also  zbesk
c***routines called  d1mach, zabs, zairy, zs1s2, zuchk, zunhj
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zunk2
c     complex ai,arg,argd,asum,asumd,bsum,bsumd,cfn,ci,cip,ck,cone,crsc,
c    *cr1,cr2,cs,cscl,csgn,cspn,csr,css,cy,czero,c1,c2,dai,phi,phid,rz,
c    *s1,s2,y,z,zb,zeta1,zeta1d,zeta2,zeta2d,zn,zr
      double precision aarg, aic, aii, air, alim, ang, aphi, argdi,
     * argdr, argi, argr, asc, ascle, asumdi, asumdr, asumi, asumr,
     * bry, bsumdi, bsumdr, bsumi, bsumr, car, cipi, cipr, cki, ckr,
     * coner, crsc, cr1i, cr1r, cr2i, cr2r, cscl, csgni, csi,
     * cspni, cspnr, csr, csrr, cssr, cyi, cyr, c1i, c1r, c2i, c2m,
     * c2r, daii, dair, elim, fmr, fn, fnf, fnu, hpi, phidi, phidr,
     * phii, phir, pi, pti, ptr, rast, razr, rs1, rzi, rzr, sar, sgn,
     * sti, str, s1i, s1r, s2i, s2r, tol, yi, yr, yy, zbi, zbr, zeroi,
     * zeror, zeta1i, zeta1r, zeta2i, zeta2r, zet1di, zet1dr, zet2di,
     * zet2dr, zi, zni, znr, zr, zri, zrr, d1mach, zabs
      integer i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, kk,
     * kode, mr, n, nai, ndai, nw, nz, idum, j, ipard, ic
      dimension bry(3), yr(n), yi(n), asumr(2), asumi(2), bsumr(2),
     * bsumi(2), phir(2), phii(2), argr(2), argi(2), zeta1r(2),
     * zeta1i(2), zeta2r(2), zeta2i(2), cyr(2), cyi(2), cipr(4),
     * cipi(4), cssr(3), csrr(3)
      external zabs
      data zeror,zeroi,coner,cr1r,cr1i,cr2r,cr2i /
     1         0.0d0, 0.0d0, 1.0d0,
     1 1.0d0,1.73205080756887729d0 , -0.5d0,-8.66025403784438647d-01 /
      data hpi, pi, aic /
     1     1.57079632679489662d+00,     3.14159265358979324d+00,
     1     1.26551212348464539d+00/
      data cipr(1),cipi(1),cipr(2),cipi(2),cipr(3),cipi(3),cipr(4),
     * cipi(4) /
     1  1.0d0,0.0d0 ,  0.0d0,-1.0d0 ,  -1.0d0,0.0d0 ,  0.0d0,1.0d0 /
c***first executable statement  zunk2
      kdflg = 1
      nz = 0
c-----------------------------------------------------------------------
c     exp(-alim)=exp(-elim)/tol=approx. one precision greater than
c     the underflow limit
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
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach(2)
      zrr = zr
      zri = zi
      if (zr.ge.0.0d0) go to 10
      zrr = -zr
      zri = -zi
   10 continue
      yy = zri
      znr = zri
      zni = -zrr
      zbr = zrr
      zbi = zri
      inu = fnu
      fnf = fnu - inu
      ang = -hpi*fnf
      car = cos(ang)
      sar = sin(ang)
      c2r = hpi*sar
      c2i = -hpi*car
      kk = mod(inu,4) + 1
      str = c2r*cipr(kk) - c2i*cipi(kk)
      sti = c2r*cipi(kk) + c2i*cipr(kk)
      csr = cr1r*str - cr1i*sti
      csi = cr1r*sti + cr1i*str
      if (yy.gt.0.0d0) go to 20
      znr = -znr
      zbi = -zbi
   20 continue
c-----------------------------------------------------------------------
c     k(fnu,z) is computed from h(2,fnu,-i*z) where z is in the first
c     quadrant. fourth quadrant values (yy.le.0.0e0) are computed by
c     conjugation since the k function is real on the positive real axis
c-----------------------------------------------------------------------
      j = 2
      do 80 i=1,n
c-----------------------------------------------------------------------
c     j flip flops between 1 and 2 in j = 3 - j
c-----------------------------------------------------------------------
        j = 3 - j
        fn = fnu + (i-1)
        call zunhj(znr, zni, fn, 0, tol, phir(j), phii(j), argr(j),
     *   argi(j), zeta1r(j), zeta1i(j), zeta2r(j), zeta2i(j), asumr(j),
     *   asumi(j), bsumr(j), bsumi(j))
        if (kode.eq.1) go to 30
        str = zbr + zeta2r(j)
        sti = zbi + zeta2i(j)
        rast = fn/zabs(str,sti)
        str = str*rast*rast
        sti = -sti*rast*rast
        s1r = zeta1r(j) - str
        s1i = zeta1i(j) - sti
        go to 40
   30   continue
        s1r = zeta1r(j) - zeta2r(j)
        s1i = zeta1i(j) - zeta2i(j)
   40   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = s1r
        if (abs(rs1).gt.elim) go to 70
        if (kdflg.eq.1) kflag = 2
        if (abs(rs1).lt.alim) go to 50
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = zabs(phir(j),phii(j))
        aarg = zabs(argr(j),argi(j))
        rs1 = rs1 + log(aphi) - 0.25d0*log(aarg) - aic
        if (abs(rs1).gt.elim) go to 70
        if (kdflg.eq.1) kflag = 1
        if (rs1.lt.0.0d0) go to 50
        if (kdflg.eq.1) kflag = 3
   50   continue
c-----------------------------------------------------------------------
c     scale s1 to keep intermediate arithmetic on scale near
c     exponent extremes
c-----------------------------------------------------------------------
        c2r = argr(j)*cr2r - argi(j)*cr2i
        c2i = argr(j)*cr2i + argi(j)*cr2r
        call zairy(c2r, c2i, 0, 2, air, aii, nai, idum)
        call zairy(c2r, c2i, 1, 2, dair, daii, ndai, idum)
        str = dair*bsumr(j) - daii*bsumi(j)
        sti = dair*bsumi(j) + daii*bsumr(j)
        ptr = str*cr2r - sti*cr2i
        pti = str*cr2i + sti*cr2r
        str = ptr + (air*asumr(j)-aii*asumi(j))
        sti = pti + (air*asumi(j)+aii*asumr(j))
        ptr = str*phir(j) - sti*phii(j)
        pti = str*phii(j) + sti*phir(j)
        s2r = ptr*csr - pti*csi
        s2i = ptr*csi + pti*csr
        str = exp(s1r)*cssr(kflag)
        s1r = str*cos(s1i)
        s1i = str*sin(s1i)
        str = s2r*s1r - s2i*s1i
        s2i = s1r*s2i + s2r*s1i
        s2r = str
        if (kflag.ne.1) go to 60
        call zuchk(s2r, s2i, nw, bry(1), tol)
        if (nw.ne.0) go to 70
   60   continue
        if (yy.le.0.0d0) s2i = -s2i
        cyr(kdflg) = s2r
        cyi(kdflg) = s2i
        yr(i) = s2r*csrr(kflag)
        yi(i) = s2i*csrr(kflag)
        str = csi
        csi = -csr
        csr = str
        if (kdflg.eq.2) go to 85
        kdflg = 2
        go to 80
   70   continue
        if (rs1.gt.0.0d0) go to 320
c-----------------------------------------------------------------------
c     for zr.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
        if (zr.lt.0.0d0) go to 320
        kdflg = 1
        yr(i)=zeror
        yi(i)=zeroi
        nz=nz+1
        str = csi
        csi =-csr
        csr = str
        if (i.eq.1) go to 80
        if ((yr(i-1).eq.zeror).and.(yi(i-1).eq.zeroi)) go to 80
        yr(i-1)=zeror
        yi(i-1)=zeroi
        nz=nz+1
   80 continue
      i = n
   85 continue
      razr = 1.0d0/zabs(zrr,zri)
      str = zrr*razr
      sti = -zri*razr
      rzr = (str+str)*razr
      rzi = (sti+sti)*razr
      ckr = fn*rzr
      cki = fn*rzi
      ib = i + 1
      if (n.lt.ib) go to 180
c-----------------------------------------------------------------------
c     test last member for underflow and overflow. set sequence to zero
c     on underflow.
c-----------------------------------------------------------------------
      fn = fnu + (n-1)
      ipard = 1
      if (mr.ne.0) ipard = 0
      call zunhj(znr, zni, fn, ipard, tol, phidr, phidi, argdr, argdi,
     * zet1dr, zet1di, zet2dr, zet2di, asumdr, asumdi, bsumdr, bsumdi)
      if (kode.eq.1) go to 90
      str = zbr + zet2dr
      sti = zbi + zet2di
      rast = fn/zabs(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = zet1dr - str
      s1i = zet1di - sti
      go to 100
   90 continue
      s1r = zet1dr - zet2dr
      s1i = zet1di - zet2di
  100 continue
      rs1 = s1r
      if (abs(rs1).gt.elim) go to 105
      if (abs(rs1).lt.alim) go to 120
c-----------------------------------------------------------------------
c     refine estimate and test
c-----------------------------------------------------------------------
      aphi = zabs(phidr,phidi)
      rs1 = rs1+log(aphi)
      if (abs(rs1).lt.elim) go to 120
  105 continue
      if (rs1.gt.0.0d0) go to 320
c-----------------------------------------------------------------------
c     for zr.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
      if (zr.lt.0.0d0) go to 320
      nz = n
      do 106 i=1,n
        yr(i) = zeror
        yi(i) = zeroi
  106 continue
      return
  120 continue
      s1r = cyr(1)
      s1i = cyi(1)
      s2r = cyr(2)
      s2i = cyi(2)
      c1r = csrr(kflag)
      ascle = bry(kflag)
      do 130 i=ib,n
        c2r = s2r
        c2i = s2i
        s2r = ckr*c2r - cki*c2i + s1r
        s2i = ckr*c2i + cki*c2r + s1i
        s1r = c2r
        s1i = c2i
        ckr = ckr + rzr
        cki = cki + rzi
        c2r = s2r*c1r
        c2i = s2i*c1r
        yr(i) = c2r
        yi(i) = c2i
        if (kflag.ge.3) go to 130
        str = abs(c2r)
        sti = abs(c2i)
        c2m = max(str,sti)
        if (c2m.le.ascle) go to 130
        kflag = kflag + 1
        ascle = bry(kflag)
        s1r = s1r*c1r
        s1i = s1i*c1r
        s2r = c2r
        s2i = c2i
        s1r = s1r*cssr(kflag)
        s1i = s1i*cssr(kflag)
        s2r = s2r*cssr(kflag)
        s2i = s2i*cssr(kflag)
        c1r = csrr(kflag)
  130 continue
  180 continue
      if (mr.eq.0) return
c-----------------------------------------------------------------------
c     analytic continuation for re(z).lt.0.0d0
c-----------------------------------------------------------------------
      nz = 0
      fmr = mr
      sgn = -dsign(pi,fmr)
c-----------------------------------------------------------------------
c     cspn and csgn are coeff of k and i functions resp.
c-----------------------------------------------------------------------
      csgni = sgn
      if (yy.le.0.0d0) csgni = -csgni
      ifn = inu + n - 1
      ang = fnf*sgn
      cspnr = cos(ang)
      cspni = sin(ang)
      if (mod(ifn,2).eq.0) go to 190
      cspnr = -cspnr
      cspni = -cspni
  190 continue
c-----------------------------------------------------------------------
c     cs=coeff of the j function to get the i function. i(fnu,z) is
c     computed from exp(i*fnu*hpi)*j(fnu,-i*z) where z is in the first
c     quadrant. fourth quadrant values (yy.le.0.0e0) are computed by
c     conjugation since the i function is real on the positive real axis
c-----------------------------------------------------------------------
      csr = sar*csgni
      csi = car*csgni
      in = mod(ifn,4) + 1
      c2r = cipr(in)
      c2i = cipi(in)
      str = csr*c2r + csi*c2i
      csi = -csr*c2i + csi*c2r
      csr = str
      asc = bry(1)
      iuf = 0
      kk = n
      kdflg = 1
      ib = ib - 1
      ic = ib - 1
      do 290 k=1,n
        fn = fnu + (kk-1)
c-----------------------------------------------------------------------
c     logic to sort out cases whose parameters were set for the k
c     function above
c-----------------------------------------------------------------------
        if (n.gt.2) go to 175
  172   continue
        phidr = phir(j)
        phidi = phii(j)
        argdr = argr(j)
        argdi = argi(j)
        zet1dr = zeta1r(j)
        zet1di = zeta1i(j)
        zet2dr = zeta2r(j)
        zet2di = zeta2i(j)
        asumdr = asumr(j)
        asumdi = asumi(j)
        bsumdr = bsumr(j)
        bsumdi = bsumi(j)
        j = 3 - j
        go to 210
  175   continue
        if ((kk.eq.n).and.(ib.lt.n)) go to 210
        if ((kk.eq.ib).or.(kk.eq.ic)) go to 172
        call zunhj(znr, zni, fn, 0, tol, phidr, phidi, argdr,
     *   argdi, zet1dr, zet1di, zet2dr, zet2di, asumdr,
     *   asumdi, bsumdr, bsumdi)
  210   continue
        if (kode.eq.1) go to 220
        str = zbr + zet2dr
        sti = zbi + zet2di
        rast = fn/zabs(str,sti)
        str = str*rast*rast
        sti = -sti*rast*rast
        s1r = -zet1dr + str
        s1i = -zet1di + sti
        go to 230
  220   continue
        s1r = -zet1dr + zet2dr
        s1i = -zet1di + zet2di
  230   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = s1r
        if (abs(rs1).gt.elim) go to 280
        if (kdflg.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 240
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = zabs(phidr,phidi)
        aarg = zabs(argdr,argdi)
        rs1 = rs1 + log(aphi) - 0.25d0*log(aarg) - aic
        if (abs(rs1).gt.elim) go to 280
        if (kdflg.eq.1) iflag = 1
        if (rs1.lt.0.0d0) go to 240
        if (kdflg.eq.1) iflag = 3
  240   continue
        call zairy(argdr, argdi, 0, 2, air, aii, nai, idum)
        call zairy(argdr, argdi, 1, 2, dair, daii, ndai, idum)
        str = dair*bsumdr - daii*bsumdi
        sti = dair*bsumdi + daii*bsumdr
        str = str + (air*asumdr-aii*asumdi)
        sti = sti + (air*asumdi+aii*asumdr)
        ptr = str*phidr - sti*phidi
        pti = str*phidi + sti*phidr
        s2r = ptr*csr - pti*csi
        s2i = ptr*csi + pti*csr
        str = exp(s1r)*cssr(iflag)
        s1r = str*cos(s1i)
        s1i = str*sin(s1i)
        str = s2r*s1r - s2i*s1i
        s2i = s2r*s1i + s2i*s1r
        s2r = str
        if (iflag.ne.1) go to 250
        call zuchk(s2r, s2i, nw, bry(1), tol)
        if (nw.eq.0) go to 250
        s2r = zeror
        s2i = zeroi
  250   continue
        if (yy.le.0.0d0) s2i = -s2i
        cyr(kdflg) = s2r
        cyi(kdflg) = s2i
        c2r = s2r
        c2i = s2i
        s2r = s2r*csrr(iflag)
        s2i = s2i*csrr(iflag)
c-----------------------------------------------------------------------
c     add i and k functions, k sequence in y(i), i=1,n
c-----------------------------------------------------------------------
        s1r = yr(kk)
        s1i = yi(kk)
        if (kode.eq.1) go to 270
        call zs1s2(zrr, zri, s1r, s1i, s2r, s2i, nw, asc, alim, iuf)
        nz = nz + nw
  270   continue
        yr(kk) = s1r*cspnr - s1i*cspni + s2r
        yi(kk) = s1r*cspni + s1i*cspnr + s2i
        kk = kk - 1
        cspnr = -cspnr
        cspni = -cspni
        str = csi
        csi = -csr
        csr = str
        if (c2r.ne.0.0d0 .or. c2i.ne.0.0d0) go to 255
        kdflg = 1
        go to 290
  255   continue
        if (kdflg.eq.2) go to 295
        kdflg = 2
        go to 290
  280   continue
        if (rs1.gt.0.0d0) go to 320
        s2r = zeror
        s2i = zeroi
        go to 250
  290 continue
      k = n
  295 continue
      il = n - k
      if (il.eq.0) return
c-----------------------------------------------------------------------
c     recur backward for remainder of i sequence and add in the
c     k functions, scaling the i sequence during recurrence to keep
c     intermediate arithmetic on scale near exponent extremes.
c-----------------------------------------------------------------------
      s1r = cyr(1)
      s1i = cyi(1)
      s2r = cyr(2)
      s2i = cyi(2)
      csr = csrr(iflag)
      ascle = bry(iflag)
      fn = inu+il
      do 310 i=1,il
        c2r = s2r
        c2i = s2i
        s2r = s1r + (fn+fnf)*(rzr*c2r-rzi*c2i)
        s2i = s1i + (fn+fnf)*(rzr*c2i+rzi*c2r)
        s1r = c2r
        s1i = c2i
        fn = fn - 1.0d0
        c2r = s2r*csr
        c2i = s2i*csr
        ckr = c2r
        cki = c2i
        c1r = yr(kk)
        c1i = yi(kk)
        if (kode.eq.1) go to 300
        call zs1s2(zrr, zri, c1r, c1i, c2r, c2i, nw, asc, alim, iuf)
        nz = nz + nw
  300   continue
        yr(kk) = c1r*cspnr - c1i*cspni + c2r
        yi(kk) = c1r*cspni + c1i*cspnr + c2i
        kk = kk - 1
        cspnr = -cspnr
        cspni = -cspni
        if (iflag.ge.3) go to 310
        c2r = abs(ckr)
        c2i = abs(cki)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 310
        iflag = iflag + 1
        ascle = bry(iflag)
        s1r = s1r*csr
        s1i = s1i*csr
        s2r = ckr
        s2i = cki
        s1r = s1r*cssr(iflag)
        s1i = s1i*cssr(iflag)
        s2r = s2r*cssr(iflag)
        s2i = s2i*cssr(iflag)
        csr = csrr(iflag)
  310 continue
      return
  320 continue
      nz = -1
      return
      end
