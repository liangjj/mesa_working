*deck zunk1
      subroutine zunk1 (zr, zi, fnu, kode, mr, n, yr, yi, nz, tol, elim,
     +   alim)
c***begin prologue  zunk1
c***subsidiary
c***purpose  subsidiary to zbesk
c***library   slatec
c***type      all (cunk1-a, zunk1-a)
c***author  amos, d. e., (snl)
c***description
c
c     zunk1 computes k(fnu,z) and its analytic continuation from the
c     right half plane to the left half plane by means of the
c     uniform asymptotic expansion.
c     mr indicates the direction of rotation for analytic continuation.
c     nz=-1 means an overflow will occur
c
c***see also  zbesk
c***routines called  d1mach, zabs, zs1s2, zuchk, zunik
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zunk1
c     complex cfn,ck,cone,crsc,cs,cscl,csgn,cspn,csr,css,cwrk,cy,czero,
c    *c1,c2,phi,phid,rz,sum,sumd,s1,s2,y,z,zeta1,zeta1d,zeta2,zeta2d,zr
      double precision alim, ang, aphi, asc, ascle, bry, cki, ckr,
     * coner, crsc, cscl, csgni, cspni, cspnr, csr, csrr, cssr,
     * cwrki, cwrkr, cyi, cyr, c1i, c1r, c2i, c2m, c2r, elim, fmr, fn,
     * fnf, fnu, phidi, phidr, phii, phir, pi, rast, razr, rs1, rzi,
     * rzr, sgn, sti, str, sumdi, sumdr, sumi, sumr, s1i, s1r, s2i,
     * s2r, tol, yi, yr, zeroi, zeror, zeta1i, zeta1r, zeta2i, zeta2r,
     * zet1di, zet1dr, zet2di, zet2dr, zi, zr, zri, zrr, d1mach, zabs
      integer i, ib, iflag, ifn, il, init, inu, iuf, k, kdflg, kflag,
     * kk, kode, mr, n, nw, nz, initd, ic, ipard, j, m
      dimension bry(3), init(2), yr(n), yi(n), sumr(2), sumi(2),
     * zeta1r(2), zeta1i(2), zeta2r(2), zeta2i(2), cyr(2), cyi(2),
     * cwrkr(16,3), cwrki(16,3), cssr(3), csrr(3), phir(2), phii(2)
      external zabs
      data zeror,zeroi,coner / 0.0d0, 0.0d0, 1.0d0 /
      data pi / 3.14159265358979324d0 /
c***first executable statement  zunk1
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
      j = 2
      do 70 i=1,n
c-----------------------------------------------------------------------
c     j flip flops between 1 and 2 in j = 3 - j
c-----------------------------------------------------------------------
        j = 3 - j
        fn = fnu + (i-1)
        init(j) = 0
        call zunik(zrr, zri, fn, 2, 0, tol, init(j), phir(j), phii(j),
     *   zeta1r(j), zeta1i(j), zeta2r(j), zeta2i(j), sumr(j), sumi(j),
     *   cwrkr(1,j), cwrki(1,j))
        if (kode.eq.1) go to 20
        str = zrr + zeta2r(j)
        sti = zri + zeta2i(j)
        rast = fn/zabs(str,sti)
        str = str*rast*rast
        sti = -sti*rast*rast
        s1r = zeta1r(j) - str
        s1i = zeta1i(j) - sti
        go to 30
   20   continue
        s1r = zeta1r(j) - zeta2r(j)
        s1i = zeta1i(j) - zeta2i(j)
   30   continue
        rs1 = s1r
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        if (abs(rs1).gt.elim) go to 60
        if (kdflg.eq.1) kflag = 2
        if (abs(rs1).lt.alim) go to 40
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = zabs(phir(j),phii(j))
        rs1 = rs1 + log(aphi)
        if (abs(rs1).gt.elim) go to 60
        if (kdflg.eq.1) kflag = 1
        if (rs1.lt.0.0d0) go to 40
        if (kdflg.eq.1) kflag = 3
   40   continue
c-----------------------------------------------------------------------
c     scale s1 to keep intermediate arithmetic on scale near
c     exponent extremes
c-----------------------------------------------------------------------
        s2r = phir(j)*sumr(j) - phii(j)*sumi(j)
        s2i = phir(j)*sumi(j) + phii(j)*sumr(j)
        str = exp(s1r)*cssr(kflag)
        s1r = str*cos(s1i)
        s1i = str*sin(s1i)
        str = s2r*s1r - s2i*s1i
        s2i = s1r*s2i + s2r*s1i
        s2r = str
        if (kflag.ne.1) go to 50
        call zuchk(s2r, s2i, nw, bry(1), tol)
        if (nw.ne.0) go to 60
   50   continue
        cyr(kdflg) = s2r
        cyi(kdflg) = s2i
        yr(i) = s2r*csrr(kflag)
        yi(i) = s2i*csrr(kflag)
        if (kdflg.eq.2) go to 75
        kdflg = 2
        go to 70
   60   continue
        if (rs1.gt.0.0d0) go to 300
c-----------------------------------------------------------------------
c     for zr.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
        if (zr.lt.0.0d0) go to 300
        kdflg = 1
        yr(i)=zeror
        yi(i)=zeroi
        nz=nz+1
        if (i.eq.1) go to 70
        if ((yr(i-1).eq.zeror).and.(yi(i-1).eq.zeroi)) go to 70
        yr(i-1)=zeror
        yi(i-1)=zeroi
        nz=nz+1
   70 continue
      i = n
   75 continue
      razr = 1.0d0/zabs(zrr,zri)
      str = zrr*razr
      sti = -zri*razr
      rzr = (str+str)*razr
      rzi = (sti+sti)*razr
      ckr = fn*rzr
      cki = fn*rzi
      ib = i + 1
      if (n.lt.ib) go to 160
c-----------------------------------------------------------------------
c     test last member for underflow and overflow. set sequence to zero
c     on underflow.
c-----------------------------------------------------------------------
      fn = fnu + (n-1)
      ipard = 1
      if (mr.ne.0) ipard = 0
      initd = 0
      call zunik(zrr, zri, fn, 2, ipard, tol, initd, phidr, phidi,
     * zet1dr, zet1di, zet2dr, zet2di, sumdr, sumdi, cwrkr(1,3),
     * cwrki(1,3))
      if (kode.eq.1) go to 80
      str = zrr + zet2dr
      sti = zri + zet2di
      rast = fn/zabs(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = zet1dr - str
      s1i = zet1di - sti
      go to 90
   80 continue
      s1r = zet1dr - zet2dr
      s1i = zet1di - zet2di
   90 continue
      rs1 = s1r
      if (abs(rs1).gt.elim) go to 95
      if (abs(rs1).lt.alim) go to 100
c-----------------------------------------------------------------------
c     refine estimate and test
c-----------------------------------------------------------------------
      aphi = zabs(phidr,phidi)
      rs1 = rs1+log(aphi)
      if (abs(rs1).lt.elim) go to 100
   95 continue
      if (abs(rs1).gt.0.0d0) go to 300
c-----------------------------------------------------------------------
c     for zr.lt.0.0, the i function to be added will overflow
c-----------------------------------------------------------------------
      if (zr.lt.0.0d0) go to 300
      nz = n
      do 96 i=1,n
        yr(i) = zeror
        yi(i) = zeroi
   96 continue
      return
c-----------------------------------------------------------------------
c     forward recur for remainder of the sequence
c-----------------------------------------------------------------------
  100 continue
      s1r = cyr(1)
      s1i = cyi(1)
      s2r = cyr(2)
      s2i = cyi(2)
      c1r = csrr(kflag)
      ascle = bry(kflag)
      do 120 i=ib,n
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
        if (kflag.ge.3) go to 120
        str = abs(c2r)
        sti = abs(c2i)
        c2m = max(str,sti)
        if (c2m.le.ascle) go to 120
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
  120 continue
  160 continue
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
      inu = fnu
      fnf = fnu - inu
      ifn = inu + n - 1
      ang = fnf*sgn
      cspnr = cos(ang)
      cspni = sin(ang)
      if (mod(ifn,2).eq.0) go to 170
      cspnr = -cspnr
      cspni = -cspni
  170 continue
      asc = bry(1)
      iuf = 0
      kk = n
      kdflg = 1
      ib = ib - 1
      ic = ib - 1
      do 270 k=1,n
        fn = fnu + (kk-1)
c-----------------------------------------------------------------------
c     logic to sort out cases whose parameters were set for the k
c     function above
c-----------------------------------------------------------------------
        m=3
        if (n.gt.2) go to 175
  172   continue
        initd = init(j)
        phidr = phir(j)
        phidi = phii(j)
        zet1dr = zeta1r(j)
        zet1di = zeta1i(j)
        zet2dr = zeta2r(j)
        zet2di = zeta2i(j)
        sumdr = sumr(j)
        sumdi = sumi(j)
        m = j
        j = 3 - j
        go to 180
  175   continue
        if ((kk.eq.n).and.(ib.lt.n)) go to 180
        if ((kk.eq.ib).or.(kk.eq.ic)) go to 172
        initd = 0
  180   continue
        call zunik(zrr, zri, fn, 1, 0, tol, initd, phidr, phidi,
     *   zet1dr, zet1di, zet2dr, zet2di, sumdr, sumdi,
     *   cwrkr(1,m), cwrki(1,m))
        if (kode.eq.1) go to 200
        str = zrr + zet2dr
        sti = zri + zet2di
        rast = fn/zabs(str,sti)
        str = str*rast*rast
        sti = -sti*rast*rast
        s1r = -zet1dr + str
        s1i = -zet1di + sti
        go to 210
  200   continue
        s1r = -zet1dr + zet2dr
        s1i = -zet1di + zet2di
  210   continue
c-----------------------------------------------------------------------
c     test for underflow and overflow
c-----------------------------------------------------------------------
        rs1 = s1r
        if (abs(rs1).gt.elim) go to 260
        if (kdflg.eq.1) iflag = 2
        if (abs(rs1).lt.alim) go to 220
c-----------------------------------------------------------------------
c     refine  test and scale
c-----------------------------------------------------------------------
        aphi = zabs(phidr,phidi)
        rs1 = rs1 + log(aphi)
        if (abs(rs1).gt.elim) go to 260
        if (kdflg.eq.1) iflag = 1
        if (rs1.lt.0.0d0) go to 220
        if (kdflg.eq.1) iflag = 3
  220   continue
        str = phidr*sumdr - phidi*sumdi
        sti = phidr*sumdi + phidi*sumdr
        s2r = -csgni*sti
        s2i = csgni*str
        str = exp(s1r)*cssr(iflag)
        s1r = str*cos(s1i)
        s1i = str*sin(s1i)
        str = s2r*s1r - s2i*s1i
        s2i = s2r*s1i + s2i*s1r
        s2r = str
        if (iflag.ne.1) go to 230
        call zuchk(s2r, s2i, nw, bry(1), tol)
        if (nw.eq.0) go to 230
        s2r = zeror
        s2i = zeroi
  230   continue
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
        if (kode.eq.1) go to 250
        call zs1s2(zrr, zri, s1r, s1i, s2r, s2i, nw, asc, alim, iuf)
        nz = nz + nw
  250   continue
        yr(kk) = s1r*cspnr - s1i*cspni + s2r
        yi(kk) = cspnr*s1i + cspni*s1r + s2i
        kk = kk - 1
        cspnr = -cspnr
        cspni = -cspni
        if (c2r.ne.0.0d0 .or. c2i.ne.0.0d0) go to 255
        kdflg = 1
        go to 270
  255   continue
        if (kdflg.eq.2) go to 275
        kdflg = 2
        go to 270
  260   continue
        if (rs1.gt.0.0d0) go to 300
        s2r = zeror
        s2i = zeroi
        go to 230
  270 continue
      k = n
  275 continue
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
      do 290 i=1,il
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
        if (kode.eq.1) go to 280
        call zs1s2(zrr, zri, c1r, c1i, c2r, c2i, nw, asc, alim, iuf)
        nz = nz + nw
  280   continue
        yr(kk) = c1r*cspnr - c1i*cspni + c2r
        yi(kk) = c1r*cspni + c1i*cspnr + c2i
        kk = kk - 1
        cspnr = -cspnr
        cspni = -cspni
        if (iflag.ge.3) go to 290
        c2r = abs(ckr)
        c2i = abs(cki)
        c2m = max(c2r,c2i)
        if (c2m.le.ascle) go to 290
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
  290 continue
      return
  300 continue
      nz = -1
      return
      end
