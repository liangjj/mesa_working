*deck zacon
      subroutine zacon (zr, zi, fnu, kode, mr, n, yr, yi, nz, rl, fnul,
     +   tol, elim, alim)
c***begin prologue  zacon
c***subsidiary
c***purpose  subsidiary to zbesh and zbesk
c***library   slatec
c***type      all (cacon-a, zacon-a)
c***author  amos, d. e., (snl)
c***description
c
c     zacon applies the analytic continuation formula
c
c         k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn)
c                 mp=pi*mr*cmplx(0.0,1.0)
c
c     to continue the k function from the right half to the left
c     half z plane
c
c***see also  zbesh, zbesk
c***routines called  d1mach, zabs, zbinu, zbknu, zmlt, zs1s2
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zacon
c     complex ck,cone,cscl,cscr,csgn,cspn,cy,czero,c1,c2,rz,sc1,sc2,st,
c    *s1,s2,y,z,zn
      double precision alim, arg, ascle, as2, azn, bry, bscle, cki,
     * ckr, coner, cpn, cscl, cscr, csgni, csgnr, cspni, cspnr,
     * csr, csrr, cssr, cyi, cyr, c1i, c1m, c1r, c2i, c2r, elim, fmr,
     * fn, fnu, fnul, pi, pti, ptr, razn, rl, rzi, rzr, sc1i, sc1r,
     * sc2i, sc2r, sgn, spn, sti, str, s1i, s1r, s2i, s2r, tol, yi, yr,
     * yy, zeror, zi, zni, znr, zr, d1mach, zabs
      integer i, inu, iuf, kflag, kode, mr, n, nn, nw, nz
      dimension yr(n), yi(n), cyr(2), cyi(2), cssr(3), csrr(3), bry(3)
      external zabs
      data pi / 3.14159265358979324d0 /
      data zeror,coner / 0.0d0,1.0d0 /
c***first executable statement  zacon
      nz = 0
      znr = -zr
      zni = -zi
      nn = n
      call zbinu(znr, zni, fnu, kode, nn, yr, yi, nw, rl, fnul, tol,
     * elim, alim)
      if (nw.lt.0) go to 90
c-----------------------------------------------------------------------
c     analytic continuation to the left half plane for the k function
c-----------------------------------------------------------------------
      nn = min(2,n)
      call zbknu(znr, zni, fnu, kode, nn, cyr, cyi, nw, tol, elim, alim)
      if (nw.ne.0) go to 90
      s1r = cyr(1)
      s1i = cyi(1)
      fmr = mr
      sgn = -dsign(pi,fmr)
      csgnr = zeror
      csgni = sgn
      if (kode.eq.1) go to 10
      yy = -zni
      cpn = cos(yy)
      spn = sin(yy)
      call zmlt(csgnr, csgni, cpn, spn, csgnr, csgni)
   10 continue
c-----------------------------------------------------------------------
c     calculate cspn=exp(fnu*pi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*sgn
      cpn = cos(arg)
      spn = sin(arg)
      cspnr = cpn
      cspni = spn
      if (mod(inu,2).eq.0) go to 20
      cspnr = -cspnr
      cspni = -cspni
   20 continue
      iuf = 0
      c1r = s1r
      c1i = s1i
      c2r = yr(1)
      c2i = yi(1)
      ascle = 1.0d+3*d1mach(1)/tol
      if (kode.eq.1) go to 30
      call zs1s2(znr, zni, c1r, c1i, c2r, c2i, nw, ascle, alim, iuf)
      nz = nz + nw
      sc1r = c1r
      sc1i = c1i
   30 continue
      call zmlt(cspnr, cspni, c1r, c1i, str, sti)
      call zmlt(csgnr, csgni, c2r, c2i, ptr, pti)
      yr(1) = str + ptr
      yi(1) = sti + pti
      if (n.eq.1) return
      cspnr = -cspnr
      cspni = -cspni
      s2r = cyr(2)
      s2i = cyi(2)
      c1r = s2r
      c1i = s2i
      c2r = yr(2)
      c2i = yi(2)
      if (kode.eq.1) go to 40
      call zs1s2(znr, zni, c1r, c1i, c2r, c2i, nw, ascle, alim, iuf)
      nz = nz + nw
      sc2r = c1r
      sc2i = c1i
   40 continue
      call zmlt(cspnr, cspni, c1r, c1i, str, sti)
      call zmlt(csgnr, csgni, c2r, c2i, ptr, pti)
      yr(2) = str + ptr
      yi(2) = sti + pti
      if (n.eq.2) return
      cspnr = -cspnr
      cspni = -cspni
      azn = zabs(znr,zni)
      razn = 1.0d0/azn
      str = znr*razn
      sti = -zni*razn
      rzr = (str+str)*razn
      rzi = (sti+sti)*razn
      fn = fnu + 1.0d0
      ckr = fn*rzr
      cki = fn*rzi
c-----------------------------------------------------------------------
c     scale near exponent extremes during recurrence on k functions
c-----------------------------------------------------------------------
      cscl = 1.0d0/tol
      cscr = tol
      cssr(1) = cscl
      cssr(2) = coner
      cssr(3) = cscr
      csrr(1) = cscr
      csrr(2) = coner
      csrr(3) = cscl
      bry(1) = ascle
      bry(2) = 1.0d0/ascle
      bry(3) = d1mach(2)
      as2 = zabs(s2r,s2i)
      kflag = 2
      if (as2.gt.bry(1)) go to 50
      kflag = 1
      go to 60
   50 continue
      if (as2.lt.bry(2)) go to 60
      kflag = 3
   60 continue
      bscle = bry(kflag)
      s1r = s1r*cssr(kflag)
      s1i = s1i*cssr(kflag)
      s2r = s2r*cssr(kflag)
      s2i = s2i*cssr(kflag)
      csr = csrr(kflag)
      do 80 i=3,n
        str = s2r
        sti = s2i
        s2r = ckr*str - cki*sti + s1r
        s2i = ckr*sti + cki*str + s1i
        s1r = str
        s1i = sti
        c1r = s2r*csr
        c1i = s2i*csr
        str = c1r
        sti = c1i
        c2r = yr(i)
        c2i = yi(i)
        if (kode.eq.1) go to 70
        if (iuf.lt.0) go to 70
        call zs1s2(znr, zni, c1r, c1i, c2r, c2i, nw, ascle, alim, iuf)
        nz = nz + nw
        sc1r = sc2r
        sc1i = sc2i
        sc2r = c1r
        sc2i = c1i
        if (iuf.ne.3) go to 70
        iuf = -4
        s1r = sc1r*cssr(kflag)
        s1i = sc1i*cssr(kflag)
        s2r = sc2r*cssr(kflag)
        s2i = sc2i*cssr(kflag)
        str = sc2r
        sti = sc2i
   70   continue
        ptr = cspnr*c1r - cspni*c1i
        pti = cspnr*c1i + cspni*c1r
        yr(i) = ptr + csgnr*c2r - csgni*c2i
        yi(i) = pti + csgnr*c2i + csgni*c2r
        ckr = ckr + rzr
        cki = cki + rzi
        cspnr = -cspnr
        cspni = -cspni
        if (kflag.ge.3) go to 80
        ptr = abs(c1r)
        pti = abs(c1i)
        c1m = max(ptr,pti)
        if (c1m.le.bscle) go to 80
        kflag = kflag + 1
        bscle = bry(kflag)
        s1r = s1r*csr
        s1i = s1i*csr
        s2r = str
        s2i = sti
        s1r = s1r*cssr(kflag)
        s1i = s1i*cssr(kflag)
        s2r = s2r*cssr(kflag)
        s2i = s2i*cssr(kflag)
        csr = csrr(kflag)
   80 continue
      return
   90 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
