*deck cacon
      subroutine cacon (z, fnu, kode, mr, n, y, nz, rl, fnul, tol, elim,
     +   alim)
c***begin prologue  cacon
c***subsidiary
c***purpose  subsidiary to cbesh and cbesk
c***library   slatec
c***type      all (cacon-a, zacon-a)
c***author  amos, d. e., (snl)
c***description
c
c     cacon applies the analytic continuation formula
c
c         k(fnu,zn*exp(mp))=k(fnu,zn)*exp(-mp*fnu) - mp*i(fnu,zn)
c                 mp=pi*mr*cmplx(0.0,1.0)
c
c     to continue the k function from the right half to the left
c     half z plane
c
c***see also  cbesh, cbesk
c***routines called  cbinu, cbknu, cs1s2, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cacon
      complex ck, cone, cs, cscl, cscr, csgn, cspn, css, csr, c1, c2,
     * rz, sc1, sc2, st, s1, s2, y, z, zn, cy
      real alim, arg, ascle, as2, bscle, bry, cpn, c1i, c1m, c1r, elim,
     * fmr, fnu, fnul, pi, rl, sgn, spn, tol, yy, r1mach
      integer i, inu, iuf, kflag, kode, mr, n, nn, nw, nz
      dimension y(n), cy(2), css(3), csr(3), bry(3)
      data pi / 3.14159265358979324e0 /
      data cone / (1.0e0,0.0e0) /
c***first executable statement  cacon
      nz = 0
      zn = -z
      nn = n
      call cbinu(zn, fnu, kode, nn, y, nw, rl, fnul, tol, elim, alim)
      if (nw.lt.0) go to 80
c-----------------------------------------------------------------------
c     analytic continuation to the left half plane for the k function
c-----------------------------------------------------------------------
      nn = min(2,n)
      call cbknu(zn, fnu, kode, nn, cy, nw, tol, elim, alim)
      if (nw.ne.0) go to 80
      s1 = cy(1)
      fmr = mr
      sgn = -sign(pi,fmr)
      csgn = cmplx(0.0e0,sgn)
      if (kode.eq.1) go to 10
      yy = -aimag(zn)
      cpn = cos(yy)
      spn = sin(yy)
      csgn = csgn*cmplx(cpn,spn)
   10 continue
c-----------------------------------------------------------------------
c     calculate cspn=exp(fnu*pi*i) to minimize losses of significance
c     when fnu is large
c-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu-inu)*sgn
      cpn = cos(arg)
      spn = sin(arg)
      cspn = cmplx(cpn,spn)
      if (mod(inu,2).eq.1) cspn = -cspn
      iuf = 0
      c1 = s1
      c2 = y(1)
      ascle = 1.0e+3*r1mach(1)/tol
      if (kode.eq.1) go to 20
      call cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc1 = c1
   20 continue
      y(1) = cspn*c1 + csgn*c2
      if (n.eq.1) return
      cspn = -cspn
      s2 = cy(2)
      c1 = s2
      c2 = y(2)
      if (kode.eq.1) go to 30
      call cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc2 = c1
   30 continue
      y(2) = cspn*c1 + csgn*c2
      if (n.eq.2) return
      cspn = -cspn
      rz = cmplx(2.0e0,0.0e0)/zn
      ck = cmplx(fnu+1.0e0,0.0e0)*rz
c-----------------------------------------------------------------------
c     scale near exponent extremes during recurrence on k functions
c-----------------------------------------------------------------------
      cscl = cmplx(1.0e0/tol,0.0e0)
      cscr = cmplx(tol,0.0e0)
      css(1) = cscl
      css(2) = cone
      css(3) = cscr
      csr(1) = cscr
      csr(2) = cone
      csr(3) = cscl
      bry(1) = ascle
      bry(2) = 1.0e0/ascle
      bry(3) = r1mach(2)
      as2 = abs(s2)
      kflag = 2
      if (as2.gt.bry(1)) go to 40
      kflag = 1
      go to 50
   40 continue
      if (as2.lt.bry(2)) go to 50
      kflag = 3
   50 continue
      bscle = bry(kflag)
      s1 = s1*css(kflag)
      s2 = s2*css(kflag)
      cs = csr(kflag)
      do 70 i=3,n
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        c1 = s2*cs
        st = c1
        c2 = y(i)
        if (kode.eq.1) go to 60
        if (iuf.lt.0) go to 60
        call cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
        nz = nz + nw
        sc1 = sc2
        sc2 = c1
        if (iuf.ne.3) go to 60
        iuf = -4
        s1 = sc1*css(kflag)
        s2 = sc2*css(kflag)
        st = sc2
   60   continue
        y(i) = cspn*c1 + csgn*c2
        ck = ck + rz
        cspn = -cspn
        if (kflag.ge.3) go to 70
        c1r = real(c1)
        c1i = aimag(c1)
        c1r = abs(c1r)
        c1i = abs(c1i)
        c1m = max(c1r,c1i)
        if (c1m.le.bscle) go to 70
        kflag = kflag + 1
        bscle = bry(kflag)
        s1 = s1*cs
        s2 = st
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        cs = csr(kflag)
   70 continue
      return
   80 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
      end
