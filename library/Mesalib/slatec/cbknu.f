*deck cbknu
      subroutine cbknu (z, fnu, kode, n, y, nz, tol, elim, alim)
c***begin prologue  cbknu
c***subsidiary
c***purpose  subsidiary to cairy, cbesh, cbesi and cbesk
c***library   slatec
c***type      all (cbknu-a, zbknu-a)
c***author  amos, d. e., (snl)
c***description
c
c     cbknu computes the k bessel function in the right half z plane
c
c***see also  cairy, cbesh, cbesi, cbesk
c***routines called  ckscl, cshch, cuchk, gamln, i1mach, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cbknu
c
      complex cch, ck, coef, cone, crsc, cs, cscl, csh, csr, css, ctwo,
     * cz, czero, f, fmu, p, pt, p1, p2, q, rz, smu, st, s1, s2, y, z,
     * zd, celm, cy
      real aa, ak, alim, ascle, a1, a2, bb, bk, bry, caz, cc, dnu,
     * dnu2, elim, etest, fc, fhs, fk, fks, fnu, fpi, g1, g2, hpi, pi,
     * p2i, p2m, p2r, rk, rthpi, r1, s, spi, tm, tol, tth, t1, t2, xx,
     * yy, gamln, r1mach, helim, elm, xd, yd, alas, as
      integer i, idum, iflag, inu, k, kflag, kk, kmax, kode, koded, n,
     * nz, i1mach, nw, j, ic, inub
      dimension bry(3), cc(8), css(3), csr(3), y(n), cy(2)
c
      data kmax / 30 /
      data r1 / 2.0e0 /
      data czero,cone,ctwo /(0.0e0,0.0e0),(1.0e0,0.0e0),(2.0e0,0.0e0)/
c
      data pi, rthpi, spi ,hpi, fpi, tth /
     1     3.14159265358979324e0,       1.25331413731550025e0,
     2     1.90985931710274403e0,       1.57079632679489662e0,
     3     1.89769999331517738e0,       6.66666666666666666e-01/
c
      data cc(1), cc(2), cc(3), cc(4), cc(5), cc(6), cc(7), cc(8)/
     1     5.77215664901532861e-01,    -4.20026350340952355e-02,
     2    -4.21977345555443367e-02,     7.21894324666309954e-03,
     3    -2.15241674114950973e-04,    -2.01348547807882387e-05,
     4     1.13302723198169588e-06,     6.11609510448141582e-09/
c
c***first executable statement  cbknu
      xx = real(z)
      yy = aimag(z)
      caz = abs(z)
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
      nz = 0
      iflag = 0
      koded = kode
      rz = ctwo/z
      inu = fnu+0.5e0
      dnu = fnu - inu
      if (abs(dnu).eq.0.5e0) go to 110
      dnu2 = 0.0e0
      if (abs(dnu).gt.tol) dnu2 = dnu*dnu
      if (caz.gt.r1) go to 110
c-----------------------------------------------------------------------
c     series for abs(z).le.r1
c-----------------------------------------------------------------------
      fc = 1.0e0
      smu = clog(rz)
      fmu = smu*cmplx(dnu,0.0e0)
      call cshch(fmu, csh, cch)
      if (dnu.eq.0.0e0) go to 10
      fc = dnu*pi
      fc = fc/sin(fc)
      smu = csh*cmplx(1.0e0/dnu,0.0e0)
   10 continue
      a2 = 1.0e0 + dnu
c-----------------------------------------------------------------------
c     gam(1-z)*gam(1+z)=pi*z/sin(pi*z), t1=1/gam(1-dnu), t2=1/gam(1+dnu)
c-----------------------------------------------------------------------
      t2 = exp(-gamln(a2,idum))
      t1 = 1.0e0/(t2*fc)
      if (abs(dnu).gt.0.1e0) go to 40
c-----------------------------------------------------------------------
c     series for f0 to resolve indeterminacy for small abs(dnu)
c-----------------------------------------------------------------------
      ak = 1.0e0
      s = cc(1)
      do 20 k=2,8
        ak = ak*dnu2
        tm = cc(k)*ak
        s = s + tm
        if (abs(tm).lt.tol) go to 30
   20 continue
   30 g1 = -s
      go to 50
   40 continue
      g1 = (t1-t2)/(dnu+dnu)
   50 continue
      g2 = 0.5e0*(t1+t2)*fc
      g1 = g1*fc
      f = cmplx(g1,0.0e0)*cch + smu*cmplx(g2,0.0e0)
      pt = cexp(fmu)
      p = cmplx(0.5e0/t2,0.0e0)*pt
      q = cmplx(0.5e0/t1,0.0e0)/pt
      s1 = f
      s2 = p
      ak = 1.0e0
      a1 = 1.0e0
      ck = cone
      bk = 1.0e0 - dnu2
      if (inu.gt.0 .or. n.gt.1) go to 80
c-----------------------------------------------------------------------
c     generate k(fnu,z), 0.0d0 .le. fnu .lt. 0.5d0 and n=1
c-----------------------------------------------------------------------
      if (caz.lt.tol) go to 70
      cz = z*z*cmplx(0.25e0,0.0e0)
      t1 = 0.25e0*caz*caz
   60 continue
      f = (f*cmplx(ak,0.0e0)+p+q)*cmplx(1.0e0/bk,0.0e0)
      p = p*cmplx(1.0e0/(ak-dnu),0.0e0)
      q = q*cmplx(1.0e0/(ak+dnu),0.0e0)
      rk = 1.0e0/ak
      ck = ck*cz*cmplx(rk,0.0)
      s1 = s1 + ck*f
      a1 = a1*t1*rk
      bk = bk + ak + ak + 1.0e0
      ak = ak + 1.0e0
      if (a1.gt.tol) go to 60
   70 continue
      y(1) = s1
      if (koded.eq.1) return
      y(1) = s1*cexp(z)
      return
c-----------------------------------------------------------------------
c     generate k(dnu,z) and k(dnu+1,z) for forward recurrence
c-----------------------------------------------------------------------
   80 continue
      if (caz.lt.tol) go to 100
      cz = z*z*cmplx(0.25e0,0.0e0)
      t1 = 0.25e0*caz*caz
   90 continue
      f = (f*cmplx(ak,0.0e0)+p+q)*cmplx(1.0e0/bk,0.0e0)
      p = p*cmplx(1.0e0/(ak-dnu),0.0e0)
      q = q*cmplx(1.0e0/(ak+dnu),0.0e0)
      rk = 1.0e0/ak
      ck = ck*cz*cmplx(rk,0.0e0)
      s1 = s1 + ck*f
      s2 = s2 + ck*(p-f*cmplx(ak,0.0e0))
      a1 = a1*t1*rk
      bk = bk + ak + ak + 1.0e0
      ak = ak + 1.0e0
      if (a1.gt.tol) go to 90
  100 continue
      kflag = 2
      bk = real(smu)
      a1 = fnu + 1.0e0
      ak = a1*abs(bk)
      if (ak.gt.alim) kflag = 3
      p2 = s2*css(kflag)
      s2 = p2*rz
      s1 = s1*css(kflag)
      if (koded.eq.1) go to 210
      f = cexp(z)
      s1 = s1*f
      s2 = s2*f
      go to 210
c-----------------------------------------------------------------------
c     iflag=0 means no underflow occurred
c     iflag=1 means an underflow occurred- computation proceeds with
c     koded=2 and a test for on scale values is made during forward
c     recursion
c-----------------------------------------------------------------------
  110 continue
      coef = cmplx(rthpi,0.0e0)/csqrt(z)
      kflag = 2
      if (koded.eq.2) go to 120
      if (xx.gt.alim) go to 290
c     blank line
      a1 = exp(-xx)*real(css(kflag))
      pt = cmplx(a1,0.0e0)*cmplx(cos(yy),-sin(yy))
      coef = coef*pt
  120 continue
      if (abs(dnu).eq.0.5e0) go to 300
c-----------------------------------------------------------------------
c     miller algorithm for abs(z).gt.r1
c-----------------------------------------------------------------------
      ak = cos(pi*dnu)
      ak = abs(ak)
      if (ak.eq.0.0e0) go to 300
      fhs = abs(0.25e0-dnu2)
      if (fhs.eq.0.0e0) go to 300
c-----------------------------------------------------------------------
c     compute r2=f(e). if abs(z).ge.r2, use forward recurrence to
c     determine the backward index k. r2=f(e) is a straight line on
c     12.le.e.le.60. e is computed from 2**(-e)=b**(1-i1mach(11))=
c     tol where b is the base of the arithmetic.
c-----------------------------------------------------------------------
      t1 = (i1mach(11)-1)*r1mach(5)*3.321928094e0
      t1 = max(t1,12.0e0)
      t1 = min(t1,60.0e0)
      t2 = tth*t1 - 6.0e0
      if (xx.ne.0.0e0) go to 130
      t1 = hpi
      go to 140
  130 continue
      t1 = atan(yy/xx)
      t1 = abs(t1)
  140 continue
      if (t2.gt.caz) go to 170
c-----------------------------------------------------------------------
c     forward recurrence loop when abs(z).ge.r2
c-----------------------------------------------------------------------
      etest = ak/(pi*caz*tol)
      fk = 1.0e0
      if (etest.lt.1.0e0) go to 180
      fks = 2.0e0
      rk = caz + caz + 2.0e0
      a1 = 0.0e0
      a2 = 1.0e0
      do 150 i=1,kmax
        ak = fhs/fks
        bk = rk/(fk+1.0e0)
        tm = a2
        a2 = bk*a2 - ak*a1
        a1 = tm
        rk = rk + 2.0e0
        fks = fks + fk + fk + 2.0e0
        fhs = fhs + fk + fk
        fk = fk + 1.0e0
        tm = abs(a2)*fk
        if (etest.lt.tm) go to 160
  150 continue
      go to 310
  160 continue
      fk = fk + spi*t1*sqrt(t2/caz)
      fhs = abs(0.25e0-dnu2)
      go to 180
  170 continue
c-----------------------------------------------------------------------
c     compute backward index k for abs(z).lt.r2
c-----------------------------------------------------------------------
      a2 = sqrt(caz)
      ak = fpi*ak/(tol*sqrt(a2))
      aa = 3.0e0*t1/(1.0e0+caz)
      bb = 14.7e0*t1/(28.0e0+caz)
      ak = (alog(ak)+caz*cos(aa)/(1.0e0+0.008e0*caz))/cos(bb)
      fk = 0.12125e0*ak*ak/caz + 1.5e0
  180 continue
      k = fk
c-----------------------------------------------------------------------
c     backward recurrence loop for miller algorithm
c-----------------------------------------------------------------------
      fk = k
      fks = fk*fk
      p1 = czero
      p2 = cmplx(tol,0.0e0)
      cs = p2
      do 190 i=1,k
        a1 = fks - fk
        a2 = (fks+fk)/(a1+fhs)
        rk = 2.0e0/(fk+1.0e0)
        t1 = (fk+xx)*rk
        t2 = yy*rk
        pt = p2
        p2 = (p2*cmplx(t1,t2)-p1)*cmplx(a2,0.0e0)
        p1 = pt
        cs = cs + p2
        fks = a1 - fk + 1.0e0
        fk = fk - 1.0e0
  190 continue
c-----------------------------------------------------------------------
c     compute (p2/cs)=(p2/abs(cs))*(conjg(cs)/abs(cs)) for better
c     scaling
c-----------------------------------------------------------------------
      tm = abs(cs)
      pt = cmplx(1.0e0/tm,0.0e0)
      s1 = pt*p2
      cs = conjg(cs)*pt
      s1 = coef*s1*cs
      if (inu.gt.0 .or. n.gt.1) go to 200
      zd = z
      if(iflag.eq.1) go to 270
      go to 240
  200 continue
c-----------------------------------------------------------------------
c     compute p1/p2=(p1/abs(p2)*conjg(p2)/abs(p2) for scaling
c-----------------------------------------------------------------------
      tm = abs(p2)
      pt = cmplx(1.0e0/tm,0.0e0)
      p1 = pt*p1
      p2 = conjg(p2)*pt
      pt = p1*p2
      s2 = s1*(cone+(cmplx(dnu+0.5e0,0.0e0)-pt)/z)
c-----------------------------------------------------------------------
c     forward recursion on the three term recursion relation with
c     scaling near exponent extremes on kflag=1 or kflag=3
c-----------------------------------------------------------------------
  210 continue
      ck = cmplx(dnu+1.0e0,0.0e0)*rz
      if (n.eq.1) inu = inu - 1
      if (inu.gt.0) go to 220
      if (n.eq.1) s1=s2
      zd = z
      if(iflag.eq.1) go to 270
      go to 240
  220 continue
      inub = 1
      if (iflag.eq.1) go to 261
  225 continue
      p1 = csr(kflag)
      ascle = bry(kflag)
      do 230 i=inub,inu
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        ck = ck + rz
        if (kflag.ge.3) go to 230
        p2 = s2*p1
        p2r = real(p2)
        p2i = aimag(p2)
        p2r = abs(p2r)
        p2i = abs(p2i)
        p2m = max(p2r,p2i)
        if (p2m.le.ascle) go to 230
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*p1
        s2 = p2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        p1 = csr(kflag)
  230 continue
      if (n.eq.1) s1 = s2
  240 continue
      y(1) = s1*csr(kflag)
      if (n.eq.1) return
      y(2) = s2*csr(kflag)
      if (n.eq.2) return
      kk = 2
  250 continue
      kk = kk + 1
      if (kk.gt.n) return
      p1 = csr(kflag)
      ascle = bry(kflag)
      do 260 i=kk,n
        p2 = s2
        s2 = ck*s2 + s1
        s1 = p2
        ck = ck + rz
        p2 = s2*p1
        y(i) = p2
        if (kflag.ge.3) go to 260
        p2r = real(p2)
        p2i = aimag(p2)
        p2r = abs(p2r)
        p2i = abs(p2i)
        p2m = max(p2r,p2i)
        if (p2m.le.ascle) go to 260
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*p1
        s2 = p2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        p1 = csr(kflag)
  260 continue
      return
c-----------------------------------------------------------------------
c     iflag=1 cases, forward recurrence on scaled values on underflow
c-----------------------------------------------------------------------
  261 continue
      helim = 0.5e0*elim
      elm = exp(-elim)
      celm = cmplx(elm,0.0)
      ascle = bry(1)
      zd = z
      xd = xx
      yd = yy
      ic = -1
      j = 2
      do 262 i=1,inu
        st = s2
        s2 = ck*s2+s1
        s1 = st
        ck = ck+rz
        as = abs(s2)
        alas = alog(as)
        p2r = -xd+alas
        if(p2r.lt.(-elim)) go to 263
        p2 = -zd+clog(s2)
        p2r = real(p2)
        p2i = aimag(p2)
        p2m = exp(p2r)/tol
        p1 = cmplx(p2m,0.0e0)*cmplx(cos(p2i),sin(p2i))
        call cuchk(p1,nw,ascle,tol)
        if(nw.ne.0) go to 263
        j=3-j
        cy(j) = p1
        if(ic.eq.(i-1)) go to 264
        ic = i
        go to 262
  263   continue
        if(alas.lt.helim) go to 262
        xd = xd-elim
        s1 = s1*celm
        s2 = s2*celm
        zd = cmplx(xd,yd)
  262 continue
      if(n.eq.1) s1 = s2
      go to 270
  264 continue
      kflag = 1
      inub = i+1
      s2 = cy(j)
      j = 3 - j
      s1 = cy(j)
      if(inub.le.inu) go to 225
      if(n.eq.1) s1 = s2
      go to 240
  270 continue
      y(1) = s1
      if (n.eq.1) go to 280
      y(2) = s2
  280 continue
      ascle = bry(1)
      call ckscl(zd, fnu, n, y, nz, rz, ascle, tol, elim)
      inu = n - nz
      if (inu.le.0) return
      kk = nz + 1
      s1 = y(kk)
      y(kk) = s1*csr(1)
      if (inu.eq.1) return
      kk = nz + 2
      s2 = y(kk)
      y(kk) = s2*csr(1)
      if (inu.eq.2) return
      t2 = fnu + (kk-1)
      ck = cmplx(t2,0.0e0)*rz
      kflag = 1
      go to 250
  290 continue
c-----------------------------------------------------------------------
c     scale by exp(z), iflag = 1 cases
c-----------------------------------------------------------------------
      koded = 2
      iflag = 1
      kflag = 2
      go to 120
c-----------------------------------------------------------------------
c     fnu=half odd integer case, dnu=-0.5
c-----------------------------------------------------------------------
  300 continue
      s1 = coef
      s2 = coef
      go to 210
  310 continue
      nz=-2
      return
      end
