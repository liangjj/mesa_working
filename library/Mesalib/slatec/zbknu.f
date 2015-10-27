*deck zbknu
      subroutine zbknu (zr, zi, fnu, kode, n, yr, yi, nz, tol, elim,
     +   alim)
c***begin prologue  zbknu
c***subsidiary
c***purpose  subsidiary to zairy, zbesh, zbesi and zbesk
c***library   slatec
c***type      all (cbknu-a, zbknu-a)
c***author  amos, d. e., (snl)
c***description
c
c     zbknu computes the k bessel function in the right half z plane.
c
c***see also  zairy, zbesh, zbesi, zbesk
c***routines called  d1mach, dgamln, i1mach, zabs, zdiv, zexp, zkscl,
c                    zlog, zmlt, zshch, zsqrt, zuchk
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zexp, zlog and zsqrt to external statement.  (rwc)
c***end prologue  zbknu
c
      double precision aa, ak, alim, ascle, a1, a2, bb, bk, bry, caz,
     * cbi, cbr, cc, cchi, cchr, cki, ckr, coefi, coefr, conei, coner,
     * crscr, csclr, cshi, cshr, csi, csr, csrr, cssr, ctwor,
     * czeroi, czeror, czi, czr, dnu, dnu2, dpi, elim, etest, fc, fhs,
     * fi, fk, fks, fmui, fmur, fnu, fpi, fr, g1, g2, hpi, pi, pr, pti,
     * ptr, p1i, p1r, p2i, p2m, p2r, qi, qr, rak, rcaz, rthpi, rzi,
     * rzr, r1, s, smui, smur, spi, sti, str, s1i, s1r, s2i, s2r, tm,
     * tol, tth, t1, t2, yi, yr, zi, zr, dgamln, d1mach, zabs, elm,
     * celmr, zdr, zdi, as, alas, helim, cyr, cyi
      integer i, iflag, inu, k, kflag, kk, kmax, kode, koded, n, nz,
     * idum, i1mach, j, ic, inub, nw
      dimension yr(n), yi(n), cc(8), cssr(3), csrr(3), bry(3), cyr(2),
     * cyi(2)
      external zabs, zexp, zlog, zsqrt
c     complex z,y,a,b,rz,smu,fu,fmu,f,flrz,cz,s1,s2,csh,cch
c     complex ck,p,q,coef,p1,p2,cbk,pt,czero,cone,ctwo,st,ez,cs,dk
c
      data kmax / 30 /
      data czeror,czeroi,coner,conei,ctwor,r1/
     1  0.0d0 , 0.0d0 , 1.0d0 , 0.0d0 , 2.0d0 , 2.0d0 /
      data dpi, rthpi, spi ,hpi, fpi, tth /
     1     3.14159265358979324d0,       1.25331413731550025d0,
     2     1.90985931710274403d0,       1.57079632679489662d0,
     3     1.89769999331517738d0,       6.66666666666666666d-01/
      data cc(1), cc(2), cc(3), cc(4), cc(5), cc(6), cc(7), cc(8)/
     1     5.77215664901532861d-01,    -4.20026350340952355d-02,
     2    -4.21977345555443367d-02,     7.21894324666309954d-03,
     3    -2.15241674114950973d-04,    -2.01348547807882387d-05,
     4     1.13302723198169588d-06,     6.11609510448141582d-09/
c***first executable statement  zbknu
      caz = zabs(zr,zi)
      csclr = 1.0d0/tol
      crscr = tol
      cssr(1) = csclr
      cssr(2) = 1.0d0
      cssr(3) = crscr
      csrr(1) = crscr
      csrr(2) = 1.0d0
      csrr(3) = csclr
      bry(1) = 1.0d+3*d1mach(1)/tol
      bry(2) = 1.0d0/bry(1)
      bry(3) = d1mach(2)
      nz = 0
      iflag = 0
      koded = kode
      rcaz = 1.0d0/caz
      str = zr*rcaz
      sti = -zi*rcaz
      rzr = (str+str)*rcaz
      rzi = (sti+sti)*rcaz
      inu = fnu+0.5d0
      dnu = fnu - inu
      if (abs(dnu).eq.0.5d0) go to 110
      dnu2 = 0.0d0
      if (abs(dnu).gt.tol) dnu2 = dnu*dnu
      if (caz.gt.r1) go to 110
c-----------------------------------------------------------------------
c     series for abs(z).le.r1
c-----------------------------------------------------------------------
      fc = 1.0d0
      call zlog(rzr, rzi, smur, smui, idum)
      fmur = smur*dnu
      fmui = smui*dnu
      call zshch(fmur, fmui, cshr, cshi, cchr, cchi)
      if (dnu.eq.0.0d0) go to 10
      fc = dnu*dpi
      fc = fc/sin(fc)
      smur = cshr/dnu
      smui = cshi/dnu
   10 continue
      a2 = 1.0d0 + dnu
c-----------------------------------------------------------------------
c     gam(1-z)*gam(1+z)=pi*z/sin(pi*z), t1=1/gam(1-dnu), t2=1/gam(1+dnu)
c-----------------------------------------------------------------------
      t2 = exp(-dgamln(a2,idum))
      t1 = 1.0d0/(t2*fc)
      if (abs(dnu).gt.0.1d0) go to 40
c-----------------------------------------------------------------------
c     series for f0 to resolve indeterminacy for small abs(dnu)
c-----------------------------------------------------------------------
      ak = 1.0d0
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
      g2 = (t1+t2)*0.5d0
      fr = fc*(cchr*g1+smur*g2)
      fi = fc*(cchi*g1+smui*g2)
      call zexp(fmur, fmui, str, sti)
      pr = 0.5d0*str/t2
      pi = 0.5d0*sti/t2
      call zdiv(0.5d0, 0.0d0, str, sti, ptr, pti)
      qr = ptr/t1
      qi = pti/t1
      s1r = fr
      s1i = fi
      s2r = pr
      s2i = pi
      ak = 1.0d0
      a1 = 1.0d0
      ckr = coner
      cki = conei
      bk = 1.0d0 - dnu2
      if (inu.gt.0 .or. n.gt.1) go to 80
c-----------------------------------------------------------------------
c     generate k(fnu,z), 0.0d0 .le. fnu .lt. 0.5d0 and n=1
c-----------------------------------------------------------------------
      if (caz.lt.tol) go to 70
      call zmlt(zr, zi, zr, zi, czr, czi)
      czr = 0.25d0*czr
      czi = 0.25d0*czi
      t1 = 0.25d0*caz*caz
   60 continue
      fr = (fr*ak+pr+qr)/bk
      fi = (fi*ak+pi+qi)/bk
      str = 1.0d0/(ak-dnu)
      pr = pr*str
      pi = pi*str
      str = 1.0d0/(ak+dnu)
      qr = qr*str
      qi = qi*str
      str = ckr*czr - cki*czi
      rak = 1.0d0/ak
      cki = (ckr*czi+cki*czr)*rak
      ckr = str*rak
      s1r = ckr*fr - cki*fi + s1r
      s1i = ckr*fi + cki*fr + s1i
      a1 = a1*t1*rak
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      if (a1.gt.tol) go to 60
   70 continue
      yr(1) = s1r
      yi(1) = s1i
      if (koded.eq.1) return
      call zexp(zr, zi, str, sti)
      call zmlt(s1r, s1i, str, sti, yr(1), yi(1))
      return
c-----------------------------------------------------------------------
c     generate k(dnu,z) and k(dnu+1,z) for forward recurrence
c-----------------------------------------------------------------------
   80 continue
      if (caz.lt.tol) go to 100
      call zmlt(zr, zi, zr, zi, czr, czi)
      czr = 0.25d0*czr
      czi = 0.25d0*czi
      t1 = 0.25d0*caz*caz
   90 continue
      fr = (fr*ak+pr+qr)/bk
      fi = (fi*ak+pi+qi)/bk
      str = 1.0d0/(ak-dnu)
      pr = pr*str
      pi = pi*str
      str = 1.0d0/(ak+dnu)
      qr = qr*str
      qi = qi*str
      str = ckr*czr - cki*czi
      rak = 1.0d0/ak
      cki = (ckr*czi+cki*czr)*rak
      ckr = str*rak
      s1r = ckr*fr - cki*fi + s1r
      s1i = ckr*fi + cki*fr + s1i
      str = pr - fr*ak
      sti = pi - fi*ak
      s2r = ckr*str - cki*sti + s2r
      s2i = ckr*sti + cki*str + s2i
      a1 = a1*t1*rak
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      if (a1.gt.tol) go to 90
  100 continue
      kflag = 2
      a1 = fnu + 1.0d0
      ak = a1*abs(smur)
      if (ak.gt.alim) kflag = 3
      str = cssr(kflag)
      p2r = s2r*str
      p2i = s2i*str
      call zmlt(p2r, p2i, rzr, rzi, s2r, s2i)
      s1r = s1r*str
      s1i = s1i*str
      if (koded.eq.1) go to 210
      call zexp(zr, zi, fr, fi)
      call zmlt(s1r, s1i, fr, fi, s1r, s1i)
      call zmlt(s2r, s2i, fr, fi, s2r, s2i)
      go to 210
c-----------------------------------------------------------------------
c     iflag=0 means no underflow occurred
c     iflag=1 means an underflow occurred- computation proceeds with
c     koded=2 and a test for on scale values is made during forward
c     recursion
c-----------------------------------------------------------------------
  110 continue
      call zsqrt(zr, zi, str, sti)
      call zdiv(rthpi, czeroi, str, sti, coefr, coefi)
      kflag = 2
      if (koded.eq.2) go to 120
      if (zr.gt.alim) go to 290
c     blank line
      str = exp(-zr)*cssr(kflag)
      sti = -str*sin(zi)
      str = str*cos(zi)
      call zmlt(coefr, coefi, str, sti, coefr, coefi)
  120 continue
      if (abs(dnu).eq.0.5d0) go to 300
c-----------------------------------------------------------------------
c     miller algorithm for abs(z).gt.r1
c-----------------------------------------------------------------------
      ak = cos(dpi*dnu)
      ak = abs(ak)
      if (ak.eq.czeror) go to 300
      fhs = abs(0.25d0-dnu2)
      if (fhs.eq.czeror) go to 300
c-----------------------------------------------------------------------
c     compute r2=f(e). if abs(z).ge.r2, use forward recurrence to
c     determine the backward index k. r2=f(e) is a straight line on
c     12.le.e.le.60. e is computed from 2**(-e)=b**(1-i1mach(14))=
c     tol where b is the base of the arithmetic.
c-----------------------------------------------------------------------
      t1 = i1mach(14)-1
      t1 = t1*d1mach(5)*3.321928094d0
      t1 = max(t1,12.0d0)
      t1 = min(t1,60.0d0)
      t2 = tth*t1 - 6.0d0
      if (zr.ne.0.0d0) go to 130
      t1 = hpi
      go to 140
  130 continue
      t1 = datan(zi/zr)
      t1 = abs(t1)
  140 continue
      if (t2.gt.caz) go to 170
c-----------------------------------------------------------------------
c     forward recurrence loop when abs(z).ge.r2
c-----------------------------------------------------------------------
      etest = ak/(dpi*caz*tol)
      fk = coner
      if (etest.lt.coner) go to 180
      fks = ctwor
      ckr = caz + caz + ctwor
      p1r = czeror
      p2r = coner
      do 150 i=1,kmax
        ak = fhs/fks
        cbr = ckr/(fk+coner)
        ptr = p2r
        p2r = cbr*p2r - p1r*ak
        p1r = ptr
        ckr = ckr + ctwor
        fks = fks + fk + fk + ctwor
        fhs = fhs + fk + fk
        fk = fk + coner
        str = abs(p2r)*fk
        if (etest.lt.str) go to 160
  150 continue
      go to 310
  160 continue
      fk = fk + spi*t1*sqrt(t2/caz)
      fhs = abs(0.25d0-dnu2)
      go to 180
  170 continue
c-----------------------------------------------------------------------
c     compute backward index k for abs(z).lt.r2
c-----------------------------------------------------------------------
      a2 = sqrt(caz)
      ak = fpi*ak/(tol*sqrt(a2))
      aa = 3.0d0*t1/(1.0d0+caz)
      bb = 14.7d0*t1/(28.0d0+caz)
      ak = (log(ak)+caz*cos(aa)/(1.0d0+0.008d0*caz))/cos(bb)
      fk = 0.12125d0*ak*ak/caz + 1.5d0
  180 continue
c-----------------------------------------------------------------------
c     backward recurrence loop for miller algorithm
c-----------------------------------------------------------------------
      k = fk
      fk = k
      fks = fk*fk
      p1r = czeror
      p1i = czeroi
      p2r = tol
      p2i = czeroi
      csr = p2r
      csi = p2i
      do 190 i=1,k
        a1 = fks - fk
        ak = (fks+fk)/(a1+fhs)
        rak = 2.0d0/(fk+coner)
        cbr = (fk+zr)*rak
        cbi = zi*rak
        ptr = p2r
        pti = p2i
        p2r = (ptr*cbr-pti*cbi-p1r)*ak
        p2i = (pti*cbr+ptr*cbi-p1i)*ak
        p1r = ptr
        p1i = pti
        csr = csr + p2r
        csi = csi + p2i
        fks = a1 - fk + coner
        fk = fk - coner
  190 continue
c-----------------------------------------------------------------------
c     compute (p2/cs)=(p2/abs(cs))*(conjg(cs)/abs(cs)) for better
c     scaling
c-----------------------------------------------------------------------
      tm = zabs(csr,csi)
      ptr = 1.0d0/tm
      s1r = p2r*ptr
      s1i = p2i*ptr
      csr = csr*ptr
      csi = -csi*ptr
      call zmlt(coefr, coefi, s1r, s1i, str, sti)
      call zmlt(str, sti, csr, csi, s1r, s1i)
      if (inu.gt.0 .or. n.gt.1) go to 200
      zdr = zr
      zdi = zi
      if(iflag.eq.1) go to 270
      go to 240
  200 continue
c-----------------------------------------------------------------------
c     compute p1/p2=(p1/abs(p2)*conjg(p2)/abs(p2) for scaling
c-----------------------------------------------------------------------
      tm = zabs(p2r,p2i)
      ptr = 1.0d0/tm
      p1r = p1r*ptr
      p1i = p1i*ptr
      p2r = p2r*ptr
      p2i = -p2i*ptr
      call zmlt(p1r, p1i, p2r, p2i, ptr, pti)
      str = dnu + 0.5d0 - ptr
      sti = -pti
      call zdiv(str, sti, zr, zi, str, sti)
      str = str + 1.0d0
      call zmlt(str, sti, s1r, s1i, s2r, s2i)
c-----------------------------------------------------------------------
c     forward recursion on the three term recursion with relation with
c     scaling near exponent extremes on kflag=1 or kflag=3
c-----------------------------------------------------------------------
  210 continue
      str = dnu + 1.0d0
      ckr = str*rzr
      cki = str*rzi
      if (n.eq.1) inu = inu - 1
      if (inu.gt.0) go to 220
      if (n.gt.1) go to 215
      s1r = s2r
      s1i = s2i
  215 continue
      zdr = zr
      zdi = zi
      if(iflag.eq.1) go to 270
      go to 240
  220 continue
      inub = 1
      if(iflag.eq.1) go to 261
  225 continue
      p1r = csrr(kflag)
      ascle = bry(kflag)
      do 230 i=inub,inu
        str = s2r
        sti = s2i
        s2r = ckr*str - cki*sti + s1r
        s2i = ckr*sti + cki*str + s1i
        s1r = str
        s1i = sti
        ckr = ckr + rzr
        cki = cki + rzi
        if (kflag.ge.3) go to 230
        p2r = s2r*p1r
        p2i = s2i*p1r
        str = abs(p2r)
        sti = abs(p2i)
        p2m = max(str,sti)
        if (p2m.le.ascle) go to 230
        kflag = kflag + 1
        ascle = bry(kflag)
        s1r = s1r*p1r
        s1i = s1i*p1r
        s2r = p2r
        s2i = p2i
        str = cssr(kflag)
        s1r = s1r*str
        s1i = s1i*str
        s2r = s2r*str
        s2i = s2i*str
        p1r = csrr(kflag)
  230 continue
      if (n.ne.1) go to 240
      s1r = s2r
      s1i = s2i
  240 continue
      str = csrr(kflag)
      yr(1) = s1r*str
      yi(1) = s1i*str
      if (n.eq.1) return
      yr(2) = s2r*str
      yi(2) = s2i*str
      if (n.eq.2) return
      kk = 2
  250 continue
      kk = kk + 1
      if (kk.gt.n) return
      p1r = csrr(kflag)
      ascle = bry(kflag)
      do 260 i=kk,n
        p2r = s2r
        p2i = s2i
        s2r = ckr*p2r - cki*p2i + s1r
        s2i = cki*p2r + ckr*p2i + s1i
        s1r = p2r
        s1i = p2i
        ckr = ckr + rzr
        cki = cki + rzi
        p2r = s2r*p1r
        p2i = s2i*p1r
        yr(i) = p2r
        yi(i) = p2i
        if (kflag.ge.3) go to 260
        str = abs(p2r)
        sti = abs(p2i)
        p2m = max(str,sti)
        if (p2m.le.ascle) go to 260
        kflag = kflag + 1
        ascle = bry(kflag)
        s1r = s1r*p1r
        s1i = s1i*p1r
        s2r = p2r
        s2i = p2i
        str = cssr(kflag)
        s1r = s1r*str
        s1i = s1i*str
        s2r = s2r*str
        s2i = s2i*str
        p1r = csrr(kflag)
  260 continue
      return
c-----------------------------------------------------------------------
c     iflag=1 cases, forward recurrence on scaled values on underflow
c-----------------------------------------------------------------------
  261 continue
      helim = 0.5d0*elim
      elm = exp(-elim)
      celmr = elm
      ascle = bry(1)
      zdr = zr
      zdi = zi
      ic = -1
      j = 2
      do 262 i=1,inu
        str = s2r
        sti = s2i
        s2r = str*ckr-sti*cki+s1r
        s2i = sti*ckr+str*cki+s1i
        s1r = str
        s1i = sti
        ckr = ckr+rzr
        cki = cki+rzi
        as = zabs(s2r,s2i)
        alas = log(as)
        p2r = -zdr+alas
        if(p2r.lt.(-elim)) go to 263
        call zlog(s2r,s2i,str,sti,idum)
        p2r = -zdr+str
        p2i = -zdi+sti
        p2m = exp(p2r)/tol
        p1r = p2m*cos(p2i)
        p1i = p2m*sin(p2i)
        call zuchk(p1r,p1i,nw,ascle,tol)
        if(nw.ne.0) go to 263
        j = 3 - j
        cyr(j) = p1r
        cyi(j) = p1i
        if(ic.eq.(i-1)) go to 264
        ic = i
        go to 262
  263   continue
        if(alas.lt.helim) go to 262
        zdr = zdr-elim
        s1r = s1r*celmr
        s1i = s1i*celmr
        s2r = s2r*celmr
        s2i = s2i*celmr
  262 continue
      if(n.ne.1) go to 270
      s1r = s2r
      s1i = s2i
      go to 270
  264 continue
      kflag = 1
      inub = i+1
      s2r = cyr(j)
      s2i = cyi(j)
      j = 3 - j
      s1r = cyr(j)
      s1i = cyi(j)
      if(inub.le.inu) go to 225
      if(n.ne.1) go to 240
      s1r = s2r
      s1i = s2i
      go to 240
  270 continue
      yr(1) = s1r
      yi(1) = s1i
      if(n.eq.1) go to 280
      yr(2) = s2r
      yi(2) = s2i
  280 continue
      ascle = bry(1)
      call zkscl(zdr,zdi,fnu,n,yr,yi,nz,rzr,rzi,ascle,tol,elim)
      inu = n - nz
      if (inu.le.0) return
      kk = nz + 1
      s1r = yr(kk)
      s1i = yi(kk)
      yr(kk) = s1r*csrr(1)
      yi(kk) = s1i*csrr(1)
      if (inu.eq.1) return
      kk = nz + 2
      s2r = yr(kk)
      s2i = yi(kk)
      yr(kk) = s2r*csrr(1)
      yi(kk) = s2i*csrr(1)
      if (inu.eq.2) return
      t2 = fnu + (kk-1)
      ckr = t2*rzr
      cki = t2*rzi
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
      s1r = coefr
      s1i = coefi
      s2r = coefr
      s2i = coefi
      go to 210
c
c
  310 continue
      nz=-2
      return
      end
